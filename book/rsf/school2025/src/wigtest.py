import m8r, numpy, sys
import matplotlib.pyplot as plt

def wiggle_plot(data, time, o2, d2, scale_factor=1.0, 
                color='black', linewidth=0.5, fill_positive=True):
    """
    Draws a wiggle plot for 2D array `data` of shape (n_samples, n_traces).
    time : 1D array of length n_samples
    o2, d2 : origin and sampling interval on the horizontal axis
    scale_factor : multiplies amplitudes before plotting
    """
    n_samples, n_traces = data.shape
    for itr in range(n_traces):
        baseline = o2 + itr * d2
        trace = data[:, itr] * scale_factor
        x = baseline + trace
        plt.plot(x, time, color=color, linewidth=linewidth)
        if fill_positive:
            plt.fill_betweenx(time, baseline, x, where=(trace > 0), 
                              facecolor=color, linewidth=0)

if __name__ == "__main__":
    # require stdin
    if sys.stdin.isatty():
        print("Usage: python wiggle.py < data.rsf [> output.pdf]", file=sys.stderr)
        sys.exit(1)

    inp = m8r.Input()
    # read RSF headers
    n1, n2    = inp.int('n1'), inp.int('n2')
    d1, d2    = inp.float('d1'), inp.float('d2')
    o1, o2    = inp.float('o1'), inp.float('o2')
    l1, l2    = inp.string('label1'), inp.string('label2')
    u1, u2    = inp.string('unit1'), inp.string('unit2')
    title     = inp.string('title')
    data      = inp.read(datatype=numpy.float32).T  # shape → (n1, n2)

    # build time (or depth) axis
    time = o1 + numpy.arange(n1) * d1

    # autoscale wiggles to occupy ~80% of trace spacing
    max_amp = numpy.max(numpy.abs(data))
    if max_amp > 0:
        scale = 0.8 * d2 / max_amp
    else:
        scale = 1.0

    # draw
    plt.figure(figsize=(8, 6))
    wiggle_plot(data, time, o2, d2, scale_factor=scale)

    # decorate
    plt.gca().invert_yaxis()  # time/depth grows downwards
    plt.xlabel(f"{l2} ({u2})")
    plt.ylabel(f"{l1} ({u1})")
    if title:
        plt.title(title)

    # output
    if sys.stdout.isatty():
        plt.show()
    else:
        plt.savefig(sys.stdout.buffer, format='pdf')
    sys.exit(0)
