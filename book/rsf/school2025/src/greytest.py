import m8r, numpy, sys
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # check stdin
    if sys.stdin.isatty():
        print("Usage: python greytest.py < data.rsf [> output.pdf]", file=sys.stderr)
        sys.exit(1)
    # read data from stdin
    inp = m8r.Input()
    # read data and axes
    n1, n2 = inp.int('n1'), inp.int('n2')
    d1, d2 = inp.float('d1'), inp.float('d2')
    o1, o2 = inp.float('o1'), inp.float('o2')
    l1, l2 = inp.string('label1'), inp.string('label2')
    u1, u2 = inp.string('unit1'), inp.string('unit2')
    title = inp.string('title')
    data = inp.read(datatype=numpy.float32).T

    # calculate clip
    clip = numpy.percentile(numpy.abs(data), 99)
    # display data
    plt.imshow(data, aspect='auto', cmap='gray',
               vmin=-clip, vmax=clip,
               extent=(o2, o2 + d2 * n2, o1 + d1 * n1, o1))
    plt.ylabel(f"{l1} ({u1})")
    plt.xlabel(f"{l2} ({u2})")
    if title: plt.title(title)
    # check stdout
    if sys.stdout.isatty(): plt.show()
    else: plt.savefig(sys.stdout.buffer, format='pdf', dpi=300)
    sys.exit(0)
    