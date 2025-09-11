import m8r, numpy, sys
import matplotlib.pyplot as plt
from rsfpy import grey3

if __name__ == "__main__":
    # check stdin
    if sys.stdin.isatty():
        print("Usage: python greytest.py < data.rsf [> output.pdf]", file=sys.stderr)
        sys.exit(1)
    # read header from stdin
    inp = m8r.Input()
    n1, n2, n3 = inp.int('n1'), inp.int('n2'),  inp.int('n3')
    d1, d2, d3 = inp.float('d1'), inp.float('d2'), inp.float('d3')
    o1, o2, o3 = inp.float('o1'), inp.float('o2'), inp.float('o3')
    l1, l2, l3 = inp.string('label1'), inp.string('label2'), inp.string('label3')
    u1, u2, u3 = inp.string('unit1'), inp.string('unit2'), inp.string('unit3')
    title = inp.string('title')

    # Read parameters from command line
    pars = m8r.Par()
    frame1, frame2, frame3 = pars.int('frame1', 0), pars.int('frame2', 0), pars.int('frame3', 0)
    point1, point2 = pars.float('point1', 0.8), pars.float('point2', 0.4)
    clip = pars.float('clip', None)
    pclip = pars.float('pclip', 99.)
    isflat = pars.bool('flat', True)
    allpos = pars.bool('allpos', False)
    colorbar = pars.bool('scalebar', False)
    color = pars.string('color', 'gray')
    title = pars.string('title', title)

    # Validate parameters
    if frame1 < 0: frame1 = 0
    if frame2 < 0: frame2 = 0
    if frame3 < 0: frame3 = 0
    if frame1 >= n1: frame1 = n1 - 1
    if frame2 >= n2: frame2 = n2 - 1
    if frame3 >= n3: frame3 = n3 - 1

    if point1 < 0 or point1 > 1: point1 = 0.8
    if point2 < 0 or point2 > 1: point2 = 0.4
    if pclip < 0 or pclip > 100: pclip = 99.

    # Read data 
    data = inp.read(datatype=numpy.float32)
    data = data.transpose((2, 1, 0))  # Transpose to match RSF's (n1, n3, n2) order

    fig = plt.figure(figsize=(8, 8))
    # display data
    grey3(data, flat=isflat, frame1=frame1, frame2=frame2, frame3=frame3,
          point1=point1, point2=point2, clip=clip, pclip=pclip,
          title=title, label1=f'{l1} ({u1})',
          label2=f'{l2} ({u2})', label3=f'{l3} ({u3})',
          d1=d1, d2=d2, d3=d3,
          o1=o1, o2=o2, o3=o3,
          color=color, allpos=allpos, colorbar=colorbar)
    # check stdout
    if sys.stdout.isatty(): fig.show()
    else: fig.savefig(sys.stdout.buffer, format='pdf', dpi=600)
    sys.exit(0)
    