#!/usr/bin/env python

import sys, re
from numpy import *

#
# Survey parameters
nxline = 188
ninline = 345

#
# Regular expressions for a number in text form
rxixl = re.compile ('\d{1,3}\,')
rxinl = re.compile ('\d{1,3}i')

#
# Array for marking existing CDP bins
mask = zeros ((ninline,nxline), dtype = float32)

#
# Execution starts here
#

# Go over text lines from standard input and
# read coordinate pairs
for line in sys.stdin.readlines ():
    ixline = int(rxixl.search (line).group ()[:-1])
    inline = int(rxinl.search (line).group ()[:-1])
    mask[inline - 1, ixline - 1] = 1.0
#
# End of standard input scanning, proceed to interpolation
#

#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#plt.imshow (mask, interpolation='nearest', cmap=cm.jet, origin='lower')
#plt.colorbar ()
#plt.show ()

mask.tofile (sys.stdout)

