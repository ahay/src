#!/usr/bin/env python

import sys, re, os
from numpy import *
from scipy.interpolate import interp1d
from scipy.interpolate import Rbf

# Dump RSF text header to standard output
def create_rsf_header (n, d, o, file):
    for i in range (size (n)):
        sys.stdout.write ('n%d=%d\n' % (i + 1, n[i]))
        sys.stdout.write ('d%d=%g\n' % (i + 1, d[i]))
        sys.stdout.write ('o%d=%g\n' % (i + 1, o[i]))
    sys.stdout.write ('data_format=\"native_float\"\n')
    sys.stdout.write ('esize=4\n')
    sys.stdout.write ('in=%s\n' % file)
    sys.stdout.flush ()

# Interpolation along time axis
def interp_from_t (t, v, tint):
    return interp1d (t, v, kind='linear')(tint).astype (float32)

# Interpolation of velocities in a constant-t slice
def interp_for_const_t (pnt, v, gx, gy):
    return Rbf(pnt[:,0], pnt[:,1], v, function='thin_plate')(gx, gy).T.astype (float32)

# Visualize one slice after interpolation
def show_slice (slice, pnt, nx, ny, title='Constant time slice'):
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    plt.imshow (slice, interpolation='nearest', cmap=cm.jet,
                extent = (1,nx,1,ny), origin='lower')
    plt.title (title)
    plt.plot (pnt[:,0], pnt[:,1], 'k.')
    plt.colorbar ()
    plt.show ()

#
# Survey parameters
nxline = 188
ninline = 345
finline = 1 # First inline
fxline = 1 # First xline
# Time axis
maxt = 3001
dt = 2
# Densely sampled time axis
tint = arange (0, maxt, dt, dtype = float32)
nt = size (tint)
# Interpolated velocity values on time axis
vint = arange (0, dtype = float32)
vint = vint.reshape (nt, 0)

#
# Regular expressions for a number in text form
rxcdp = re.compile ('\d{1,5}')
rxtvel = re.compile ('\d{1,5}\.\d{2,2}')

#
# Arrays for storing found CDP locations and time-velocity pairs
ts = []
vs = []
ils = []
xls = []
# Current CDP location in the search below
cdp = 1
inline = 1
xline = 1

#
# Execution starts here
#

# Go over text lines from standard input and
# find time-velocities pairs for each CDP
for line in sys.stdin.readlines ():
    if line.find ('CDP=') != -1:
        # Store previous CDP location
        if len (ts) > 0:
           vint = append (vint, interp_from_t (ts, vs, tint).reshape (nt, 1),
                          axis = 1)
        # Remember current CDP location
        scdp = rxcdp.search (line)
        cdp = int(scdp.group ())
        # Inline
        ils.append (cdp / nxline + finline)
        # Xline
        xls.append (cdp % nxline + fxline)
        ts = []
        vs = []
    elif line.find ('TIME=') != -1 and \
         line.find ('VEL=') != -1:
        # Get time and velocity
        tvel = rxtvel.findall (line)
        ts.append (float(tvel[0]))
        vs.append (float(tvel[1]))
# Store last CDP location
vint = append (vint, interp_from_t (ts, vs, tint).reshape (nt, 1),
               axis = 1)

#
# End of standard input scanning, proceed to interpolation
#

# Read survey mask
mask = ones ((ninline, nxline))
maskfile = 'surveymask.dat'
if (len(sys.argv) > 1):
    maskfile = sys.argv[1]
if (os.path.exists (maskfile)):
    mask = fromfile (maskfile, dtype = float32).reshape (ninline, nxline)

# Number of points per time slice
np = len (ils)

# Inline/xline grid for a time slice
gx, gy = mgrid[1:nxline+1:1, 1:ninline+1:1]
# Array of x,y points on the time slice
pnt = array (xls, dtype = float32).reshape (np, 1)
pnt = append (pnt, array (ils, dtype = float32).reshape (np, 1), axis = 1)

# Get one slice and visualize for test
#slice = mask*interp_for_const_t (pnt, vint[500,:], gx, gy)
#show_slice (slice, pnt, nxline, ninline, title="t=1.0s")

#
# Veloicty volume
vvol = zeros ((nt, ninline, nxline), dtype = float32)

# Gor along time slices and do triangulation
for it in range (nt):
    print >> sys.stderr, "Processing slice", it + 1, "of", nt
    vvol[it] = mask*interp_for_const_t (pnt, vint[it,:], gx, gy)

# Finally, dump trace by trace to binary file
binfile = 'vvol.dat' # Binary part
bfid = open (binfile, 'w+')
print vvol.shape
for il in range (ninline):
    print >> sys.stderr, "Saving inline", il + 1, "of", ninline
    for xl in range (nxline):
        vvol[:,il,xl].tofile (bfid)
bfid.flush ()
bfid.close ()

# Madagascar header
create_rsf_header (array ([nt, nxline, ninline], dtype = int32),
                   array ([dt, 1.0, 1.0], dtype = float32),
                   array ([0.0, fxline, finline], dtype = float32),
                   binfile)


