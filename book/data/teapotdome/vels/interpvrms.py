#!/usr/bin/python

import sys, re
from numpy import *
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from scipy.interpolate import Rbf
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Interpolation along time axis
def interp_from_t (t, v):
    return interp1d (t, v, kind='linear')

# Interpolation of velocities in a constant-t slice
def interp_for_const_t (p, v, gx, gy):
    # Do triangulation and mark points outside of the convex hull with NaNs
    trslice = griddata (p, v, (gx, gy), method='linear')
    # Do radial-basis function interpolation
    slice = Rbf(pnt[:,0], pnt[:,1], v, function='thin_plate')(gx, gy)
    # Mark points outside of the convex hull with NaNs
    slice[where (isnan (trslice) == True)] = nan
    return slice
#    return griddata (p, v, (gx, gy), method='linear', fill_value = -1.0)

# Survey parameters
nxline = 188
ninline = 345
# Time axis
maxt = 3001
dt = 2
# Densely sampled time axis
tint = arange (0, maxt, dt)
nt = size (tint)

# Regular expressions for number
rxcdp = re.compile ('\d{1,5}')
rxtvel = re.compile ('\d{1,5}\.\d{2,2}')

# Arrays for storing found CDP locations and time-velocity pairs
ts = []
vs = []
ils = []
xls = []
# Current CDP location in the search below
cdp = 1
inline = 1
xline = 1

# Interpolated velocity values on time axis
vint = arange (0)
vint = vint.reshape (nt, 0)

#
# Execution starts here
#

# Go over text lines from standard input and find time-velocities pairs for each CDP
for line in sys.stdin.readlines ():
    if line.find ('CDP=') != -1:
        # Store previous CDP location
        if len (ts) > 0:
           vint = append (vint, interp_from_t (ts, vs)(tint).reshape (nt, 1),
                          axis=1)
        # Remember current CDP location
        scdp = rxcdp.search (line)
        cdp = int(scdp.group ())
        # Inline
        ils.append (cdp / nxline)
        # Xline
        xls.append (cdp % nxline)
        ts = []
        vs = []
    elif line.find ('TIME=') != -1 and \
         line.find ('VEL=') != -1:
        # Get time and velocity
        tvel = rxtvel.findall (line)
        ts.append (float(tvel[0]))
        vs.append (float(tvel[1]))
# Store last CDP location
vint = append (vint, interp_from_t (ts, vs)(tint).reshape (nt, 1),
               axis=1)

#
# End of standard input scanning, proceed to interpolation
#

# Number of points per time slice
np = len (ils)

# Inline/xline grid for a time slice
gx, gy = mgrid[1:nxline+1:1, 1:ninline+1:1]

# Array of x,y points on the time slice
pnt = array (xls).reshape (np, 1)
pnt = append (pnt, array (ils).reshape (np, 1), axis=1)

# Get one slice and visualize for test
slice = interp_for_const_t (pnt, vint[500,:], gx, gy)
plt.imshow (slice.T, interpolation='nearest', cmap=cm.jet,
            extent = (1,nxline,1,ninline), origin='lower')
plt.title ('t = 1.0s, Thin-Plate Spline interpolation')
#plt.title ('t = 1.0s, Cubic')
#plt.title ('t = 1.0s, Delaunay triangulation + linear interpolation')
plt.plot (pnt[:,0], pnt[:,1], 'k.')
plt.colorbar ()
plt.show ()

# Gor along time slices and do triangulation
#for it in range (nt):

