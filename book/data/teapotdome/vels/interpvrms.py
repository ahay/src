#!/usr/bin/python

import sys, re
from numpy import *
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Interpolation along time axis
def interp_from_t (t, v):
    return interp1d (t, v, kind='linear')
def interp_for_const_t (p, v, gx, gy):
    return griddata (p, v, (gx, gy), method='nearest')

# Survey parameters
nxline = 188
ninline = 345
# Time axis
maxt = 3001
dt = 2
tint = arange (0, maxt, dt)
nt = size (tint)

# Regular expressions for number
rxcdp = re.compile ('\d{1,5}')
rxtvel = re.compile ('\d{1,5}\.\d{2,2}')

cdp = 1
inline = 1
xline = 1

ts = []
vs = []
ils = []
xls = []

# Values for interpolation on time axis
vint = arange (0)
vint = vint.reshape (nt, 0)

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

#print vint[0,:]
# Number of points per time slice
np = len (ils)

# Inline/xline grid for a time slice
gx, gy = mgrid[1:nxline+1:1, 1:ninline+1:1]

# Array of x,y points on the time slice
pnt = array (xls).reshape (np, 1)
pnt = append (pnt, array (ils).reshape (np, 1), axis=1)
slice = interp_for_const_t (pnt, vint[1000,:], gx, gy)

plt.imshow (slice, interpolation='bilinear', cmap=cm.gray)
plt.show ()
# Gor along time slices and do triangulation
#for it in range (nt):

