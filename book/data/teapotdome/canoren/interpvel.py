#!/usr/bin/env python

from __future__ import print_function
import sys, re, os
from numpy import *
from scipy.interpolate import interp1d
from scipy.interpolate import Rbf
import rsf.api as rsf

# This Python code is a modified version of what already exists in the Madagascar repository "https://github.com/ahay/src/tree/master/book/data/teapotdome/vels" and is provided with the manuscript entitled "An Overview of Reproducible 3D Seismic Data Processing and Imaging Using Madagascar" published in Geophysics 83, 2 (2018), pp. F9-F20 by Can Oren and Robert L. Nowack. 

######################
# This SEG version of the code may be found at:
# http://software.seg.org/2018/0002
# It is governed by the SEG copyright, which may be found
# in the file "disclaimer.txt" or at
# http://software.seg.org/disclaimer.txt
# Note the version in the Madagascar repository has
# its own separate copyright statement.
# The SEG copyright only applies to the version of the code
# downloaded from the SEG website.
######################

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

# Get parameters from command line
command_line_par=rsf.Par()

# Survey parameters
nxline=command_line_par.int("nxline")
if nxline==None:
    print("nxline is a required parameter in interpvel")
    sys.exit(2)

ninline=command_line_par.int("ninline")
if ninline==None:
    print("ninline is a required parameter in interpvel")
    sys.exit(2)

fxline=command_line_par.int("fxline")
if fxline==None:
    print("fxline is a required parameter in interpvel")
    sys.exit(2)

finline=command_line_par.int("finline")
if finline==None:
    print("finline is a required parameter in interpvel")
    sys.exit(2)

# Time axis
maxt=command_line_par.float("maxt")
if maxt==None:
    print("maxt is a required parameter in interpvel")
    sys.exit(2)

dt=command_line_par.float("dt")
if dt==None:
    print("dt is a required parameter in interpvel")
    sys.exit(2)

print("nxline=",nxline," fxline=",fxline,file=sys.stderr)
print("ninline=",ninline," finline=",finline,file=sys.stderr)
print("maxt=",maxt," dt=",dt,file=sys.stderr)

maskfile =command_line_par.string("mask")
print("maskfile=",maskfile,file=sys.stderr)
# Densely sampled time axis
tint = arange (0, maxt, dt, dtype = float32)
nt = size (tint)
# Interpolated velocity values on time axis
vint = arange (0, dtype = float32)
vint = vint.reshape (nt, 0)

# Regular expressions for a number in text form
rxcdp = re.compile ('\d{1,5}')
rxtvel = re.compile ('\d{1,6}\.\d{2,3}')

# Arrays for storing found CDP locations and time-velocity pairs
ts = []
vs = []
ils = []
xls = []

# Current CDP location in the search below
cdp = 1
inline = 1
xline = 1

# Execution starts here #

# Go over text lines from standard input and find time-velocity pairs for each CDP
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

# End of standard input scanning, proceed to interpolation

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

# Velocity volume
vrms3d = zeros ((nt, ninline, nxline), dtype = float32)

# Go along time slices and do triangulation
for it in range (nt):
    if it<3 or it%50==0:
        print("Processing slice", it + 1, "of", nt,file=sys.stderr)
    if it==2:
        print("now print status for every 50 slices",file=sys.stderr)
    vrms3d[it] = mask*interp_for_const_t (pnt, vint[it,:], gx, gy)
print("Completed all time slices",file=sys.stderr)

# Finally, dump trace by trace to binary file
binfile = 'vrms3d.dat' # Binary part
bfid = open (binfile, 'w+')
print(vrms3d.shape)
for il in range (ninline):
    if il<3 or il%50==0:
        print("Saving inline", il + 1, "of", ninline,file=sys.stderr)
    if il==2:
        print("now print status for every 50 lines",file=sys.stderr)
    for xl in range (nxline):
        vrms3d[:,il,xl].tofile (bfid)
bfid.flush ()
bfid.close ()

# Madagascar header
create_rsf_header (array ([nt, nxline, ninline], dtype = int32),
                   array ([dt, 1.0, 1.0], dtype = float32),
                   array ([0.0, fxline, finline], dtype = float32),
                   binfile)
