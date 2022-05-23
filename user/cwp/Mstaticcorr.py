#!/usr/bin/env python
'''
Static Corrections for Redatuming SHARAD data to MOLA surface
Input: 2D correlation between SHARAD data and clutter simulation
Output: Smoothed maximum correlation picks representing residual ionosphere time delay
Created on: Mar 9, 2021
'''
import rsf.api as rsf
import numpy as np
import sys
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve,inv
from scipy.signal import medfilt
from scipy import stats


par = rsf.Par()

verb = par.bool('verb',False) # verbosity flag

# Local Functions
def idnop1D(n):
    # Creates a square matrix of size n x n with ones along the main diagonal
    e = np.ones( n, dtype='float')
    Lopr = spdiags([e], [0], n,n)
    return Lopr

def lapop1D(n):
    # Creates a square matrix of size n x n with:
    # -2 on the main diagonal
    # 1 on the 1st upper diagonal
    # 1 on the 1st lower diagonal
    e = np.ones( n, dtype='float')
    Lopr = spdiags([e, -2*e, e], [-1, 0, +1], n,n)
    return Lopr

# Read standard input
Fin = rsf.Input()             # input 2D correlation file
nt = Fin.int("n1")            # number of time samples
dt = Fin.float("d1")
ot = Fin.float("o1")
maxt = ot + (dt*nt)
tax = np.arange(ot,maxt,dt) # time axis

nr = Fin.int('n2') # number of records in orbit
cor_input = np.zeros((nt,nr),dtype=float)

# Read data to numpy
j2 = 1
dn1 = np.zeros(nt,'f')
i=0
while i < nr:
	Fin.read(dn1)
	cor_input[:,i] = dn1[:]
	i+=j2

# Find maximum correlation for each record
picks = np.zeros(nr)
for i in range(nr):
	picks[i] = tax[np.argmax(cor_input[:,i])]

# Define model and data parameters
xx = np.arange(1,nr+1,1)
nm = nr
nd = nr

Gop  = idnop1D(nd)	# mapping (data domain to model domain) operator
Rop  = lapop1D(nm)	# Laplacian operator

# ------------------------------------------------------------
# model domain
xmin = np.min(xx)
xmax = np.max(xx)

ox = xmin
nx = nm
dx = 1 #np.abs(xmax-xmin)/(nm-1), hardcoding this as one for speed
ax = (nx,ox,dx, xmax)

# ------------------------------------------------------------
# Index of first 2% of samples
ind_2pc = int(round(nx*0.02))
tS = np.median(picks[0:ind_2pc])
xS = xx[0]

tE = np.median(picks[nd-1-ind_2pc:nd-1])
xE = xx[nd-1]

# reference model
#mbar = tS + (xx-xS)/(xE-xS) * (tE-tS)
picks_mode = stats.mode(picks)[0][0]
mbar = np.ones(nr)*picks_mode

# ------------------------------------------------------------
# data domain
dbar = picks

# data uncertainty
sigd = 1e-1
#dsig = sigd + 5e5 * np.abs(Rop * tin)

# Number of samples for 5% and 3% kernel
ker_5pc = int(round(nx*0.05))
if ker_5pc%2==0:
	ker_5pc+=1
ker_3pc = int(round(nx*0.03))
if ker_3pc%2==0:
	ker_3pc+=1

tmd = medfilt(picks,ker_5pc)
#dsig = sigd + 1e3 * np.abs(picks - medfilt(picks,ker_3pc))
dsig = sigd + 1e5 * np.abs(picks - picks_mode) + 1e4 * np.abs(picks - medfilt(picks,ker_3pc))

WDop = spdiags(1/dsig,[0],nd,nd) #  data weight

# ------------------------------------------------------------
# Model fitting
ne = 61                             # number of ?
ee = np.linspace(-2.0, +2.0, ne)
se = np.power(10,ee)

ME = np.zeros( (nx,ne), dtype='float')
rd = np.zeros(     ne,  dtype='float')
rm = np.zeros(     ne,  dtype='float')

for ie in range(ne):
	# setup regularization
	WMop =  Rop / se[ie]

	# solve IP
	modE = spsolve( (WDop*Gop).T * (WDop*Gop) + WMop.T * WMop , (WDop*Gop).T * WDop * dbar + WMop.T * WMop * mbar)

	# store model
	ME[:,ie] = modE

	# compute residual norms
	rd[ie] = np.linalg.norm( WDop * (Gop * modE - dbar))
	rm[ie] = np.linalg.norm(  Rop * (      modE - mbar))

# L-curve optimization
# normalized residual norm
rdn = rd / np.max(rd)
rmn = rm / np.max(rm)

# distance from origin in residual space
rrn = rdn**2 + rmn**2

# index of optimal model
je = np.argmin(rrn)

# optimal curve
tou = ME[:,je]

# Write standard output
Fou = rsf.Output()
Fou.put('n1',nr)
Fou.put('o1',0)
Fou.put('d1',1)
Fou.put('n2',1)
Fou.put('d2',1)
Fou.put('o2',0)
Fou.write(tou)

#------------------------------------------------
Fin.close()
Fou.close()
