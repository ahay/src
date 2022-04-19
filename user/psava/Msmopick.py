#!/usr/bin/env python
'''
smooth picker on the first dimension
'''
import rsf.api as rsf
import numpy as np
from scipy.sparse        import spdiags
from scipy.sparse.linalg import spsolve
import sys

#import idnUtil as IDN
#import lapUtil as LAP

# ------------------------------------------------------------
def idnop1D(n):
    e = np.ones( n, dtype='float')
    Lopr = spdiags([e], [0], n,n)
    return Lopr

def idnop2D(nx,nz):
    n = nx*nz
    e = np.ones( n, dtype='float')

    Lopr = spdiags([e], [0], n,n)
    return Lopr
# ------------------------------------------------------------
def lapop1D(n):
    e = np.ones( n, dtype='float')
    Lopr = spdiags([e, -2*e, e], [-1, 0, +1], n,n)
    return Lopr

def lapop2D(nx,nz):
    n = nx*nz
    e = np.ones( n, dtype='float')

    u = e.copy()
    u[::nz]=0

    l = e.copy()
    l[nz-1::nz]=0

    Lopr = spdiags([e, l, -4*e, u,e], [-nz, -1, 0, +1, +nz], n,n)
    return Lopr
# ------------------------------------------------------------

par = rsf.Par()

verb = par.bool('verb',False) # verbosity flag

# ------------------------------------------------------------
Fin = rsf.Input()             # input file
n1 = Fin.int  ("n1")
o1 = Fin.float("o1")
d1 = Fin.float("d1")

nd = Fin.size(1);             # number of traces
nm = nd

# ------------------------------------------------------------
Fou = rsf.Output()            # output file
Fou.put("n1",nd)
Fou.put("o1",0.0)
Fou.put('d1',1.0)

Fou.put("n2",3) # temporary

Fou.put("n3",1)
Fou.put("n4",1)

# ------------------------------------------------------------
# pick max value and weight
# ------------------------------------------------------------

din = np.zeros(n1,'f')
pck = np.zeros(nd,'f') # picks
wgh = np.zeros(nd,'f') # weights

for i in range(nd):
    Fin.read(din)
    pck[i] = o1 + np.argmax(din) * d1
    wgh[i] = np.max(din)

wgh /= np.max(wgh)

# ------------------------------------------------------------
# smooth picks
# ------------------------------------------------------------

xx = np.arange(nd)
xS = xx[0]
xE = xx[nd-1]

nb = np.max([1,nd//10])
#tS = np.median(pck[      0:  nb])     # first 1% of picks
#tE = np.median(pck[nd-1-nb:nd-1])     #  last 1% of picks
#mbar = tS + (xx-xS)/(xE-xS) * (tE-tS) #   reference model
#mbar = pck*0 + np.median(pck)
mbar = pck*0

dbar = pck                            # rough picks

wpo = par.float('wpo',1.0)
wgh = np.power(wgh-np.min(wgh),wpo)
WDop = spdiags(wgh,[0],nd,nd)         # data weight operator
Gop  = idnop1D(nd)                # mapping operator
Rop  = lapop1D(nm)                # Laplacian operator

# ------------------------------------------------------------
# L-curve
ne = par.int('ne',1)
oe = par.float('oe',0.0)
de = par.float('de',+0.1)
ee = np.arange(oe, oe+ne*de, de)
#print(ee,file=sys.stderr)
se = np.power(10,ee)

ME = np.zeros( (nm,ne), dtype='float')
rd = np.zeros(     ne,  dtype='float')
rm = np.zeros(     ne,  dtype='float')

for ie in range(ne):
    WMop =  Rop / se[ie] # setup regularization

    # solve IP
    modE = spsolve( (WDop*Gop).T * (WDop*Gop)  + WMop.T * WMop , \
                    (WDop*Gop).T * WDop * dbar + WMop.T * WMop * mbar)
    ME[:,ie] = modE # store model

    # compute residual norms
    rd[ie] = np.linalg.norm( WDop * (Gop * modE - dbar))
    rm[ie] = np.linalg.norm(  Rop * (      modE - mbar))

rdn = rd / np.max(rd) # normalized residual norm
rmn = rm / np.max(rm)
rrn = rdn**2 + rmn**2 # distance from origin
je = np.argmin(rrn)   # index of the optimal model
spk = ME[:,je]        # smooth picks

# ------------------------------------------------------------
# write picks
# ------------------------------------------------------------
Fou.write(pck) # rough picks
Fou.write(spk) # smooth picks
Fou.write(wgh)

Fin.close()
Fou.close()
