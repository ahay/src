#!/usr/bin/env python
'''
smooth picker on the first dimension
'''
import rsf.api as rsf
import numpy as np
from scipy.sparse        import spdiags
from scipy.sparse.linalg import spsolve
import sys

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
wpo = par.float('wpo',1.0)    # weight exponent

# ------------------------------------------------------------
Fin = rsf.Input()             # input file
n1 = Fin.int  ("n1")
o1 = Fin.float("o1")
d1 = Fin.float("d1")
l1 = 0.5 * n1*d1              # xcorr time window
z1 = int( abs(o1)/d1)

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

    # raw picks
    #pck[i] = o1 + np.argmax(din) * d1
    #wgh[i] = np.max(din)

    din += 1.0 # avoid zero weights
    pck[i] = o1 + z1 * d1
    wgh[i] = din[z1]
    for j in range(n1):
        if  wgh[i] < din[j]:
            pck[i] = o1 + j * d1
            wgh[i] = din[j]

# bias weights
pckmed = np.median(pck)
for i in range(nd):
    #wgh[i]*= np.exp(- abs(pck[i]) )
    wgh[i] *= np.cos( 0.5*np.pi * abs((pck[i]-pckmed)/l1) )

wgh /= np.max(wgh)

# ------------------------------------------------------------
# smooth picks
# ------------------------------------------------------------

xx = np.arange(nd)
xS = xx[0]
xE = xx[nd-1]

nb = np.max([1,nd//10])
#mbar = pck*0                          # reference model
mbar = pck*0 + pckmed
dbar = pck                            # rough picks

wgh = np.power(wgh-np.min(wgh),wpo)

WDop = spdiags(wgh,[0],nd,nd)         # data weight operator
Gop  = idnop1D(nd)                    # mapping operator
Rop  = lapop1D(nm)                    # Laplacian operator

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

if rd > 0:
    rdn = rd / np.max(rd) # normalized residual norm
else:
    rdn = 1.0
if rm > 0:
    rmn = rm / np.max(rm)
else:
    rmn = 1.0
rrn = rdn**2 + rmn**2 # distance from origin
je = np.argmin(rrn)   # index of the optimal model
spk = ME[:,je]        # smooth picks

# ------------------------------------------------------------
# write picks
# ------------------------------------------------------------
Fou.write(spk) # smooth picks
Fou.write(pck) # rough picks
Fou.write(wgh)

Fin.close()
Fou.close()
