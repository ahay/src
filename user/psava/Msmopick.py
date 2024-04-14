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
def myMean(a,n):
    
    N = np.size(a)   # number of samples
    m = int((n-1)/2) # assume an odd window size
    
    b = np.zeros(N)
    for i in range(N):
        wLO = np.max([i-m,  0])  #  low index
        wHI = np.min([i+m,N-1])  # high index
        b[i] = np.mean( a[wLO:wHI] )
    
    return b

# ------------------------------------------------------------
def myMedian(a,n):
    
    N = np.size(a)   # number of samples
    m = int((n-1)/2) # assume an odd window size
    
    b = np.zeros(N)
    for i in range(N):
        wLO = np.max([i-m,  0])  #  low index
        wHI = np.min([i+m,N-1])  # high index  
        b[i] = np.median( a[wLO:wHI] )
    
    return b
# ------------------------------------------------------------

# ------------------------------------------------------------
par = rsf.Par()

verb = par.bool('verb',False) # verbosity flag
mode = par.int ('mode',0)     # smoothing mode

if mode == 0 or mode == 1:
    nwin = par.int('nwin',1) # window size (mean & median)

# ------------------------------------------------------------
Fin = rsf.Input()             # input file
n1 = Fin.int  ("n1")
o1 = Fin.float("o1")
d1 = Fin.float("d1")

n2 = Fin.int  ("n2")
o2 = Fin.float("o2")
d2 = Fin.float("d2")

nd = n2                       # number of samples
nm = nd

# ------------------------------------------------------------
Fou = rsf.Output()            # output file

Fou.put("n1",n2)
Fou.put("o1",o2)
Fou.put('d1',d2)

Fou.put("n2",3)
Fou.put("n3",1)
Fou.put("n4",1)

# ------------------------------------------------------------
# pick max value and weight
# ------------------------------------------------------------
din = np.zeros(n1,'f') # xcorrelation
pck = np.zeros(nd,'f') # picks
wgh = np.zeros(nd,'f') # weights

# ------------------------------------------------------------
l1 = 0.5 * n1*d1              # xcorr time window

for i in range(nd):
    Fin.read(din)

    pck[i] = o1 + np.argmax(din) * d1  # raw picks
    wgh[i] = np.max(din)               # data weight

wgh /= np.max(wgh)
pckmed = np.median(pck)
#pckmed = np.median(pck[pck < 0.9*l1 ])
#print(pckmed,file=sys.stderr)

wgh = np.where(np.abs(pck) < 0.9*l1, wgh, 0.0)
pck = np.where(np.abs(pck) < 0.9*l1, pck, 0.0)

if   mode == 0: # use mean
    print("USE MEAN",mode,file=sys.stderr)
    spk =   myMean(pck,nwin)

elif mode == 1: # use median
    print("USE MEDIAN",mode,file=sys.stderr)
    spk = myMedian(pck,nwin)

else:           # solve inverse problem
    print("USE INVERSION",mode,file=sys.stderr)

    #print(np.mean(pck),file=sys.stderr)
    #print(np.mean(wgh),file=sys.stderr)

    # ------------------------------------------------------------
    mbar = pck*0 + pckmed                 # reference model
    dbar = pck                            # rough picks

    WDop = spdiags(wgh,[0],nd,nd)         # data weight operator
    Gop  = idnop1D(nd)                    # mapping operator
    Rop  = lapop1D(nm)                    # regularization operator

    # ------------------------------------------------------------
    # L-curve
    # ------------------------------------------------------------
    ne = par.int  ('ne',1)
    oe = par.float('oe',0.0)
    de = par.float('de',+0.1)
    ee = np.arange(oe, oe+ne*de, de)
    se = np.power(10,ee)

    ME = np.zeros( (nm,ne), dtype='float')
    rd = np.zeros(     ne,  dtype='float')
    rm = np.zeros(     ne,  dtype='float')

    for ie in range(ne):

        # scale the regularization operator
        WMop =  Rop / se[ie] 

        # solve the IP
        modE = spsolve( (WDop*Gop).T * (WDop*Gop)  + WMop.T * WMop , \
                        (WDop*Gop).T * WDop * dbar + WMop.T * WMop * mbar)
    
        # store the model
        ME[:,ie] = modE 

        # compute residual norms
        rd[ie] = np.linalg.norm( WDop * (Gop * modE - dbar))
        rm[ie] = np.linalg.norm(  Rop * (      modE - mbar))

    rdn = rd / np.max(rd) # normalized the  data residual norm
    rmn = rm / np.max(rm) # normalized the model residual norm
    rrn = rdn**2 + rmn**2 # compute the distance from origin

    je = np.argmin(rrn)   # index of the optimal model
    spk = ME[:,je]        #              optimal model


# ------------------------------------------------------------
# write picks and weights
# ------------------------------------------------------------
Fou.write(spk) # smooth picks
Fou.write(pck) # rough picks
Fou.write(wgh) # xcor pick weight

Fin.close()
Fou.close()
