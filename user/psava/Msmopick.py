#!/usr/bin/env python
'''
smooth picker on the first dimension
'''
import rsf.api as rsf
import numpy as np
import sys

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


par = rsf.Par()

verb = par.bool('verb',False) # verbosity flag
mean = par.bool('mean',False) # verbosity flag
nmed = par.int('ns',1)        # number of filter samples

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
for i in range(nd):
    Fin.read(din)

    pck[i] = o1 + np.argmax(din) * d1  # raw picks
    wgh[i] = np.max(din)               # data weight

#spk = signal.medfilt(pck,2751)
if mean:
    spk = myMean(pck,nmed)
else:
    spk = myMedian(pck,nmed)

# ------------------------------------------------------------
# write picks and weights
# ------------------------------------------------------------
Fou.write(spk) # smooth picks
Fou.write(pck) # rough picks
Fou.write(wgh) # xcor pick weight

Fin.close()
Fou.close()
