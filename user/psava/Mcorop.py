#!/usr/bin/env python
'''
    Correlation operator w/ adjoint
    wfl [file] : is taken from stdin
    opr [file] : is taken from  "opr"
    Requires both files to have the same dimensions
'''
import rsf.api as rsf
import numpy as np
import sys

par = rsf.Par()

verb = par.bool('verb',False) # verbosity flag
adj  = par.bool( 'adj',False) # adjoint flag

# Operator
Fopr = rsf.Input("opr")
nt = Fopr.int  ("n1")
ot = Fopr.float("o1")
dt = Fopr.float("d1")
lt = Fopr.string("label1")
ut = Fopr.string("unit1")

nn = Fopr.size(1) # number of traces

opr = np.zeros(nt,'f') # allocate opr array

# ------------------------------------------------------------
# setup output
if adj==1:
    Fcor = rsf.Input()         # input cor

    nc = Fcor.int("n1")
    ncor = (nc-1)/2

    Fwfl = rsf.Output()        # output wfl
    Fwfl.put("n1",nt)
    Fwfl.put("o1",ot)
    Fwfl.put('d1',dt)
    Fwfl.put('label1',lt)
    Fwfl.put('unit1',ut)

else:
    Fwfl = rsf.Input()         #  input wfl

    ncor = par.int("ncor",100)
    nc=2*ncor+1 # number of correlation lags

    Fcor = rsf.Output()        # output cor
    Fcor.put("n1",2*ncor+1)
    Fcor.put("o1", -ncor*dt)
    Fcor.put('d1',       dt)

# ------------------------------------------------------------
if adj==1:
    print  >> sys.stderr,"ADJT op"
    cor = np.zeros(2*nt-1,'f') # allocate "full" cor array
    for i in range(nn):
        Fopr.read(opr)
        Fcor.read(cor[nt-ncor:nt+ncor+1]) # insert cor in "full" array
        wfl = np.convolve(opr,cor,mode='valid')
        Fwfl.write(wfl)
else:
    print  >> sys.stderr,"FORW op"
    wfl = np.zeros(nt,'f')
    for i in range(nn):
        Fopr.read(opr)
        Fwfl.read(wfl)
        cor = np.correlate(wfl,opr,mode='full')
        Fcor.write(cor[nt-ncor:nt+ncor+1])

# ------------------------------------------------------------
Fopr.close()
Fwfl.close()
Fcor.close()
