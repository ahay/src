#!/usr/bin/env python
'''
convert amplitude to decibels
'''
import rsf.api as rsf
import numpy as np
import sys,math

par = rsf.Par()

verb = par.bool('verb',False) # verbosity flag
aref = par.float('aref',1.0)  # reference amplitude
inv = par.bool('inv',False) # invere transform (decibel to amplitude)

Fin = rsf.Input()         # input
n1 = Fin.int("n1")
o1 = Fin.float("o1")
d1 = Fin.float("d1")

nn = Fin.size(1); # number of traces

Fou = rsf.Output()        # output 
Fou.put("n1",n1)
Fou.put("o1",o1)
Fou.put('d1',d1)

din = np.zeros(n1,'f') 
dou = np.zeros(n1,'f')

for ii in range(nn):
    Fin.read(din)

    if inv:
        for i1 in range(n1):
            dou[i1] = aref * pow(10,0.05*din[i1])
    else:
        for i1 in range(n1):
            if math.fabs(din[i1])>1e-6:
                dou[i1] = 20 * math.log10( math.fabs(din[i1]/aref) )
            else:
                dou[i1] = -120

    Fou.write(dou)
    
# ------------------------------------------------------------
Fin.close()
Fou.close()

