#!/usr/bin/env python
'''
find max value in a file
'''
import rsf.api as rsf
import numpy as np
import sys

par = rsf.Par()

verb = par.bool('verb',False) # verbosity flag

Fin = rsf.Input()         # input
n1 = Fin.int  ("n1")
nn = Fin.size(1); # number of traces
din = np.zeros(n1,'f') 

Fou = rsf.Output()        # output 
Fou.put("n1",1)
Fou.put("o1",0)
Fou.put('d1',1)
Fou.put("n2",1)
Fou.put("n3",1)
Fou.put("n4",1)

dou = np.zeros(1,'f')

mymax=0
for i in range(nn):
    Fin.read(din)

    for i1 in range(n1):
        if abs(din[i1])>abs(mymax):
            mymax=din[i1]
        
dou[0]=mymax
Fou.write(dou)

print  >> sys.stderr,'max=',mymax

# ------------------------------------------------------------
Fin.close()
Fou.close()

