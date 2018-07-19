#!/usr/bin/env python
'''
Update sonic well log 

Use the initial sonic well log and an updated TDR to generate and updated sonic log and TDR
'''

import sys
import numpy as np
import rsf.api as rsf
import string
import math as mt

par = rsf.Par()
sonica = rsf.Input() # Input initial sonic well log

assert 'float' == sonica.type
sdepth = sonica.int("n1")

sonic = np.zeros(sdepth,'f')
sonica.read(sonic)

ms = par.int("ms") # (0 = Units of sonic in s); (1 = Units of sonic in ms)

switch = par.float("stretch") # (0 = Output TDR from input sonic log); (1 = Output updated sonic and TDR)
assert switch
dels = par.float("dels") # Depth step (units of m or ft)
assert dels

tdr = np.zeros(sdepth)

if (ms == 1):
    for i in range(1, sdepth):
        tdr[i] = tdr[i-1] + sonic[i]*(2*dels/1000000)
if (ms == 0):
    for i in range(1, sdepth):
        tdr[i] = tdr[i-1] + sonic[i]*(2*dels)

if (switch == 0):
    tdrFo = rsf.Output()
    tdrFo.write(np.array(tdr));

if (switch == 1):
    tdr_name = par.string('tdrNew')
    tdrNewa = rsf.Input(tdr_name) # Updated TDR
    assert 'float' == tdrNewa.type

    tdrNew = np.zeros(sdepth,'f')
    tdrNewa.read(tdrNew)

    sonicOld = np.zeros(sdepth)
    sonicNew = np.zeros(sdepth)
    warp = np.zeros(sdepth)

    sonicOut = np.zeros(sdepth)

    warp[0] = 1
    sonicOut[0] = sonic[0]

    for i in range(1, sdepth):
        sonicOld[i] = (tdr[i] - tdr[i-1])/(2*dels)
        sonicNew[i] = (tdrNew[i] - tdrNew[i-1])/(2*dels)

    for i in range(1, sdepth):
        warp[i] = 1 - (sonicOld[i]-sonicNew[i])/sonicOld[i]

        if (abs(warp[i] - warp[i-1]) > 0.1): #0.05
            warp[i] = warp[i-1]

        sonicOut[i] = warp[i]*sonic[i]

    tdrOut = np.zeros(sdepth)
    
    tdrOut[0] = tdrNew[0]

    if (ms == 1):
        for i in range(1, sdepth):
            tdrOut[i] = tdrOut[i-1] + sonicOut[i]*(2*dels/1000000)
    if (ms == 0):
        for i in range(1, sdepth):
            tdrOut[i] = tdrOut[i-1] + sonicOut[i]*(2*dels)

    tdrFo = rsf.Output() # Output updated TDR from updated sonic log
    sonicFo = rsf.Output('sonicFo') # Output updated sonic log

    tdrFo.write(np.array(tdrOut));
    sonicFo.write(np.array(sonicOut));
