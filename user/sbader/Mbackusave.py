#!/usr/bin/env python
'C-Wave Backus Averaging (See Marion et al., 1994)'

import sys
import numpy as np
import rsf.api as rsf
import string
import math as mt

# Inputs
par = rsf.Par()
slowness_name = par.string('slowness') # Slowness from Logs
rhob_name = par.string('density') # Density from Logs

deptha = rsf.Input() # Depth from Logs
slowa = rsf.Input(slowness_name)
rhoba = rsf.Input(rhob_name)

assert 'float' == deptha.type
assert 'float' == slowa.type
assert 'float' == rhoba.type

sdepth = deptha.int("n1")
sslow = slowa.int("n1")
srhob = rhoba.int("n1")

depth = np.zeros(sdepth,'f')
slow = np.zeros(sdepth,'f')
rhob = np.zeros(sdepth,'f')

deptha.read(depth)
slowa.read(slow)
rhoba.read(rhob)

ratio = par.float("ratio") # Percent of dom wavelength
peak_f = par.float("peak_f") # Dom wavelength
depthsample = par.float("depthsample") # Depth Sampling

depth_bkn = rsf.Output() # Output depth sampling
vel_bkn = rsf.Output('vel_bk') # Backus Avg. velocity
slow_bkn = rsf.Output('slow_bk') # Backus Avg. slowness
rhob_bkn = rsf.Output('rhob_bk') # Averaged density

assert sdepth == sslow
assert sdepth == srhob
assert sslow == srhob

count = 0
mask = np.zeros(sdepth)
for i in range(sdepth):
    if ((rhob[i] > 0) & (slow[i] > 0)):
        mask[i] = 1;
        count  = count + 1

slow = [a*b for a,b, in zip(slow,mask)]
rhob = [a*b for a,b, in zip(rhob,mask)]

vel = np.zeros(sdepth)

mean = 0
temp = 0
for i in range(sdepth):
    if (slow[i] > 0):
        vel[i] = 1000000/slow[i]
        mean = mean + vel[i]
        temp = temp + 1
    else:
        vel[i] = 0;

mean = mean/temp

wave = mean/peak_f
N = round((1/depthsample)*ratio*wave);

if (mt.floor(N/2) != mt.ceil(N/2)):
    N = N + 1

fk = 1/N

N = (int)(N)

velf = np.trim_zeros(vel)
depthf = np.zeros(sdepth)
mask2 = np.zeros(sdepth)

for i in range(count):
    for j in range(i, sdepth):
        if (velf[i] == vel[j]):
            depthf[j] = depth[j]

minval = np.min([x for x in depthf if x !=0])
maxval = np.max([x for x in depthf if x !=0])

for i in range(sdepth):
    if(depthf[i] == minval):
        mini = i
    if(depthf[i] == maxval):
        maxi = i

for i in range(mini+(N/2),maxi-(N/2)):
    mask2[i] = 1

pave = np.zeros(sdepth)
term = np.zeros(sdepth)

for i in range((N/2),sdepth-(N/2)):
    for j in range((i-(N/2)),i+(N/2)-1):
        pave[i] = pave[i] + fk*rhob[j]
    for j in range((i-(N/2)),i+(N/2)-1):
        if ((rhob[j] > 0) & (vel[j] > 0)):
                term[i] = term[i] + fk*(1/(rhob[j]*vel[j]*vel[j]))

pave = [a*b for a,b, in zip(pave,mask2)]
term = [a*b for a,b, in zip(term,mask2)]

vEMT_bk = np.zeros(sdepth)
sEMT_bk = np.zeros(sdepth)

for i in range(sdepth):
    sEMT_bk[i] = mt.sqrt(pave[i]*term[i])
    if (term[i] != 0):
        vEMT_bk[i] = 1/mt.sqrt(pave[i]*term[i])

depth_bkn.write(np.array(depthf))
vel_bkn.write(np.array(vEMT_bk))
slow_bkn.write(np.array(sEMT_bk))
rhob_bkn.write(np.array(pave))


depth_bkn.close()
vel_bkn.close()
slow_bkn.close()
rhob_bkn.close()
