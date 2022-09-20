#!/usr/bin/env python

import sys
import numpy as np
import m8r

def slow_median(data):
    'find median by slow sorting, changes data'

    n = len(data)
    for k in range(n):
        item1 = data[k]

        # assume everything up to k is sorted
        for k2 in range(k,-1,-1):
            item2 = data[k2-1]
            if item1 >= item2:
                break
            data[k2] = item2

        data[k2] = item1
    
    return data[n/2]

# initialization
par = m8r.Par()
inp = m8r.Input()
out = m8r.Output()

# get data dimensions
n1 = inp.int('n1')
n2 = inp.int('n2')
n3 = inp.leftsize(2)

# input and output 
data = np.zeros([n2,n1],'f')
signal = np.zeros([n2,n1],'f')

# sliding window
w1 = par.int('w1',5)
w2 = par.int('w2',5)

nw = w1*w2
win = np.zeros([w2,w1],'f')

how = par.string('how','fast')
# what to compute 
# (fast median, slow median, mean) 

for i3 in range(n3):
    # read data plane
    inp.read(data)

    for i2 in range(n2):
        sys.stderr.write("\r%d of %d" % (i2,n2))
        s2 = max(0,min(n2-w2,i2-w2/2-1))
        for i1 in range(n1):
            s1 = max(0,min(n1-w1,i1-w1/2-1))

	    # copy window
            win = data[s2:s2+w2,s1:s1+w1]

            if how[0] == 'f': # fast median
                signal[i2,i1] = np.median(win)
            elif how[0] == 's': # slow median
                signal[i2,i1] = slow_median(win.flatten())
            elif how[0] == 'm': # mean
                # !!! ADD CODE !!! 
                pass
            else:
                sys.stderr.write(
                    "Unknown method \"%s\"\n" % how)
                sys.exit(1)
    sys.stderr.write("\n") 	
    # write out
    out.write(signal)

sys.exit(0)
