#!/usr/bin/env python

import numpy
import rsf.api as rsf
try:
    from builtins import range
except:
    from __builtin__ import range

# input and output files
par = rsf.Par()
input  = rsf.Input()
output = rsf.Output()

#check that the input type is float
assert 'float' == input.type

#get parameters from the rsf head
n1 = input.int("n1")
d1 = input.float("d1")
n2 = input.size(1)
assert n1

#set parameters for output file
output.put("n1",n1)
output.put("d1",d1)
output.put("label1","Time")
output.put("unit1","s")

#get parameters from command line
upper = par.float("upper")
assert upper

lower = par.float("lower")
assert upper

# call numpy
trace = numpy.zeros(n1,'f')

# loop over traces
for i2 in range(n2):
    input.read(trace)
    trace = numpy.clip(trace,lower,upper)
    output.write(trace)

output.close()
input.close()
