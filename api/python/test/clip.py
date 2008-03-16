#!/usr/bin/env python

import numpy
import rsf

par = rsf.Par()
input  = rsf.Input()
output = rsf.Output()
assert 'float' == input.type

n1 = input.int("n1")
n2 = input.size(1)
assert n1

clip = par.float("clip")
assert clip

trace = numpy.zeros(n1,'f')

for i2 in xrange(n2): # loop over traces
    input.read(trace)
    trace = numpy.clip(trace,-clip,clip)
    output.write(trace)

