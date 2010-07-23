#!/usr/bin/env python

import numpy
import m8r

par = m8r.Par()
input  = m8r.Input()
output = m8r.Output()

n1 = input.int("n1") # trace length
n2 = input.size(1)   # number of traces

clip = par.float("clip")

trace = numpy.zeros(n1,'f')
for i2 in xrange(n2): # loop over traces
    input.read(trace)
    trace = numpy.clip(trace,-clip,clip)
    output.write(trace)

