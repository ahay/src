#!/usr/bin/env python

import sys
import numpy
import m8r

par = m8r.Par()
inp  = m8r.Input()
output = m8r.Output()
assert 'float' == inp.type

n1 = inp.int("n1")
n2 = inp.size(1)
assert n1

clip = par.float("clip")
assert clip

trace = numpy.zeros(n1,'f')

for i2 in range(n2): # loop over traces
    inp.read(trace)
    trace = numpy.clip(trace,-clip,clip)
    output.write(trace)

sys.exit(0)
