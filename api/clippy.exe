#!/usr/bin/env python
import sys, os
import numarray

root = os.environ.get('RSFROOT')
libdir = os.path.join(root,'lib')
sys.path.append(libdir)

import rsf

par = rsf.Par()
input = rsf.Input()
output = rsf.Output()

assert 'float' == input.type

n1 = input.int("n1")
assert n1
n2 = input.size(1)

clip = par.float("clip")

trace = numarray.zeros(n1,'f')

for i2 in xrange(n2):
    input.read(trace)
    trace = numarray.clip(trace,-clip,clip)
    output.write(trace)

