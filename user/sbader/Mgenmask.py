#!/usr/bin/env python

import sys
import numpy as np
import rsf.api as rsf
import string
import math as mt

par = rsf.Par()
dataa = rsf.Input()

fnt = par.int("beg")
bak = par.int("end")

assert 'float' == dataa.type

length = dataa.int("n1")

data = np.zeros(length,'f')
dataa.read(data)

data2 = np.array(data)

data2[0:fnt] = 0
data2[(length - bak):(length - 0)] = 0

log_eo = rsf.Output()
log_eo.write(np.array(data2));
