#!/usr/bin/env python

import numpy as np
import m8r
import sys

par = m8r.Par()
inp  = m8r.Input()
output = m8r.Output("out")
assert 'float' == inp.type

n1 = inp.int("n1")
n2 = inp.leftsize(1)

indx_fldr=inp.get_segy_keyindx('fldr')

sys.stderr.write('indx_fldr=%s/n'%repr(indx_fldr))

while True:
    eof_get_tah, intrace, inheader = inp.get_tah()
    if eof_get_tah:
        break

    fldr=inheader[indx_fldr]

    if fldr<700:
        fldr+=200
    inheader[indx_fldr]=fldr

    output.put_tah(intrace,inheader)
