#!/usr/bin/env python
'Second order statistics'

##   Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
##  
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##  
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##  
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

try:
	from numpy import *
	import rsf.api as rsf
except Exception, e:
	import sys
	print 'ERROR : need numpy'
	sys.exit(1)

par = rsf.Par()
pin = rsf.Input()
pout= rsf.Output()

mode=par.int("mode",0)
#	0	autocorrelation matrix \n	1	covariance matrix

nn=pin.shape()
nd=len(nn)

n1 = nn[-1]
n2 = nn[-2]
n3 = 1
for id in range(-nd,-2):
	n3 = n3 * nn[id]	

o1 = pin.float("o1")
d1 = pin.float("d1")

pout.put("n2",n1)
pout.put("o2",o1)
pout.put("d2",d1)

u1 = zeros((n2,n1),'f')


for i3 in range(n3):
	pin.read(u1)
	A = matrix(u1)
	if mode == 1:
		A = A - A.mean(1)
	u2 = A * A.H
	pout.write(array(u2).astype('f'))


