#!/usr/bin/env python
'Fast PCA by iterative algorithm'

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

import sys
try:
	from numpy import *
	import rsf.api as rsf
except Exception, e:
	print 'ERROR : need numpy'
	sys.exit(1)

par = rsf.Par()
pin = rsf.Input()
pout= rsf.Output()

nc=par.int("nc",1)
#	component number
niter = par.int("niter",1)
#	iterations for each component
seed=par.int("seed",2012)
#	seed for random number

nn=pin.shape()
nd=len(nn)

n1 = nn[-1]
n2 = nn[-2]
n3 = 1
for id in range(-nd,-3):
	n3 = n3 * nn[id]	


pout.put("n2",nc)
pout.put("o2",0)
pout.put("d2",1)

u1 = zeros((n2,n1),'f')
u2 = zeros((nc,n1),'f')

random.seed(seed)

for i3 in range(n3):
	pin.read(u1)
	A = matrix(u1)
	A = A - A.mean(1)
	for ic in range(nc):
		p = matrix(random.randn(n1)).T
		for it in range(niter):
			t = matrix(zeros(n1,'f')).T
			for i2 in range(n2):
				t1 = A[i2,:].T*p
				t += A[i2,:]*t1 
			p = t/(sqrt(t.T*t))
		B[ic,:] = p
	u2 = array(B)
	sys.stderr.write('size=%d\n'%(u2.size))
	pout.write(u2)


