#!/usr/bin/env python
'shuffle the data'
'This code is the same one as RSFSRC/user/chen/Mshuffle.py. It resides in here '
'to test the configuration of SConstruct for python environment is correct or not.'
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

par=rsf.Par()
pi=rsf.Input()
po=rsf.Output()

nn=pi.shape()
nd=len(nn)

axis=par.int("axis",2)
if axis>nd:
	sys.stderr.write('axis=%d greater than nd=%d'%(axis, nd))
n2=nn[-axis]
seed=par.int("seed", n2)
inv=par.bool("inv", False)

n1=1
n3=1
for i1 in range(axis-1):
	n1=n1*nn[-1-i1]
for i1 in range(nd-axis):
	n3=n3*nn[i1]

sys.stderr.write('n1=%d n2=%d n3=%d\n'%(n1, n2, n3))

u1=zeros((n2,n1),'f')
u2=zeros((n2,n1),'f')

random.seed(seed)
j1=range(n2)
j2=j1
random.shuffle(j1)

for i3 in range(n3):
	pi.read(u1)
	for i2 in j2:
		if inv:
			u2[j1[i2],:]=u1[i2,:]
		else:	
			u2[i2,:]=u1[j1[i2],:]
	po.write(u2)

po.close()


