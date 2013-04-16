#!/usr/bin/env python
'Add random noise using python.'

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
n1=nn[-1]
n2=nn[-2]
n3=1
for i1 in range(nd-2):
	n3=n3*nn[i1]

rangea=par.float("rangea",-1)# noise rangea (default=-1)
rangeb=par.float("rangeb",1)#  noise rangeb (default= 1)
seed=par.int("seed", n2) # random seed (default=n2)
rep=par.bool("rep",False) # if y, replace data with noise

sys.stderr.write('n1=%d n2=%d n3=%d\n'%(n1, n2, n3))

u1=zeros((n2,n1),'f')
u2=zeros((n2,n1),'f')

random.seed(seed)
for i3 in range(n3):
    pi.read(u1)
    for i2 in range(n2):
	for i1 in range(n1):	
	    if rep==False :    
		u2[i2,i1]=u1[i2,i1]+random.uniform(rangea,rangeb)
	    else:
		u2[i2,i1]=random.uniform(rangea,rangeb)
    po.write(u2)

po.close()


