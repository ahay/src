#!/usr/bin/env python
'shuffle the data'

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

from __future__ import print_function
import sys
try:
        from numpy import *
        import rsf.api as rsf
except Exception as e:
        print('ERROR : need numpy')
        sys.exit(1)

par=rsf.Par()
pi=rsf.Input()
pi2=rsf.Input("pi2")

po=rsf.Output()
po2=rsf.Output("po2")

nn=pi.shape()
nd=len(nn)

axis=par.int("axis",2)
if axis>nd:
        sys.stderr.write('axis=%d greater than nd=%d'%(axis, nd))
#n2=nn[-axis]
seed=par.int("seed")
#inv=par.bool("inv", False)

n1 = pi.int("n1")
n2 = pi.int("n2")
#n3 = pi.int("n3")

po2.put('n1',1)
po2.put('n2',n1)

po.put('n2',n1)
po.put('n1',n2)

#n1=1
#n3=1
#for i1 in range(axis-1):
        #n1=n1*nn[-1-i1]
#for i1 in range(nd-axis):
        #n3=n3*nn[i1]

#sys.stderr.write('n1=%d n2=%d n3=%d\n'%(n1, n2, n3))

u1=zeros((n2,n1),'f')
u2=zeros((1,n1),'f')

#sys.stderr.write('n1=%d'%(pi2.shape[0]))
#sys.stderr.write('n2=%d'%(pi1.shape[1]))

random.seed(seed)
#j1=list(range(n2))
#j2=j1

pi.read(u1)
pi2.read(u2)

u1=transpose(u1,(1,0))
u2=transpose(u2,(1,0))

data = list(zip(u1,u2))
random.shuffle(data)

#for i3 in range(n3):
        #pi.read(u1)
        #for i2 in j2:
                #if inv:
                        #u2[j1[i2],:]=u1[i2,:]
                #else:
                        #u2[i2,:]=u1[j1[i2],:]
        #po.write(u2)


for xi,yi in data:
	#sys.stderr.write('n1=%d'%(yi.shape[0]))
	po.write(xi)
	po2.write(yi)

po.close()
po2.close()
