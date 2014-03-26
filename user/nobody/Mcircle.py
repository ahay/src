#!/usr/bin/env python
'generate circle wave snapshot'

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

par=rsf.Par()
output=rsf.Output()


n1=par.int("n1",500)
o1=par.float("o1",-250)
d1=par.float("d1",1)
n2=par.int("n2",500)
o2=par.float("o2",-250)
d2=par.float("d2",1)
f=par.float("f",0.05)
# frequency (circles per sampling interval)
df=par.float("df",0.0)
# chirp frequceny shift
a=par.float("a",0.0)
# amplitude attenuation

output.put("n1",n1)
output.put("o1",o1)
output.put("d1",d1)
output.put("n2",n2)
output.put("o2",o2)
output.put("d2",d2)

dat=zeros((n1,n2),'f')

for i1 in range(n1):
	x1=i1*d1+o1
	for i2 in range(n2):
		x2=i2*d2+o2
		r=sqrt(x1*x1+x2*x2)
		p=2*math.pi*(f+df*r)*r
		dat[i1,i2]=exp(-a*r)*cos(p)

output.write(dat)
#output.close()


