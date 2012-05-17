#!/usr/bin/env python
'generate circle wave snapshot'

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
f=par.float("f",5)

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
		p=2*math.pi*f/180*sqrt(x1*x1+x2*x2)
		dat[i1,i2]=cos(p)

output.write(dat)
output.close()


