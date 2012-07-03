#!/usr/bin/env python
'Random sort'

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

id=par.int("id",2)
if id>nd:
	sys.stderr.write('id=%d greater than nd=%d'%(id, nd))
id=id-1
n2=nn[id]
seed=par.int("seed", n2)
inv=par.bool("inv", False)

n1=1
n3=1
for i1 in range(id):
	n1=n1*nn[i1]
for i1 in range(nd-id-1):
	n3=n3*nn[id+i1+1]


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


