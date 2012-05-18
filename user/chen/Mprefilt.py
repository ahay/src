#!/usr/bin/env python
'prefilter bank of pwd'

try:
	from numpy import *
	import rsf.api as rsf
	import rsf.pcmf as mf
except Exception, e:
	import sys
	print 'ERROR: numpy needed'
	sys.exit(1)


def fir(p,x):
	np=p.size
	n=(p.size-1)/2
	n1=x.size
	ny=n1-np
	y=zeros(n1,'f')
	for i1 in range(ny):
		y[i1+n]=inner(p,x[i1:i1+np])
	return y

par=rsf.Par()
input=rsf.Input()
output=rsf.Output()

nf=par.int("nf",1)
n1=input.int("n1")
n2=input.int("n2")

n3=2*nf+1
output.put("n3",2*nf+1)

c=mf.pcmf(nf)

#print c

x=zeros(n1*n2,'f')
input.read(x)
xx=x.reshape(n2,n1)

for i3 in range(n3):
	p=c[:,n3-i3-1]
	for i2 in range(n2):
		y=fir(p,xx[i2,:])
		output.write(y)


input.close()
output.close()


