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

def dip_lop(n1, n2, p0, c, iter):
	p=p0
	for it in range(iter):
		for i2 in range(n2):
			for i1 in range(n1):
				r=polyval(c[:,i2,i1],p[i2,i1])
				d=polyval(polyder(c[:,i2,i1]),p[i2,i1])
				dp=r/d
				p[i2,i1]=p[i2,i1]-0.5*dp
	return p

par=rsf.Par()
input=rsf.Input()
dip=rsf.Output("dip")
output=rsf.Output()

nf=par.int("nf",1)
iter=par.int("iter",5)
n1=input.int("n1")
n2=input.int("n2")

n3=2*nf+1
output.put("n3",2*nf+1)
dip.put("n3",2*nf)

c=mf.pcmf(nf)

#print c

x=zeros(n1*n2,'f')
y=zeros((n3,n2,n1),'f')
input.read(x)
xx=x.reshape(n2,n1)


for i3 in range(n3):
	p=c[:,n3-i3-1]
	for i2 in range(n2):
		x=xx[i2,:]
		yy=fir(p,x)
		y[i3,i2,:]=yy

output.write(y)

po=1
# dip steps
for i3 in range(n3):
	for i2 in range(n2-1):
		y[i3,i2,:]=y[i3,i2,:]-y[i3,i2+1,:]*po
	po=-po

p=-y[1,:,:]/y[0,:,:]
dip.write(p)
for i3 in range(n3-2):
	k=i3+3
	c=y[0:k,:,:]
	p=dip_lop(n1, n2, p, c, iter)
	dip.write(p)

input.close()
dip.close()
output.close()


