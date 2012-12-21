#!/usr/bin/env python
'prefilter bank of pwd'

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
	import rsf.user.pcmf as mf
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

def divn(a,b,eta):
	n1=a.size
	x=zeros(n1,'f')
	na=0.0
	for i1 in range(n1):
#		na=((a[i1]*a[i1]+b[i1]*b[i1]))*eta
		na=(abs(a[i1])+abs(b[i1]))*eta
		x[i1]=b[i1]/(a[i1]+na*sign(a[i1]))
	return x

def dip_lop(n1, n2, p0, c, iter, eta):
	p=p0
	for it in range(iter):
		for i2 in range(n2):
			for i1 in range(n1):
				rr=polyval(c[:,i2,i1],p[i2,i1])
				dd=polyval(polyder(c[:,i2,i1]),p[i2,i1])
				na=(rr*rr+dd*dd)*eta*sign(dd)
#				na=(abs(rr)+abs(dd))*eta*sign(dd)
				dp=rr/(dd+na)
				p[i2,i1]=p[i2,i1]-0.5*dp
	return p

par=rsf.Par()
input=rsf.Input()
dip=rsf.Output()
coef=rsf.Output("coef")
pf=rsf.Output("pf")

nf=par.int("nf",1)
iter=par.int("iter",5)
eta=par.float("eta",0.05)
n1=input.int("n1")
n2=input.int("n2")

n3=2*nf+1
pf.put("n3",2*nf+1)
dip.put("n3",2*nf)
dip.put("n2",n2-1)
coef.put("n3",n3)
coef.put("n2",n2-1)

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

pf.write(y)

po=1
# dip steps
for i3 in range(n3):
	for i2 in range(n2-1):
		y[i3,i2,:]=y[i3,i2,:]-y[i3,i2+1,:]*po
	po=-po

coef.write(y[:,0:n2-1,:])

#p=-y[0,0:n2-1,:]/y[1,0:n2-1,:]
p=zeros((n2-1,n1),'f')

for i2 in range(n2-1):
	p[i2,:]=divn(-y[1,i2,:], y[0,i2,:], eta)

dip.write(p)
for i3 in range(n3-2):
	k=i3+3
	cc=y[0:k,0:n2-1,:]
	c=cc[::-1,:,:]
	p=dip_lop(n1, n2-1, p, c, iter, eta)
	dip.write(p)

input.close()
dip.close()
coef.close()
pf.close()


