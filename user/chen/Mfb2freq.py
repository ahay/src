#!/usr/bin/env python
'Frequency response of OPWD filter bank'

try:
	from numpy import *
	import rsf.api as rsf
	import rsf.pcmf as mf
except Exception, e:
	import sys
	print 'ERROR: numpy needed'
	sys.exit(1)

par=rsf.Par()
outp=rsf.Output()

n1=par.int("n1",100) # frequency samples
d1=0.5/n1
o1=-d1*n1

nf=par.int("nf",1) # filter order
nn=2*nf+1
outp.put("n1", 2*n1+1)
outp.put("o1", o1)
outp.put("d1", d1)
outp.put("n2", 2*n1+1)
outp.put("o2", o1)
outp.put("d2", d1)
outp.put("n3",nn)
outp.put("o3", 0)
outp.put("d3", 1)
outp.put("n4",nn)
outp.put("o4", 0)
outp.put("d4", 1)

c=mf.pcmf(nf)
f=o1+d1*arange(2*n1+1)
z=exp(-2.0*pi*f*1j)

for i2 in range(nn):
	d2=pow(z,-nf)*polyval(c[:,i2],z)
	for i1 in range(nn):
		d1=pow(z,-nf)*polyval(c[:,i1],z)
		dd=outer(d1,d2)
		dd=abs(dd).astype('f')
		outp.write(dd)


