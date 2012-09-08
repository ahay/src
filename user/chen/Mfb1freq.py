#!/usr/bin/env python
'Frequency response of PWD filter bank'

try:
	from numpy import *
	import rsf.api as rsf
	import rsf.user.pcmf as mf
except Exception, e:
	import sys
	print 'ERROR: numpy needed'
	sys.exit(1)

par=rsf.Par()
outp=rsf.Output()

n1=par.int("n1",100) # frequency samples

nf=par.int("nf",1) # filter order
n2=2*nf+1

outp.put("n1",n1)
outp.put("o1",0.0)
outp.put("d1",0.5/(n1-1))
outp.put("n2",n2)
#outp.settype('complex')

c=mf.pcmf(nf)
f=0.5*arange(n1)/(n1-1)
z=exp(-2.0*pi*f*1j)

for i2 in range(n2):
	p=poly1d(c[0:n2,i2],r=0)
	d=pow(z,-nf)*polyval(p,z)
	dd=abs(d).astype('f')
	outp.write(dd)


