import numpy as np

def factorial(n):
	m=1
	for i1 in range(n):
		m=m*(i1+1)
	return m


def pcmf(n):
	N=2*n+1	
	c=np.zeros((N,N),dtype=float)
	c0=1
	for i1 in range(N):
		rt=np.zeros(2*n,dtype=float)
		k=i1-n
		c1=1.0*factorial(2*n)
		c1=c0*c1*c1/factorial(4*n)/factorial(n+k)/factorial(n-k)
		c0=-c0
		for i2 in range(n-k):
			rt[i2]=2*n-i2
		for i2 in range(n+k):
			rt[i2+n-k]=i2-2*n
#		print rt
		p=np.poly1d(rt,r=1)
		p=p*c1
#		print c1
		c[i1,0:N]=p
#	print c
	return c



