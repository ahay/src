
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

import numpy as np

def factorial(n):
	m=1
	for i1 in range(n):
		m=m*(i1+1)
	return m

# B(Z) = c(-n:n, 0:2n)p^mZ^k
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
		c[N-1-i1,::-1]=p
#	print c
	return c



