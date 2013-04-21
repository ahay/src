/* Fast PCA */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>

void fpca1(float **a, int n1, int n2, float *p, float *t, int niter)
{
	int i1, i2, it, inc=1;
	float t1;
	for(it=0; it<niter; it++)
	{
		for(i1=0; i1<n1; i1++) 	t[i1] = 0.0;
		for(i2=0; i2<n2; i2++)
		{
			t1 = sdot_(&n1, a[i2], &inc, p, &inc);
			saxpy_(&n1, &t1, a[i2], &inc, t, &inc);
		}
		t1 = sqrtf(sdot_(&n1, t, &inc, t, &inc) );
		for(i1=0; i1<n1; i1++) 	p[i1] = t[i1]/t1;
	} 
}


void fpca(float **a, int n1, int n2, float **b, int m2, int niter)
/*< fpca >*/
{
	float *v1, *p, *t;
	int i1, i2, j2;

	v1 = sf_floatalloc(n1);
	p = sf_floatalloc(n1);
	t = sf_floatalloc(n1);

	// centering
	for(i1=0, v1[i1]=0.0; i1<n1; i1++)
	{
		for(i2=0; i2<n2; i2++) v1[i1] += a[i2][i1];
		v1[i1] /= n2;
	}
	for(i2=0; i2<n2; i2++) 
	for(i1=0; i1<n1; i1++) 
		a[i2][i1] -= v1[i1];

	for(j2=0; j2<m2; j2++)
	{
		for(i1=0; i1<n1; i1++) 
			p[i1] = 2.0*rand()/RAND_MAX-1.0;
		fpca1(a, n1, n2, p, t, niter);
		for(i2=0; i2<n2; i2++) 
		for(i1=0; i1<n1; i1++) 
			a[i2][i1] -= p[i1];
	}
}


