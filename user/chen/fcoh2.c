/* Fast C2 Coherence */

/*
  Copyright (C) 2013 Zhonghuan Chen, Tsinghua University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranyw of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <rsf.h>

float para_intp_max(float *u, int n, int d)
{
	int imax, i;
	float a, b, x, fx; /* f(x) = a x^2 + b x + c */

	for(i=1, imax=0; i<n; i++)	if(u[i*d] > u[imax*d])  imax = i;

	if(n<3){
		fx = u[imax*d];
		u[0] = imax; 
		return fx;
	}
	if(imax == 0) imax = 1;
	if (imax == n-1) imax = n - 2;

	
/*	c = u[imax]; */
	a = (u[(imax+1)*d] + u[(imax-1)*d]) / 2 - u[imax*d];
	b = u[(imax+1)*d] - u[imax*d] - a;
	x = -0.5 * b*a / (a*a+0.00000000001);
	if(x+imax<0)
	{
		x = 0.0;
		fx = u[0];
	}else if(x+imax>n-1)
	{
		x = n-1;
		fx = u[(n-1)*d];
	}else{
		fx = u[imax*d] + 0.5*b*x; 
		x += imax;
	}
	u[0] = x;
	return fx;
}

typedef struct
{
	int n1, n2, min, max, ntw;
	float ***d1, ***d2;
}fcoh2;

void* fcoh2_init(int n1, int n2, int pmin, int pmax, int nw)
/*< initialize >*/
{
	fcoh2 *p;

	p = sf_alloc(1,sizeof(fcoh2));
	p->n1 = n1;
	p->n2 = n2;
	p->min = pmin;
	p->max = pmax;
	p->ntw = nw;

	p->d1 = sf_floatalloc3(n1, n2, pmax-pmin+1);
	p->d2 = sf_floatalloc3(n1, n2, pmax-pmin+1);
	return p;
}

void fcoh2_close(void *h)
/*< release memory >*/
{
	fcoh2 *p = (fcoh2*) h;
	free(p->d1[0][0]);
	free(p->d1[0]);
	free(p->d1);
	free(p->d2[0][0]);
	free(p->d2[0]);
	free(p->d2);
}

void tx2txp(float **u1, float ***d1, int nxw, int n1, int n2, int m3, int n3)
{
	int i1, i2, i3, it;
/*	for(i3=m3; i3<=n3; i3++)
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		d1[i3][i2][i1] = 0.0;
		for(it=-nxw; it<=nxw; it++)
		if(i2+it>=0 && i2+it < n2 && i1+it*i3 >=0 && i1+it*i3 <n1)
		d1[i3][i2][i1] += u1[i2+it][i1+it*i3];
	}
	return ;
*/
	for(i3=m3; i3<=n3; i3++)
	for(i1=0; i1<n1; i1++)
	{
		d1[i3][0][i1] = u1[0][i1];
		for(i2=1; i2<=nxw; i2++)
		{
			it = i1 + i3*i2;
			if(it>=0 && it <n1)
				d1[i3][0][i1] += u1[i2][it];	
		}
	}
	for(i3=m3; i3<=n3; i3++)
	for(i2=1; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		it = i1 - i3;
		if(it>=0 && it < n1){
			d1[i3][i2][i1] = d1[i3][i2-1][it];	
			if(i2-nxw-1 >=0 && i2-nxw-1 < n2 &&
				i1-(nxw+1)*i3 >=0 && i1-(nxw+1)*i3 < n1 ) 
				d1[i3][i2][i1] -= u1[i2-nxw-1][i1-(nxw+1)*i3];	
			if(i2+nxw >=0 && i2+nxw < n2 &&
				i1+nxw*i3 >=0 && i1+nxw*i3 < n1 ) 
				d1[i3][i2][i1] += u1[i2+nxw][i1+nxw*i3];
		}else{ 
			d1[i3][i2][i1] = 0.0;	
			for(it=-nxw; it<=nxw; it++)
			if(i2+it>=0 && i2+it < n2 && i1+it*i3 >=0 && i1+it*i3 <n1)
			d1[i3][i2][i1] += u1[i2+it][i1+it*i3];
		}
	}

}

void rmean(float *d, int n, int nw, float *buf)
{
	int i1, i2, old, it;
	float x;

	buf[nw] = d[0];
	for(i2=1; i2<=nw; i2++)
	{
		buf[nw-i2] = 0.0;
		buf[nw+i2] = d[i2];
		d[0] += d[i2];
	}
	old=0;
	for(i1=1; i1<n; i1++, old++)
	{
		if(old==2*nw+1) old=0;
		x = d[i1-1] - buf[old];
		d[i1] = x;
		it = i1 + nw;
		if(it>=0 && it < n)
		{
			d[i1] += d[it];
			buf[old] = d[it];
		}else buf[old] = 0.0;
	}
}

void fcoh2_2d(void*h, float **u1, float **u2, int nxw)
/*< fast c1 for two neighbour trace:
input:		u1	0 	
output:		c1  p1
>*/
{
	int i1, i2, i3, m1;
	fcoh2 *p;

	p = (fcoh2*) h;

	m1 = p->max - p->min + 1;

	for(i2=0; i2<p->n2; i2++)
	for(i1=0; i1<p->n1; i1++)
		u2[i2][i1] = u1[i2][i1] * u1[i2][i1];

	tx2txp(u1, p->d1-p->min, nxw, p->n1, p->n2, p->min, p->max);
	tx2txp(u2, p->d2-p->min, nxw, p->n1, p->n2, p->min, p->max);

	for(i3=0; i3<m1; i3++)
	for(i2=0; i2<p->n2; i2++)
	{
		for(i1=0; i1<p->n1; i1++)
			p->d1[i3][i2][i1] = p->d1[i3][i2][i1] * p->d1[i3][i2][i1];
		rmean(p->d1[i3][i2], p->n1, p->ntw, u1[0]);
		rmean(p->d2[i3][i2], p->n1, p->ntw, u1[0]);
		for(i1=0; i1<p->n1; i1++)
			p->d1[i3][i2][i1] /= (p->d2[i3][i2][i1]*m1);
	}
	for(i2=0; i2<p->n2; i2++)
	for(i1=0; i1<p->n1; i1++)
	{
		u1[i2][i1] = para_intp_max(p->d1[0][i2]+i1, m1, p->n1*p->n2);
		u2[i2][i1] = p->d1[0][i2][i1]+p->min;
	}
}


