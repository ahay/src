/* Fast C1 Coherence */

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

void fcoh1_acorr(float *d, float *v, int n, int nw)
/*< autorelation >*/
{
    int i, j, it;

    v[0] = 0.0;
    for(j=0; j<=nw; j++)
	v[0] += d[j]*d[j];

    for(i=1; i<n; i++)
    {
	v[i] = v[i-1];
	it = i - nw - 1;
	if(it >=0 && it < n)	
	    v[i] -= d[it] * d[it];
	it = i + nw;
	if(it >=0 && it < n)	
	    v[i] += d[it] * d[it];
    }
}

void xcorr(float *d1, float *d2, float **c, int n, int cmin, int cmax, int nw)
{/* c[n][m]  */ 
    int i1, i2, m, it, nb, old;
    float *buf, x;

    nb = 2*nw+1;
    m = cmax-cmin+1;
    buf = d1;

    for(i2=0; i2<n; i2++)
	for(i1=cmin; i1<=cmax; i1++)
	{
	    it = i1 + i2;
	    if(it>=0 && it < n)	c[i2][i1-cmin] = d1[i2] * d2[it];
	    else c[i2][i1-cmin] = 0.0;
	}
    for(i1=0; i1<m; i1++)
    {
	/* recursive mean filter */
	buf[nw] = c[0][i1];
	for(i2=1; i2<=nw; i2++)
	{
	    buf[nw-i2] = 0.0;
	    buf[nw+i2] = c[i2][i1];
	    c[0][i1] += c[i2][i1];
	}
	old = 0;
	for(i2=1; i2<n; i2++, old++)
	{
	    if(old==nb) old=0;
	    x = c[i2-1][i1] - buf[old];
	    c[i2][i1] = x;
	    it = i2 + nw;
	    if(it>=0 && it < n)
	    {
		c[i2][i1] += c[it][i1];
		buf[old] = c[it][i1];
	    }else buf[old] = 0.0;
	}
    }
}


float para_intp_max(float *u, int n)
{
    int imax, i;
    float a, b, x, fx; /* f(x) = a x^2 + b x + c */

    for(i=1, imax=0; i<n; i++)	if(u[i] > u[imax])  imax = i;

    if(n<3){
	fx = u[imax];
	u[0] = imax; 
	return fx;
    }
    if(imax == 0) imax = 1;
    if (imax == n-1) imax = n - 2;

	
/*	c = u[imax]; */
    a = (u[imax+1] + u[imax-1]) / 2 - u[imax];
    b = u[imax+1] - u[imax] - a;
    x = -0.5 * b*a / (a*a+0.00000000001);
    if(x+imax<0)
    {
	x = 0.0;
	fx = u[0];
    }else if(x+imax>n-1)
    {
	x = n-1;
	fx = u[n-1];
    }else{
	fx = u[imax] + 0.5*b*x; 
	x += imax;
    }
    u[0] = x;
    return fx;
}

typedef struct
{
    int n1, min, max, ntw;
    float **c;
}fcoh1;

void* fcoh1_init(int n, int pmin, int pmax, int nw)
/*< initialize >*/
{
    fcoh1 *p;

    p = sf_alloc(1,sizeof(fcoh1));
    p->n1 = n;
    p->min = pmin;
    p->max = pmax;
    p->ntw = nw;

    p->c = sf_floatalloc2(pmax-pmin+1, n);
    return p;
}

void fcoh1_close(void *h)
/*< release memory >*/
{
    fcoh1 *p = (fcoh1*) h;
    free(p->c[0]);
    free(p->c);
}

void fcoh1_tr(void *h, float *u1, float *u2, float *v1, float *v2, float *c)
/*< fast c1 for two neighbour trace:
  input:		u1	u2 	v1	v2
  output:		c1	u2	p1	v2
  >*/
{
    int i1, j1, it, m2;
    float var, coh;
    fcoh1 *p = (fcoh1*) h;

    m2 = p->max - p->min + 1;

    xcorr(u1, u2, p->c, p->n1, p->min, p->max, p->ntw);

    /* search integer maximum */
    for(i1=0; i1<p->n1; i1++)
    {
	for(j1=p->min; j1<=p->max; j1++)
	{
	    it = i1 + j1;
	    var = v1[i1];
	    if(it<0) var *= v2[0];
	    else if(it>=p->n1) var *= v2[p->n1-1];
	    else var *= v2[it];
	    p->c[i1][j1-p->min] /= sqrtf(var+0.000000000001);
	}
	coh = para_intp_max(p->c[i1], m2);
	v1[i1] = p->c[i1][0]+p->min;
	if(c==NULL) u1[i1] = coh;
	else c[i1]*=coh;
    }
}


