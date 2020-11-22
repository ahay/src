/* edge-preserving smoothing by mean filter */

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
#include "sfilt.h"

typedef struct{
	int n;
	int nfw;
	float *v, *e, *u;
}EPSMEAN;

void* epsmean_init(int n1, int nfw)
/*< initialize >*/
{
	EPSMEAN *p;

	p = sf_alloc(1, sizeof(EPSMEAN));
	p->n = n1;
	p->nfw = nfw;
	p->u = sf_floatalloc(n1);
	p->v = sf_floatalloc(n1);
	p->e = sf_floatalloc(n1);

	return p;
}

static float mean_var(int n, float* x, float* mean)
{
	int i;
	float u, v=0.;
	
	for(i=0, u=0.0; i<n; i++)
		u += x[i];
	u /= n;
	*mean = u;
	for(i=0, u=0.0; i<n; i++)
	{
		u = x[i] - mean[0];
		v = u*u;
	}
	return v/n;
}	


#define MY_MAX(a,b) ((a)<(b) ? (b) : (a))
#define MY_MIN(a,b) ((a)>(b) ? (b) : (a))
void epsmean(void *h, float *x, int d)
/*< eps by mean filter >*/
{
	int i1, j1, j2, min, max, l;
	EPSMEAN *p;
	p = (EPSMEAN*) h;

	for(i1=0; i1 < p->n; i1++)
	p->u[i1] = x[i1*d];

	for(i1=0; i1 < p->n; i1++)
	{
		min = i1;
		max = MY_MIN(i1+p->nfw, p->n-1);
		l = max - min + 1;
		p->e[i1] = mean_var(l, p->u+min, p->v+i1);
	}

	// selection
	for(i1=0; i1 < p->n; i1++)
	{
		min = MY_MAX(i1-p->nfw, 0);
		max = MY_MIN(i1, p->n-p->nfw-1);
		j2 = min;
		for(j1=min+1; j1 <= max; j1++) 
		if(p->e[j1] < p->e[j2]) j2 = j1;
		x[i1*d] = p->v[j2];
	}
}

