/* edge-preserving smoothing by polynomial filter */

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
#include "polyfit.h"

typedef struct{
	int n;
	int nfw;
	int order;
	float *c, **u, *e;
}EPSPOLY;

void* epspoly_init(int n1, int nfw, int order)
/*< initialize >*/
{
	EPSPOLY *p;

	if(order > nfw) 
	sf_error("order=%d > window size in epspoly", order);
	p = sf_alloc(1, sizeof(EPSPOLY));
	p->n = n1;
	p->nfw = nfw;
	p->order = order;
	p->c = sf_floatalloc(order);
	p->e = sf_floatalloc(n1);
	p->u = sf_floatalloc2(nfw+1, n1);

	return p;
}

#define MY_MAX(a,b) ((a)<(b) ? (b) : (a))
#define MY_MIN(a,b) ((a)>(b) ? (b) : (a))
void epspoly(void *h, float *x, int d)
/*< eps by polynomial fitting >*/
{
	int i1, j1, j2, min, max;
	EPSPOLY *p;
	float t, *pv;
	void *h1;

	p = (EPSPOLY*) h;
	h1 = polyfit_init(p->nfw+1, p->order, 0, 0);

	for(i1=0; i1 < p->n; i1++)
	{
		pv = p->u[i1];
		for(j1=0; j1<=p->nfw; j1++)
		if(i1+j1 < p->n)	pv[j1] = x[(j1+i1)*d];
		else pv[j1] = 0.0;
		polyfit_coef(h1, pv, p->c);
		polyfit(h1, p->c, pv);
		p->e[i1] = 0.0;
		for(j1=0; j1<=p->nfw; j1++)
		{
			t = pv[0] - x[(i1+j1)*d];
			p->e[i1] += t*t;
		}
		p->e[i1] /= p->nfw;
	}

	polyfit_close(h1);

	for(i1=0; i1 < p->n; i1++)
	{
		min = MY_MAX(i1-p->nfw, 0);
		max = MY_MIN(i1, p->n-p->nfw-1);
		j2 = min;
        for(j1=min+1; j1 <= max; j1++)
        if(p->e[j1] < p->e[j2]) j2 = j1;

		pv = p->u[j2];
		x[i1*d] = pv[i1-j2];
	}
}

