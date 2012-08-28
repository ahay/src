/* edge-preserving smoothing by polynomial filter */

/*
  Copyright (C) 2012 University of Texas at Austin
  
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
	float *c, **u, **e;
	void *h;
}EPSPOLY;

void* epspoly_init(int n1, int nfw, int order)
/*< initialize >*/
{
	EPSPOLY *p;

	if(order >= nfw) sf_error("order=%d >= nfw=%d in epspoly", order, nfw);
	p = sf_alloc(1, sizeof(EPSPOLY));
	p->n = n1;
	p->nfw = nfw;
	p->order = order;
	p->c = sf_floatalloc(order);
	p->u = sf_floatalloc2(nfw, n1);
	p->e = sf_floatalloc2(nfw, n1);

	p->h = polyfit_init(nfw, order, 0, 0);
	return p;
}

void epspoly(void *h, float *x, int d)
/*< eps by polynomial fitting >*/
{
	int i1, j1, j2;
	EPSPOLY *p;
	p = (EPSPOLY*) h;
	float t;

	for(i1=0; i1 < p->n - p->nfw; i1++)
	{
		for(j1=0; j1<p->nfw; j1++) p->u[i1][j1] = x[(i1+j1)*d];
		polyfit_coef(p->h, p->u[i1], p->c);
		polyfit(p->h, p->c, p->u[i1]);
		for(j1=0; j1<p->nfw; j1++) 
		{
			t = p->u[i1][j1]-x[(i1+j1)*d];
			p->e[i1][j1] = t*t;
		}
	}

	for(i1=0; i1 < p->nfw; i1++)
	{
		j2 = 0;
		for(j1=0; j1<=i1; j1++)
		if(p->e[j1][i1-j1] < p->e[j2][i1-j2]) j2 = j1;
		x[i1*d] = p->u[j1][i1-j2];
	}
	for(i1=p->nfw; i1 < p->n - p->nfw; i1++)
	{
		j2 = 0;
		for(j1=0; j1<p->nfw; j1++)
		if(p->e[j1+i1-p->nfw][p->nfw-j1-1] < p->e[j2+i1-p->nfw][p->nfw-j2-1]) 
			j2 = j1;
		x[i1*d] = p->u[j2+i1-p->nfw][p->nfw-j2-1];
	}
}

