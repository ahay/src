/* edge-preserving smoothing by mean filter */

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


typedef struct{
	int n;
	int nfw;
	float *v, *e;
}EPSMEAN;

void* epsmean_init(int n1, int nfw)
/*< initialize >*/
{
	EPSMEAN *p;

	p = sf_alloc(1, sizeof(EPSMEAN));
	p->n = n1;
	p->nfw = nfw;
	p->v = sf_floatalloc(n1);
	p->e = sf_floatalloc(n1);

	return p;
}

void epsmean(void *h, float *x, int d)
/*< eps by mean filter >*/
{
	int i1, j1, j2;
	EPSMEAN *p;
	p = (EPSMEAN*) h;
	float t;
	for(i1=0; i1 < p->n; i1++)
	{
		p->v[i1] = 0.0;
		p->e[i1] = 0.0;
	}

	for(i1=0; i1 < p->n; i1++)
	{
		for(j1=0; j1 < p->nfw; j1++)
			p->v[i1] += x[(j1+i1)*d];
		p->v[i1] /= p->nfw;
		for(j1=0; j1 < p->nfw; j1++)
		{
			t = x[(i1+j1)*d] - p->v[i1];
			p->e[i1] += t*t;
		}
	}

	for(i1=p->nfw; i1 < p->n - p->nfw; i1++)
	{
		j2 = -p->nfw;
		for(j1=-p->nfw; j1<=p->nfw; j1++)
			if(p->e[i1+j1] < p->e[i1+j2]) j2 = j1;
		x[i1*d] = p->v[i1+j2];
	}
}

