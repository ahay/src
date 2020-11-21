/* edge-preserving smoothing by median filter */

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
	float *u, *v, *e;
}EPSMEDIAN;

void* epsmedian_init(int n1, int nfw)
/*< initialize >*/
{
	EPSMEDIAN *p;

	p = sf_alloc(1, sizeof(EPSMEDIAN));
	p->n = n1;
	p->nfw = nfw;
	p->u = sf_floatalloc(nfw+1);
	p->v = sf_floatalloc(n1);
	p->e = sf_floatalloc(n1);
	return p;
}

#define MY_MAX(a,b) ((a)<(b) ? (b) : (a))
#define MY_MIN(a,b) ((a)>(b) ? (b) : (a))
void epsmedian(void *h, float *x, int d)
/*< eps by median filter >*/
{
	int i1, j1, j2, min, max, l;
	EPSMEDIAN *p;
	p = (EPSMEDIAN*) h;
	float t;

    for(i1=0; i1 < p->n; i1++)
    {
        min = i1;
        max = MY_MIN(i1+p->nfw, p->n-1);
        l = max - min + 1;
    	for(j1=min; j1 <= max; j1++)
    	p->u[j1-min] = x[j1*d];
        p->v[i1] = sfilt_median(l, p->u);
    	for(j1=min, p->e[i1]=0.0 ; j1 <= max; j1++)
		{
			t = p->v[i1] - p->u[j1-min];
			p->e[i1] += t*t;
		}
		p->e[i1] /= l;
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


