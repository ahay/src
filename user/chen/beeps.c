/* Bi-Exponential Edge-Preserving Smoothing */

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

#include "kernel.h"
/*^*/

typedef struct{
	int n;
	float *p, *q;
	kernel op;
	float sigma, lamda;
}BEEPS;

void* beeps_init(int n, kernel op, float sigma, float lamda)
/*< initialize >*/
{
	BEEPS *h;
	h = (BEEPS*)sf_alloc(1, sizeof(BEEPS));
	h->p = sf_floatalloc(n);
	h->q = sf_floatalloc(n);
	h->n = n;
	h->op = op;
	h->sigma = sigma;
	h->lamda = lamda;
	return h;
}

void beeps_close(void *p)
/*< release memory >*/
{
	BEEPS *h;

	h = (BEEPS*)p;
	free(h->p);
	free(h->q);
	free(h);
}

void beeps(void *p, float *x, int inc)
/*< beeps >*/
{
	int i, k;
	BEEPS *h;
	float r, xk;

	h = (BEEPS*)p;

	h->p[0] = x[0];
	h->q[h->n-1] = x[h->n-1];
	for(k=1; k<h->n; k++)
	{
		xk = x[k*inc];
		r = ( xk - h->p[k-1] ) / h->sigma;
		r = h->op(r) * h->lamda;
		h->p[k] =(1-r)*xk + r*h->p[k-1];

		i = h->n - 1 -k;
		xk = x[i*inc];
		r = (xk - h->q[i+1]) / h->sigma;
		r = h->op(r) * h->lamda;
		h->q[i] =(1-r)*xk + r*h->q[i+1];
	}

	for(k=0; k<h->n; k++)
	{
		i = k*inc;
		x[i] = (h->p[k] + h->q[k] - (1-h->lamda)*x[i]) / (1+h->lamda);
	}
}


