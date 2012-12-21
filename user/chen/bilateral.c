/* Bilateral filter */

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


typedef struct{
	int n1, n2;
	float **u, *v;
}BILATERAL;

void* bilateral_init(int n1, int n2)
/*< initialize >*/
{
	BILATERAL *h;
	int n, i, j;

	h = (BILATERAL*)sf_alloc(1, sizeof(BILATERAL));
	n = n1*n2;
	h->u = sf_floatalloc2(n, n);
	h->v = sf_floatalloc(n);
	h->n1 = n1;
	h->n2 = n2;

	for(i=0; i<n; i++)
	for(j=0; j<n; j++)	h->u[i][j] = 1.0;
	return h;
}

void bilateral_close(void *p)
/*< release memory >*/
{
	BILATERAL *h;

	h = (BILATERAL*)p;
	free(h->u[0]);
	free(h->u);
	free(h->v);
	free(h);
}

void bilateral_norm(void *p)
/*< normalization >*/
{
	int i, j, n;
	float s;
	BILATERAL *h;
	h = (BILATERAL*)p;

	n = h->n1*h->n2;
	for(j=0; j < n; j++)
	{
		for(i=0, s=0.0; i < n; i++) s += h->u[i][j];
		for(i=0; i < n; i++) h->u[i][j] /= s;
	}
}

void bilateral_df(void *p, float sigma)
/*< domain filter only >*/
{
	BILATERAL *h;
	int i1, j1, i2, j2, i, j;
	float x1, x2, s1;

	h = (BILATERAL*)p;
	s1 = -2.0*sigma*sigma;
	for(j2=0; j2 < h->n2; j2++)
	for(j1=0; j1 < h->n1; j1++)
	for(i2=0, j = j2*h->n1+j1; i2 <= j2; i2++)
	for(i1=0; i1 <= j1; i1++)
	{
		i = i2*h->n1+i1;
		x1 = j1-i1;
		x2 = j2-i2;
		h->u[i][j] *= expf((x1*x1+x2*x2)/s1);
		if(i != j) h->u[j][i] *= h->u[i][j];
	}
}

void bilateral_rf(void *p, float sigma, float *x)
/*< range filter only >*/
{
	BILATERAL *h;
	int i1, j1, i2, j2, i, j;
	float x1, s2;

	h = (BILATERAL*)p;
	s2 = -2.0*sigma*sigma;

	for(j2=0; j2 < h->n2; j2++)
	for(j1=0; j1 < h->n1; j1++)
	for(i2=0, j = j2*h->n1+j1; i2 <= j2; i2++)
	for(i1=0; i1 <= j1; i1++)
	{
		i = i2*h->n1+i1;
		x1 = x[i]-x[j];
		h->u[i][j] *= expf((x1*x1)/s2);
		if(i != j) h->u[j][i] *= h->u[i][j];
	}
}

void bilateral_ar(void *p, float sigma, float *x)
/*< adaptive range filter >*/
{
	BILATERAL *h;
	int i1, j1, i2, j2, i, j;
	float x1, s2;

	h = (BILATERAL*)p;
	s2 = -2.0*sigma*sigma;

	for(j2=0; j2 < h->n2; j2++)
	for(j1=0; j1 < h->n1; j1++)
	for(i2=0, j = j2*h->n1+j1; i2 <= j2; i2++)
	for(i1=0; i1 <= j1; i1++)
	{
		i = i2*h->n1+i1;
		x1 = x[i]-x[j];
		h->u[i][j] *= expf((x1*x1)/s2);
		if( i != j ) h->u[j][i] *= h->u[i][j];
	}
}


void bilateral(void *p, float *x)
/*< apply the filter >*/
{
	BILATERAL *h;
	int n, i, j;

	h = (BILATERAL*)p;

	n = h->n1 * h->n2;
	for(i=0; i<n; i++) h->v[i] = x[i];

	for(i=0; i<n; i++)
	for(j=0, x[i]=0.0; j<n; j++)
	x[i] += h->u[i][j] * h->v[j];
}



