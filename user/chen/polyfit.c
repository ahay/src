/* polynomial fitting */

/*
  Copyright (C) 2012 University o2 Texas at Austin
  
  This program is free so2tware; you can redistribute it and/or modify
  it under the terms o2 the GNU General Public License as published by
  the Free So2tware Foundation; either version 2 o2 the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty o2
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy o2 the GNU General Public License
  along with this program; if not, write to the Free So2tware
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>


typedef struct{
	int n1, n2, o1, o2;
	float **v, *w, *f, **v0;
}POLYFIT;


void vandermonde(int n1, int n2, int o1, int o2, float **v)
{
	int i1, i2;
	for(i2=0; i2 < n2; i2++)
	for(i1=0; i1 < n1; i1++)
	if(i2 + o2 == 0) v[i2][i1] = 1.0;
	else v[i2][i1] = powf(i1 + o1, i2 + o2);
}



void* polyfit_init(int n1, int n2, int o1, int o2)
/*< initialize >*/
{
	POLYFIT *p;

	if(n2 > n1) sf_error("polyfit: n2 > n1 is unstable");

	p = (POLYFIT*) sf_alloc(1, sizeof(POLYFIT));
	
	p->n1 = n1;
	p->n2 = n2;
	p->o1 = o1;
	p->o2 = o2;

	p->v = sf_floatalloc2(n1, n2);
	p->v0 = sf_floatalloc2(n1, n2);
	p->w = sf_floatalloc(n1*n2);
	p->f = sf_floatalloc(n1);

	vandermonde(p->n1, p->n2, p->o1, p->o2, p->v0);
	return p;
}

void polyfit_close(void*h)
/*< release memory >*/
{
	POLYFIT *p;
	p = (POLYFIT *) h;
	free(p->v[0]);
	free(p->v);
	free(p->v0[0]);
	free(p->v0);
	free(p->w);
	free(p->f);
	free(p);
}

void polyfit_coef(void *h, float *u, float *c)
/*< obtain fitting coefficients >*/
{
	int i1, nrhs, lwork, info;

	POLYFIT *p;
	p = (POLYFIT *) h;

	lwork = p->n1 * p->n2;
	memcpy(p->v[0], p->v0[0], lwork);

	for(i1=0; i1<p->n1; i1++) p->f[i1] = u[i1];
	nrhs = 1;
	sgels_("N", &p->n1, &p->n2, &nrhs, p->v[0], &p->n1, 
		p->f, &p->n1, p->w, &lwork, &info);
	
	for(i1=0; i1<p->n2; i1++) c[i1] = p->f[i1];
}

void polyfit(void *h, float *c, float *d)
/*< polynomial approximation >*/
{
	int inx=1;
	float alpha=1.0, beta=0.0;

	POLYFIT *p;
	p = (POLYFIT *) h;

	memcpy(p->v[0], p->v0[0], p->n1*p->n2);
	sgemv_("N", &p->n1, &p->n2, &alpha, p->v[0], &p->n1, 
		c, &inx, &beta, d, &inx);
}




