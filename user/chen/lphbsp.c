/* Linear phase filter by B-Spline interpolation */

#include <rsf.h>
#include "poly.h"

typedef struct tag_lphbsp
{
	int n;
	int n0;
	float **r;
}lphbsp;

void * lphbsp_init(int n, bool causal)
/*< initiliaze >*/
{
	lphbsp *p;
	int m, i1, i2, j;
	float *a, *b;

	if(!causal) sf_warning("Noncausal B-spline not support");

	p = (lphbsp*) sf_alloc(1, sizeof(lphbsp));
	p->n = n;
	p->n0 = -n;

	m = 2*n;
	p->r = sf_floatalloc2(m+1, m+1);
	a = sf_floatalloc(m+1);
	b = sf_floatalloc(m+1);
	
	a[0] = 1.0;
	b[0] = 1.0;
	for(i2=1; i2<=2*n; i2++)
	{
		a[i2] =-a[i2-1] * (m+2-i2) / i2;
		b[i2] = b[i2-1] * (m+1-i2) / i2;
	}

	for(i1=-n; i1<=n; i1++)
	{
		for(i2=0; i2<=m; i2++)
		{
			p->r[n-i1][i2] = 0.0;
			for(j=0; j<=n-i1; j++)
				p->r[n-i1][i2] += a[j]*powf(n+0.5-i1-j, m-i2);
			p->r[n-i1][i2] *= b[i2];
		}
	}
	free(a);
	free(b);
	return p;
}

void lphbsp_filt(void *h, float delay, float *out)
/*< filter >*/
{
	lphbsp *p;
	int m, i;

	p = (lphbsp*) h;

	m = p->n-p->n0;
	for(i=0; i<=m; i++)
		out[i] = creal(poly_val(m, p->r[i], delay));	
}


void lphbsp_dfilt(void *h, float delay, float *out)
/*< derivative filter >*/
{
	lphbsp *p;
	int m, i;

	p = (lphbsp*) h;

	m = p->n-p->n0;
	for(i=0; i<=m; i++)
		out[i] = creal(poly_dval(m, p->r[i], delay));	
}


void lphbsp_freq(void *h, float delay, int nk, sf_complex *out)
/*< frequency response >*/
{
	lphbsp *p;
	int m, i;
	float *r;
	sf_complex c1, c2;

	p = (lphbsp*) h;
	m = p->n-p->n0;

	r = sf_floatalloc(m+1);
	lphbsp_filt(h, delay, r);
	
	for(i=-nk; i<=nk; i++)
	{
		c1 = cexpf(sf_cmplx(0., -2.0*SF_PI*i*p->n0/nk)) *
		poly_val(m, r, cexpf(sf_cmplx(0.,-2.0*SF_PI*i/nk)));
		out[i+nk] = c1;
	}


	free(r);
}

void lphbsp_close(void *h)
/*< release memory >*/
{
	lphbsp *p;
	p = (lphbsp*) h;
	free(p->r[0]);
	free(p->r);
	free(p);
}

