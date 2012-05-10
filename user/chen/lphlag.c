/* Linear phase filter by Lagrange approximation */

#include <rsf.h>
#include "poly.h"

typedef struct tag_lphlag
{
	int n;
	int n0;
	float **r;
}lphlag;

void * lphlag_init(int n, bool causal)
/*< initiliaze >*/
{
	lphlag *p;
	int m, m1, i1, i2;

	p = (lphlag*) sf_alloc(1, sizeof(lphlag));
	p->n = n;

	if(causal)	m1 = 0;
	else	m1 = -n;

	m = n -m1;
	p->r = sf_floatalloc2(m+1, m+1);
	for(i1=m1; i1<=n; i1++)
	{
		p->r[i1+m1][m] = 1.0;
		for(i2=m1; i2<i1; i2++)
		{
			p->r[i1+m1][m] /= (i1-i2);
			p->r[i1+m1][i2-m1] = i2;
		}
		for(i2=i1+1; i2<=n; i2++)
		{
			p->r[i1+m1][m] /= (i1-i2);
			p->r[i1+m1][i2-m1-1] = i2;
		}
	}
	p->n0 = m1;
	return p;
}

void lphlag_filt(void *h, float delay, float *out)
/*< filter >*/
{
	lphlag *p;
	int m, i;

	p = (lphlag*) h;

	m = p->n-p->n0;
	for(i=0; i<m; i++)
		out[i] = creal(plyr_val(m, p->r[i], delay));	
}


void lphlag_dfilt(void *h, float delay, float *out)
/*< derivative filter >*/
{
	lphlag *p;
	int m, i;

	p = (lphlag*) h;

	m = p->n-p->n0;
	for(i=0; i<m; i++)
		out[i] = creal(plyr_dval(m, p->r[i], delay));	
}


void lphlag_freq(void *h, float delay, int nk, sf_complex *out)
/*< frequency response >*/
{
	lphlag *p;
	int m, i;
	float *r;

	p = (lphlag*) h;
	m = p->n-p->n0;

	r = sf_floatalloc(m);
	lphlag_filt(h, delay, r);

	for(i=-nk; i<nk; i++)
		out[i+nk] = cexpf(sf_cmplx(0.,2*SF_PI*i*p->n0/nk)) *
			poly_val(m, r, cexpf(sf_cmplx(0.,2*SF_PI*i/nk)));

	free(r);
}

void lphlag_close(void *h)
/*< release memory >*/
{
	lphlag *p;
	p = (lphlag*) h;
	free(p->r[0]);
	free(p->r);
	free(p);
}

