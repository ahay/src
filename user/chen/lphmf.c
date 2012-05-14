/* Linear phase filter by maxflat approximation */

#include <rsf.h>
#include "poly.h"

typedef struct tag_lphmf
{
	int n;
	int n0;
	float **r;
}lphmf;

void * lphmf_init(int n, bool causal)
/*< initiliaze >*/
{
	lphmf *p;
	int m, i1, i2;

	p = (lphmf*) sf_alloc(1, sizeof(lphmf));
	p->n = n;

	if(causal)
	{
		p->n0 = 0;
		m = n;
		p->r = sf_floatalloc2(m+1, m+1);
		for(i1=0; i1<=n; i1++)
		{
			p->r[i1][m] = ((i1%2==0)? 1.0:-1.0);
			for(i2=0; i2<i1; i2++)
			{
				p->r[i1][m] /= (i1+1);
				p->r[i1][i2] = -i2;
			}
			for(i2=i1+n+1; i2<=2*n; i2++)
			{
				p->r[i1][m] /= (i2-i1-n);
				p->r[i1][i2-n-1] = -i2;
			}
		}
	}else{
		p->n0 = -n;
		m = 2*n;
		p->r = sf_floatalloc2(m+1, m+1);
		for(i1=-n; i1<=n; i1++)
		{
			p->r[i1+n][m] = ((n+i1)%2==0?1.0:-1.0);
			for(i2=0; i2<i1+n; i2++)
			{
				p->r[i1+n][m] /= (i2+1);
				p->r[i1+n][i2] = -(i2-m);
			}
			for(i2=i1+n+1; i2<=m; i2++)
			{
				p->r[i1+n][m] /= (i2-i1-n);
				p->r[i1+n][i2-1] = -i2;
			}
		}
	}
	return p;
}

void lphmf_filt(void *h, float delay, float *out)
/*< filter >*/
{
	lphmf *p;
	int m, i;

	p = (lphmf*) h;

	m = p->n-p->n0;
	for(i=0; i<=m; i++)
		out[i] = creal(plyr_val(m, p->r[i], delay));	
}


void lphmf_dfilt(void *h, float delay, float *out)
/*< derivative filter >*/
{
	lphmf *p;
	int m, i;

	p = (lphmf*) h;

	m = p->n-p->n0;
	for(i=0; i<=m; i++)
		out[i] = creal(plyr_dval(m, p->r[i], delay));	
}


void lphmf_freq(void *h, float delay, int nk, sf_complex *out)
/*< frequency response >*/
{
	lphmf *p;
	int m, i;
	float *r;
	sf_complex c1, c2;

	p = (lphmf*) h;
	m = p->n-p->n0;

	r = sf_floatalloc(m+1);
	lphmf_filt(h, delay, r);
	
	for(i=-nk; i<=nk; i++)
	{
		c1 = cexpf(sf_cmplx(0., -2.0*SF_PI*i*p->n0/nk)) *
		poly_val(m, r, cexpf(sf_cmplx(0.,-2.0*SF_PI*i/nk)));
		c2 = cexpf(sf_cmplx(0., 2.0*SF_PI*i*p->n0/nk)) *
		poly_val(m, r, cexpf(sf_cmplx(0., 2.0*SF_PI*i/nk)));
		out[i+nk] = c1/c2;
	}


	free(r);
}

void lphmf_close(void *h)
/*< release memory >*/
{
	lphmf *p;
	p = (lphmf*) h;
	free(p->r[0]);
	free(p->r);
	free(p);
}

