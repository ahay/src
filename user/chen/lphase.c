/* Linear phase filters */

#include <rsf.h>
#include "poly.h"

typedef struct tag_lphase
{
	int n;
	int n0;
	int itp;	// 0: maxflat 	1: Lagrange 2: bspline
	float **r;
}lphase;


static void lphase_maxflat(int n, bool causal, float **r)
{
	int m, i1, i2;

	if(causal)
	{
		m = n;
		for(i1=0; i1<=n; i1++)
		{
			r[i1][m] = ((i1%2==0)? 1.0:-1.0);
			for(i2=0; i2<i1; i2++)
			{
				r[i1][m] /= (i1+1);
				r[i1][i2] = -i2;
			}
			for(i2=i1+n+1; i2<=2*n; i2++)
			{
				r[i1][m] /= (i2-i1-n);
				r[i1][i2-n-1] = -i2;
			}
		}
	}else{
		m = 2*n;
		for(i1=-n; i1<=n; i1++)
		{
			r[i1+n][m] = ((n+i1)%2==0?1.0:-1.0);
			for(i2=0; i2<i1+n; i2++)
			{
				r[i1+n][m] /= (i2+1);
				r[i1+n][i2] = -(i2-m);
			}
			for(i2=i1+n+1; i2<=m; i2++)
			{
				r[i1+n][m] /= (i2-i1-n);
				r[i1+n][i2-1] = -i2;
			}
		}
	}
}


static void lphase_lagrange(int n, bool causal, float **r)
{
	int m, m1, i1, i2;

	if(causal)  m1 = 0;
	else	m1 = -n;

	m = n -m1;
	for(i1=m1; i1<=n; i1++)
	{
		r[n-i1][m] = 1.0;
		for(i2=m1; i2<i1; i2++)
		{
			r[n-i1][m] /= (i1-i2);
			r[n-i1][i2-m1] = i2;
		}
		for(i2=i1+1; i2<=n; i2++)
		{
			r[n-i1][m] /= (i1-i2);
			r[n-i1][i2-m1-1] = i2;
		}
	}
}


void * lphase_init(int n, bool causal, int interp)
/*< initiliaze >*/
{
	lphase *p;
	
	p = (lphase*) sf_alloc(1, sizeof(lphase));
	p->n = n;
	p->itp = interp;
	if(causal)
	{
		p->n0 = 0;
		p->r = sf_floatalloc2(n+1, n+1);
	}else{
		p->n0 = -n;
		p->r = sf_floatalloc2(2*n+1, 2*n+1);
	}
	switch(interp)
	{
	case 1:
		lphase_lagrange(n, causal, p->r);
		break;
	default:
		lphase_maxflat(n, causal, p->r);
		break;
	}
	return p;
}

void lphase_filt(void *h, float delay, float *out, bool der)
/*< filter >*/
{
	lphase *p;
	int m, i;

	p = (lphase*) h;

	m = p->n-p->n0;
	if(der)
	for(i=0; i<m; i++)
		out[i] = creal(plyr_dval(m, p->r[i], delay));	
	else
	for(i=0; i<m; i++)
		out[i] = creal(plyr_val(m, p->r[i], delay));	
}


void lphase_freq(void *h, float delay, int nk, sf_complex *out)
/*< frequency response >*/
{
	lphase *p;
	int m, i;
	float *r;
	sf_complex c1, c2;

	p = (lphase*) h;
	m = p->n-p->n0;

	r = sf_floatalloc(m+1);
	lphase_filt(h, delay, r, false);

	switch(p->itp)
	{
	case 1:
		for(i=-nk; i<nk; i++)
		{
			out[i+nk] = cexpf(sf_cmplx(0., -2.0*SF_PI*i*p->n0/nk)) *
			poly_val(m, r, cexpf(sf_cmplx(0.,-2.0*SF_PI*i/nk)));
		}
		break;
	default:
		for(i=-nk; i<nk; i++)
		{
			c1 = cexpf(sf_cmplx(0., -2.0*SF_PI*i*p->n0/nk)) *
			poly_val(m, r, cexpf(sf_cmplx(0.,-2.0*SF_PI*i/nk)));
			c2 = cexpf(sf_cmplx(0., 2.0*SF_PI*i*p->n0/nk)) *
			poly_val(m, r, cexpf(sf_cmplx(0., 2.0*SF_PI*i/nk)));
			out[i+nk] = c1/c2;
		}
		break;
	}

	free(r);
}

void lphase_close(void *h)
/*< release memory >*/
{
	lphase *p;
	p = (lphase*) h;
	free(p->r[0]);
	free(p->r);
	free(p);
}

