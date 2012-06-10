/* Linear phase filter by Lagrange approximation */

#include <rsf.h>
#include "poly.h"


float **lphlag(int n)
/*< Linear phase filter by noncausal Lagrange approximation
  [	0.0   0.5p   0.5p^2    Z
	1.0   0.0p  -1.0p^2	 
	0.0  -0.5p   0.5p^2    Z^{-1} ]
 >*/
{
	float **p;
	int i1, i2, m;

	m = 2*n;
	p = sf_floatalloc2(m+1, m+1);
	for(i1=-n; i1<=n; i1++)
	{
		p[n-i1][m] = 1.0;
		for(i2=-n; i2<i1; i2++)
		{
			p[n-i1][m] /= (i1-i2);
			p[n-i1][i2+n] = i2;
		}
		for(i2=i1+1; i2<=n; i2++)
		{
			p[n-i1][m] /= (i1-i2);
			p[n-i1][i2+n-1] = i2;
		}
		poly(m, p[n-i1]);
	}
	return p;
}


void lphlag_filt(int n, float **c, float delay, float *out)
/*< filter >*/
{
	int i;

	for(i=0; i<=2*n; i++)
		out[i] = crealf(poly_val(2*n, c[i], delay));	
}


void lphlag_dfilt(int n, float **c, float delay, float *out)
/*< derivative filter >*/
{
	int i;

	for(i=0; i<=2*n; i++)
		out[i] = crealf(poly_dval(2*n, c[i], delay));	
}


void lphlag_freq(int n, float **c, float delay, int n1, sf_complex *out)
/*< frequency response >*/
{
	int i1;
	float *b1;
	sf_complex c0, c1;

	b1 = sf_floatalloc(2*n+1);
	lphlag_filt(n, c, delay, b1);

	for(i1=-n1; i1<=n1; i1++)
	{
		c0 = cexpf(sf_cmplx(0., SF_PI*i1*n/n1));
		c1 = cexpf(sf_cmplx(0., -SF_PI*i1/n1));
		out[i1+n1] = c0*poly_val(2*n, b1, c1);
	}

	free(b1);
}



