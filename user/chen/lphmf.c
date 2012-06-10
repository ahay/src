/* Linear phase filter by maxflat approximation */

#include <rsf.h>
#include "poly.h"


float **lphmf(int n)
/*< Linear phase filter by noncausal maxflat approximation 
  [	1.0   1.5p   0.5p^2    Z
	4.0   0.0p  -1.0p^2	 
	1.0  -1.5p   0.5p^2    Z^{-1} ]
 >*/
{
	float **p;
	int i1, i2, m;

	m = 2*n;
	p = sf_floatalloc2(m+1, m+1);
	for(i1=-n; i1<=n; i1++)
	{
		p[i1+n][m] = ((n+i1)%2==0?1.0:-1.0);
		for(i2=0; i2<i1+n; i2++)
		{
			p[i1+n][m] /= (i2+1);
			p[i1+n][i2] = -(i2-m);
		}
		for(i2=i1+n+1; i2<=m; i2++)
		{
			p[i1+n][m] /= (i2-i1-n);
			p[i1+n][i2-1] = -i2;
		}
		poly(m, p[i1+n]);
	}
	return p;
}


void lphmf_filt(int n, float **c, float delay, float *out)
/*< filter >*/
{
	int i;

	for(i=0; i<=2*n; i++)
		out[i] = crealf(poly_val(2*n, c[i], delay));	
}


void lphmf_dfilt(int n, float **c, float delay, float *out)
/*< derivative filter >*/
{
	int i;

	for(i=0; i<=2*n; i++)
		out[i] = crealf(poly_dval(2*n, c[i], delay));	
}


void lphmf_freq(int n, float **c, float delay, int n1, sf_complex *out)
/*< frequency response >*/
{
	int i1;
	float *b1;
	sf_complex c0, c1, c2;

	b1 = sf_floatalloc(2*n+1);
	lphmf_filt(n, c, delay, b1);

	for(i1=-n1; i1<=n1; i1++)
	{
		c0 = cexpf(sf_cmplx(0., SF_PI*i1*2.0*n/n1));
		c1 = cexpf(sf_cmplx(0., -SF_PI*i1/n1));
		c2 = poly_val(2*n, b1, c1);
		out[i1+n1] = c0*c2/conj(c2);
	}

	free(b1);
}



