/* Linear phase filter by B-Spline interpolation */

#include <rsf.h>
#include "poly.h"


float **lphbsp(int n)
/*< Linear phase filter by noncausal B-Spline approximation
  [	0.25  1.0p   1.0p^2    Z
	1.5   0.0p  -2.0p^2	 
	0.25 -1.0p   1.0p^2    Z^{-1} ]
 >*/
{
	float **p, *a, *b;
	int i1, i2, j, m;

	m = 2*n;
	p = sf_floatalloc2(m+1, m+1);
	a = sf_floatalloc(m+1);
	b = sf_floatalloc(m+1);
	
	a[0] = 1.0;
	b[0] = 1.0;
	for(i2=1; i2<=m; i2++)
	{
		a[i2] =-a[i2-1] * (m+2-i2) / i2;
		b[i2] = b[i2-1] * (m+1-i2) / i2;
	}

	for(i1=-n; i1<=n; i1++)
	{
		for(i2=0; i2<=m; i2++)
		{
			p[n-i1][i2] = 0.0;
			for(j=0; j<=n-i1; j++)
				p[n-i1][i2] += a[j]*powf(n+0.5-i1-j, m-i2);
			p[n-i1][i2] *= b[i2];
		}
	}
	free(a);
	free(b);
	return p;
}


void lphbsp_filt(int n, float **c, float delay, float *out)
/*< filter >*/
{
	int i;

	for(i=0; i<=2*n; i++)
		out[i] = crealf(poly_val(2*n, c[i], delay));	
}


void lphbsp_dfilt(int n, float **c, float delay, float *out)
/*< derivative filter >*/
{
	int i;

	for(i=0; i<=2*n; i++)
		out[i] = crealf(poly_dval(2*n, c[i], delay));	
}


void lphbsp_freq(int n, float **c, float delay, int n1, sf_complex *out)
/*< frequency response >*/
{
	int i1;
	float *b1, *b2;
	sf_complex c0, c1;

	b1 = sf_floatalloc(2*n+1);
	b2 = sf_floatalloc(2*n+1);
	lphbsp_filt(n, c, delay, b1);

	for(i1=0; i1<=2*n; i1++) 
		b2[i1] = c[i1][0];

	for(i1=-n1; i1<=n1; i1++)
	{
		c0 = cexpf(sf_cmplx(0., SF_PI*i1*n/n1));
		c1 = cexpf(sf_cmplx(0., -SF_PI*i1/n1));
		out[i1+n1] = poly_val(2*n, b1, c1)/poly_val(2*n, b2, c1);
	}

	free(b1);
}

