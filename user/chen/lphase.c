/* Linear phase filters */

#include <rsf.h>
#include "poly.h"

float **lphase_lag(int n)
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

float **lphase_mf(int n)
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

float **lphase_bsp(int n)
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



