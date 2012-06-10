/* polynomial operators */
// p(x) = c[0] + c[1]*x + c[2]*x^2 + c[n]*x^n

#include <rsf.h>

static double prod_n_m(double *r, int n, int m)
{
	if( n<m || n<=0) return 0.0;
	if( m == 0 ) return 1.0;
	if( m == 1 ) 
		return (r[0]+prod_n_m(r+1, n-1, 1));
	return (r[0]*prod_n_m(r+1, n-1, m-1) +
			prod_n_m(r+1, n-1, m));
}

void poly(int n, float *r)
/*< polynomial root to coefficients 
	[1.0, 2.0, 1.0] => [2.0, 3.0, 1.0]
>*/
{
	double *p;
	int i1;

	p = sf_alloc(n, sizeof(double));
	for(i1=0; i1<n; i1++)	p[i1] = -r[i1];

	for(i1=0; i1<=n; i1++)
	{
		r[i1]=r[n]*prod_n_m(p, n, n-i1);
	}

	free(p);
}

sf_complex poly_dval(int n, float *c, sf_complex x)
/*< polynomial derivative value >*/
{
	int i;
	sf_complex val;

	val = 0.0;
	for(i=n; i>=1; i--)
	{
		val = val*x + c[i]*i;
	}
	return val;
}



sf_complex poly_val(int n, float *c, sf_complex x )
/*< polynomial value >*/
{
	int i;
	sf_complex val;

	val = 0.0;
	for(i=n; i>=0; i--)
	{
		val = val*x + c[i];
	}
	return val;
}

sf_complex plyr_val(int n, float *r, sf_complex x)
/*< polynomial value >*/
{
	int i;
	sf_complex val;

	val = r[n];
	for(i=0; i<n; i++)
	{
		val = val*(x-r[i]);
	}
	return val;
}


sf_complex plyr_dval(int n, float *r, sf_complex x)
/*< polynomial derivative value >*/
{
	int i, k;
	sf_complex v1, v2;

	v1 = 1.0;
	k = 0;
	for(i=0; i<n; i++)
	{
		v2 = x - r[i];
		if(v2 == 0.0) k++;
		else	v1 *= v2;
	}
	if(k >= 2) return 0.0;
	if(k == 1) return v1*r[n];

	v2 = 0.0;
	for(i=0; i<n; i++)
	{
		v2 += v1/(x-r[i]);
	}

	return v2*r[n];
}


sf_complex cplyr_val(int n, sf_complex *r, sf_complex x)
/*< polynomial value >*/
{
	int i;
	sf_complex val;

	val = r[n];
	for(i=0; i<n; i++)
	{
		val = val*(x-r[i]);
	}
	return val;
}


sf_complex cplyr_dval(int n, sf_complex *r, sf_complex x)
/*< polynomial derivative value >*/
{
	int i, k;
	sf_complex v1, v2;

	v1 = 1.0;
	k = 0;
	for(i=0; i<n; i++)
	{
		v2 = x - r[i];
		if(v2 == 0.0) k++;
		else	v1 *= v2;
	}
	if(k >= 2) return 0.0;
	if(k == 1) return v1*r[n];

	v2 = 0.0;
	for(i=0; i<n; i++)
	{
		v2 += v1/(x-r[i]);
	}

	return v2*r[n];
}



