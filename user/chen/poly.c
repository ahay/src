/* polynomial operators */
// p(x) = c[0] + c[1]*x + c[2]*x^2 + c[n]*x^n

#include <rsf.h>

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



