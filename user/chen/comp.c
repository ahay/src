/* vector compare modules */
#include <rsf.h>
//#include <math.h>

float comp_ncc(float *a, float *b, int n)
/*< normalized cross-correlation >*/
{
	double e, pa, pb;
	int i;

	pa = 0.0; 
	pb = 0.0;
	e  = 0.0;
	for(i=0; i<n; i++)
	{
		pa += a[i]*a[i];
		pb += b[i]*b[i];
		e  += a[i]*b[i];
	}
	e = e/sqrt(pa*pb);

	return (float)( e );
}



float comp_cor(float *a, float *b, int n)
/*< cross-correlation >*/
{
	double e, pa, pb;
	int i;

	pa = 0.0; 
	pb = 0.0;
	e  = 0.0;
	for(i=0; i<n; i++)
	{
		pa += a[i]*a[i];
		pb += b[i]*b[i];
		e  += 2*a[i]*b[i];
	}
	e = e/(pa+pb);

	return (float)( e );
}

float comp_mae(float *a, float *b, int n)
/*< mean absolute error >*/
{
	double e;
	int i;
	e=0.0;
	for(i=0; i<n; i++)
		e += fabs(a[i]-b[i]);
	return (float)(e/n);
}

float comp_mse(float *a, float *b, int n)
/*< mean square error >*/
{
	double e,t;
	int i;
	e=0.0;
	for(i=0; i<n; i++)
	{
		t = a[i]-b[i];
		e += t*t;
	}
	return (float)(e/n);
}


