/* vector compare modules */
#include <rsf.h>
//#include <math.h>

static int n;
static double *ref, pr;

void comp_init(int m, float *r)
/*< initialize >*/
{
	int i;

	n = m;

	ref = (double*) sf_alloc(n,sizeof(double));
	for(i=0, pr=0.0; i<n; i++)
	{
		pr += r[i]*r[i];
		ref[i] = r[i];
	}
}

float comp_xcor(float *d)
/*< unit x-correlation >*/
{
	double e, po;
	int i;

	po = 0.0;
	e  = 0.0;
	for(i=0; i<n; i++)
	{
		po += d[i]*d[i];
		e  += d[i]*ref[i];
	}
	e = 1.0 - e/sqrt(po*pr);

	return (float)( e );
}


void comp_close()
/*< free the space >*/
{
	free(ref);
}
