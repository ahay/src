/* dip estimation by line-interpolating PWD */

#include <rsf.h>
#include "lphmf.h"

static int n1, n2, nf;
static float **c, **u1, **u2;

void odip2_init(int mf, int interp,
	int m1, int m2, int niter, bool verb)
/*< initialize >*/
{
	float a;
	int i1, i2, n;

	nf = mf;
	n1 = m1;
	n2 = m2;

	n = 2*nf+1;
	u1 = sf_floatalloc2(n, n);
	u2 = sf_floatalloc2(n, n);

	c = lphmf(mf);
	for(i2=0; i2<2*nf+1; i2++)
	for(i1=0; i1<i2; i1++)
	{
		a = c[i2][i1];
		c[i2][i1] = c[i1][i2];
		c[i1][i2] = a;
	}
}



void odip2_close()
/*< release memory >*/
{
	free(u1[0]);
	free(u1);
	free(u2[0]);
	free(u2);
	free(c[0]);
	free(c);
}

static void fir2(float*in, float **out)
{
	int i1, i2, j1, j2, n;
	n = 2*nf + 1;

	for(i2=0; i2<n; i2++)
	for(i1=0; i1<n; i1++)
	for(j2=-nf, out[i2][i1]=0.0; j2<=nf; j2++)
	for(j1=-nf; j1<=nf; j1++)
		out[i2][i1] += in[j2*n1+j1]*c[i1][nf-j1]*c[i2][nf-j2];
}

static float odip2_lop(int nit)
{
	int i1, i2, n, b1, b2;

	n = 2*nf + 1;
	b1 = 1;
	b2 = 1;

	for(i2=0; i2<n; i2++, b2=-b2)
	for(i1=0; i1<n; i1++, b1=-b1)
		u2[i2][i1] = (1-b1*b2)*u1[i2][i1];
	return 0.0;
}

#define divn(a, b)  (a*b/(b*b+0.0001))

void odip2(float **in, float **dip, int nit)
/*< dip estimation >*/
{
	int i1, i2;

	for(i2=nf; i2<n2-nf; i2++)
	for(i1=nf; i1<n1-nf; i1++)
	{
		fir2(in[i2]+i1, u1);
		dip[i2][i1] = odip2_lop(nit);
	}
}


