/* omnidirectional dip estimation by PWD */

#include <rsf.h>
#include "opwd.h"

static int n1, n2, nf;
static float **u1, **u2;
static double rad;

void odip2_init(float r,  int mf, int interp,
	int m1, int m2, int niter, bool verb)
/*< initialize >*/
{
	nf = mf;
	n1 = m1;
	n2 = m2;
	rad=r;
	u1 = sf_floatalloc2(n1,n2);
	u2 = sf_floatalloc2(n1,n2);

	opwd_init(interp, nf, n1, n2, r);
}

void odip2_close()
/*< release memory >*/
{
	free(u1[0]);
	free(u2[0]);
	free(u1);
	free(u2);
	opwd_close();
}

#define divn(a, b)  (a*b/(b*b+0.0001))


void odip2(float **in, float **dip, int nit)
/*< omnidirectional dip estimation >*/
{
	int it, i1;
	double eta, u3, s1, c1;

	for(i1=0; i1<n1*n2; i1++)
	{
		dip[0][i1] = 0.0;
	}

	for (it=0; it<nit; it++)
	{
		eta=1.0/(1.0+it*it);	
		opwd(in, u1, dip);
		opwdpd(in, u2, dip, 0);
		for(i1=0; i1<n1*n2; i1++)
		{
			u3 = divn(u1[0][i1], u2[0][i1]);
			s1=rad*sin(dip[0][i1]) - u3;
			c1=rad*cos(dip[0][i1]);
			dip[0][i1] = atan2(s1*c1, c1*c1+0.00001);
		}

		opwd(in, u1, dip);
		opwdpd(in, u2, dip, 1);
		for(i1=0; i1<n1*n2; i1++)
		{
			u3 = divn(u1[0][i1], u2[0][i1]);
			s1=rad*sin(dip[0][i1]);
			c1=rad*cos(dip[0][i1]) - u3;
			dip[0][i1] = atan2(s1*c1, c1*c1+0.00001);
		}
	}

	for(i1=0; i1<n1*n2; i1++)
	{
		while(dip[0][i1]>SF_PI/2) dip[0][i1] -= SF_PI/2;
		while(dip[0][i1]<-SF_PI/2) dip[0][i1] += SF_PI/2;
	}
}



