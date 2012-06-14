/* omnidirectional dip estimation by PWD */

#include <rsf.h>
#include "opwd.h"

static int n1, n2, nf;
static float **u1, **u2, **d1, **d2;

void odip1_init(float r,  int mf, int interp,
	int m1, int m2, int niter, bool verb)
/*< initialize >*/
{
	int n, nn[4];
	nf = mf;
	n1 = m1;
	n2 = m2;

	n = n1*n2;
	nn[0] = n1;
	nn[1] = n2;

	u1 = sf_floatalloc2(n1,n2);
	u2 = sf_floatalloc2(n1,n2);
	d1 = sf_floatalloc2(n1,n2);
	d2 = sf_floatalloc2(n1,n2);

	opwd_init(interp, nf, n1, n2, r);
}

void odip1_close()
/*< release memory >*/
{
	free(u1[0]);
	free(u1);
	free(u2[0]);
	free(u2);
	free(d1[0]);
	free(d1);
	free(d2[0]);
	free(d2);
	opwd_close();
}



void odip1(float **in, float **dip, int nit)
/*< omnidirectional dip estimation >*/
{
	int it, i1;
	double eta, u3;

	for(i1=0; i1<n1*n2; i1++)
	{
		dip[0][i1]=0.0;
	}

	for (it=0; it<nit; it++)
	{
		eta=1.0/(1.0+it*it);	
		opwd(in, u1, dip);
		opwdd(in, u2, dip);

		for(i1=0; i1<n1*n2; i1++)
		{
			u3=u1[0][i1]*u2[0][i1]/
				(u2[0][i1]*u2[0][i1]+0.0001);
			dip[0][i1] -= eta*u3;
			while(dip[0][i1]>SF_PI/2) dip[0][i1] -= SF_PI/2;
			while(dip[0][i1]<-SF_PI/2) dip[0][i1] += SF_PI/2;
		}
	}
}



