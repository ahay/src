/* omnidirectional dip estimation by PWD */

#include <rsf.h>
#include "lpwd.h"

static int n1, n2, nf;
static float **u1, **u2, **u3;
static bool verb;

void ldip_init(int mf, int interp,
	int m1, int m2, int *rect, int niter, bool vb)
/*< initialize >*/
{
	int n, nn[4];
	nf = mf;
	n1 = m1;
	n2 = m2;
	verb=vb;

	n = n1*n2;
	nn[0] = n1;
	nn[1] = n2;

	u1 = sf_floatalloc2(n1,n2);
	u2 = sf_floatalloc2(n1,n2);
	u3 = sf_floatalloc2(n1,n2);

	lpwd_init(interp, nf, n1, n2);
	sf_divn_init (2, n, nn, rect, niter, false);
}

void ldip_close()
/*< release memory >*/
{
	free(u1[0]);
	free(u1);
	free(u2[0]);
	free(u2);
	free(u3[0]);
	free(u3);
	lpwd_close();
	sf_divn_close();
}

#define divn(a, b)  (a*b/(b*b+0.0001))

void ldip(float **in, float **dip, int nit)
/*< omnidirectional dip estimation >*/
{
	int it, i1;
	double eta, norm;

	for(i1=0; i1<n1*n2; i1++)
		dip[0][i1]=0.0;

	for (it=0; it<nit; it++)
	{
		lpwd(in, u1, dip);
		lpwdd(in, u2, dip);

		if(verb)
		{
			for(i1=0, norm=0.0; i1<n1*n2; i1++)
				norm += (u1[0][i1]*u1[0][i1]);
			sf_warning("%d %g", it+1, sqrtf(norm/n1/n2));
		}
		for(i1=0, norm=0.0; i1<n1*n2; i1++)
			norm += (u2[0][i1]*u2[0][i1]);
		norm=sqrtf(norm);
		for(i1=0; i1<n1*n2; i1++)
		{
			u1[0][i1] /= norm;
			u2[0][i1] /= norm;
		}

		sf_divn(u1[0], u2[0], u3[0]);

		for(i1=0; i1<n1*n2; i1++)
			dip[0][i1] -= u3[0][i1];
	}
//    for(i1=0; i1<n1*n2; i1++)
 //   {
   //     while(dip[0][i1]>SF_PI/2) dip[0][i1] -= SF_PI/2;
     //   while(dip[0][i1]<-SF_PI/2) dip[0][i1] += SF_PI/2;
//    }

	for(i1=0; i1<n1*n2; i1++)
		dip[0][i1] = atanf(dip[0][i1]);
}



