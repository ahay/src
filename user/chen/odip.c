/* omnidirectional dip estimation by PWD */

#include <rsf.h>
#include "lpwd.h"
#include "opwd.h"

static int n1, n2, nf;
static float **u1, **u2, **u3;
static bool od;

void odip_init(bool opwd, float r,  int mf,
	int m1, int m2, int *rect, int niter, bool verb)
/*< initialize >*/
{
	int n, nn[4];
	nf = mf;
	n1 = m1;
	n2 = m2;
	od = opwd;

	n = n1*n2;
	nn[0] = n1;
	nn[1] = n2;

	u1 = sf_floatalloc2(n1,n2);
	u2 = sf_floatalloc2(n1,n2);
	u3 = sf_floatalloc2(n1,n2);

	if(od)	opwd_init(0, nf, n1, n2, r);
	else lpwd_init(0, nf, n1, n2);
	sf_divn_init (2, n, nn, rect, niter, verb);
}

void odip_close()
/*< release memory >*/
{
	free(u1[0]);
	free(u1);
	free(u2[0]);
	free(u2);
	free(u3[0]);
	free(u3);
	if(od) opwd_close();
	else lpwd_close();
	sf_divn_close();
}

void odip(float **in, float **dip, int nit)
/*< omnidirectional dip estimation >*/
{
	int it, i1;
	double norm, eta;

	for(i1=0; i1<n1*n2; i1++) 
		dip[0][i1] = 0.0; // SF_PI_2*rand()/(RAND_MAX+1.0);

	for (it=0; it<nit; it++)
	{	
		eta = 1.0/(1.0+it*it);
		if(od)
		{
			opwd(in, u1, dip, false);
			opwd(in, u2, dip, true);
		}else{
			lpwd(in, u1, dip, false);
			lpwd(in, u2, dip, true);
		}

/*		norm = 0.0;
		for(i1=0; i1<n1*n2; i1++) 
			norm += u2[0][i1]*u2[0][i1];
		if(norm == 0.0) return;
		norm = sqrtf(norm/n1/n2);

		for(i1=0; i1<n1*n2; i1++) 
		{
			u1[0][i1] /= norm;
			u2[0][i1] /= norm;
		}

		sf_divn(u1[0], u2[0], u3[0]);
*/		for(i1=0; i1<n1*n2; i1++)
		{
			u3[0][i1]=u1[0][i1]*u2[0][i1]/
				(u2[0][i1]*u2[0][i1]+0.0001);
			dip[0][i1] -= eta*u3[0][i1];
			while(dip[0][i1]>SF_PI/2) dip[0][i1] -= SF_PI/2;
			while(dip[0][i1]<-SF_PI/2) dip[0][i1] += SF_PI/2;
		}
	}
}



