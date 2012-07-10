/* omnidirectional dip estimation by PWD */

#include <rsf.h>
#include "opwd.h"

static int n1, n2, nf;
static float **u1, **u2, **u3;
static double rad;
static bool verb;

void odip_init(int interp, int mf, float r,
	int m1, int m2, int *rect, int niter, bool vb)
/*< initialize >*/
{
	int n, nn[2];
	nf = mf;
	n1 = m1;
	n2 = m2;
    n = n1*n2;
    nn[0] = n1;
    nn[1] = n2;
	verb = vb;
	rad=r;
	u1 = sf_floatalloc2(n1, n2);
	u2 = sf_floatalloc2(n1, n2);
	u3 = sf_floatalloc2(n1, n2);

	opwd_init(interp, nf, r);
	sf_divn_init (2, n, nn, rect, niter, false);
}

void odip_close()
/*< release memory >*/
{
	free(u1[0]);
	free(u2[0]);
	free(u3[0]);
	free(u1);
	free(u2);
	free(u3);
	opwd_close();
	sf_divn_close();
}

#define divn(a, b)  (a*b/(b*b+10E-15))


void odip(float ****in, float **dip, int nit, float eta)
/*< omnidirectional dip estimation >*/
{
	int it, i1;
	double  norm, s1, c1;


	for(i1=0; i1<n1*n2; i1++)
	{
		dip[0][i1] = SF_PI/4;
	}

	for (it=0; it<nit; it++)
	{
		opwd(n1, n2, in, dip, u1);
		opwdpd(n1, n2, in, dip, u2, 0);
		opwdpd(n1, n2, in, dip, u3, 1);

        if(verb)
        {
            for(i1=0, norm=0.0; i1<n1*n2; i1++)
                norm += (u1[0][i1]*u1[0][i1]);
            sf_warning("res1 %d %g", it+1, sqrtf(norm/n1/n2));
        }
/*

        for(i1=0, norm=0.0; i1<n1*n2; i1++)
            norm += (u2[0][i1]*u2[0][i1]);
        norm=sqrtf(norm/(n1*n2));
        for(i1=0; i1<n1*n2; i1++)
        {
            u1[0][i1] /= norm;
            u2[0][i1] /= norm;
        }
		sf_divn(u1[0], u2[0], u3[0]);
*/
		for(i1=0; i1<n1*n2; i1++)
		{
			s1=rad*sin(dip[0][i1]);
			s1 -= eta * divn(u1[0][i1], u2[0][i1]);
			c1=rad*cos(dip[0][i1]);
			c1 -= eta * divn(u1[0][i1], u3[0][i1]);
			dip[0][i1] = atan2(s1, c1);
		}

	}
/*	for(i1=0; i1<n1*n2; i1++)
	{
		while(dip[0][i1]>SF_PI/2) dip[0][i1] -= SF_PI/2;
		while(dip[0][i1]<-SF_PI/2) dip[0][i1] += SF_PI/2;
	}
*/
	for(i1=0; i1<n1*n2; i1++)
		dip[0][i1] = atan(tan(dip[0][i1]));
}



