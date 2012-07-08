/* omnidirectional dip estimation by PWD */

#include <rsf.h>
#include "opwd.h"

static int n1, n2, nf;
static float **u1, **u2, **u3, ****m;
static double rad;
static bool verb;

void odip_init(float r,  int mf, int interp,
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
	m = sf_floatalloc4(2*nf+1, 2*nf+1, n1, n2); // m[i2][i1][j2][j1]

	opwd_init(interp, nf, r);
	sf_divn_init (2, n, nn, rect, niter, verb);
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
	free(m[0][0][0]);
	free(m[0][0]);
	free(m[0]);
	free(m);
	opwd_close();
	sf_divn_close();
}

#define divn(a, b)  (a*b/(b*b+10E-15))


void odip(float **in, float **dip, int nit)
/*< omnidirectional dip estimation >*/
{
	int it, i1;
	double eta, norm, s1, c1;

	for(i1=0; i1<n1*n2; i1++)
	{
		dip[0][i1] = 0.0;
	}
	opwd_fbank(n1, n2, in, m);

	for (it=0; it<nit; it++)
	{
		eta=1.0/(1.0+it*it);	
		opwd(n1, n2, m, dip, u1);
		opwdpd(n1, n2, m, dip, u2, 0);

        if(verb)
        {
            for(i1=0, norm=0.0; i1<n1*n2; i1++)
                norm += (u1[0][i1]*u1[0][i1]);
            sf_warning("aaa %d %g", it+1, sqrtf(norm/n1/n2));
        }


        for(i1=0, norm=0.0; i1<n1*n2; i1++)
            norm += (u2[0][i1]*u2[0][i1]);
        norm=sqrtf(norm/(n1*n2));
        sf_warning("bbb %d %g", it+1, sqrtf(norm/n1/n2));
        for(i1=0; i1<n1*n2; i1++)
        {
            u1[0][i1] /= norm;
            u2[0][i1] /= norm;
        }
		sf_divn(u1[0], u2[0], u3[0]);

		for(i1=0; i1<n1*n2; i1++)
		{
			s1=rad*sin(dip[0][i1]) - u3[0][i1];
			c1=rad*cos(dip[0][i1]);
			dip[0][i1] = atan2(s1, c1);
		}
		opwd(n1, n2, m, dip, u1);
		opwdpd(n1, n2, m, dip, u2, 1);

        if(verb)
        {
            for(i1=0, norm=0.0; i1<n1*n2; i1++)
                norm += (u1[0][i1]*u1[0][i1]);
            sf_warning("aaa %d %g", it+1, sqrtf(norm/n1/n2));
        }


        for(i1=0, norm=0.0; i1<n1*n2; i1++)
            norm += (u2[0][i1]*u2[0][i1]);
        norm=sqrtf(norm/(n1*n2));
        for(i1=0; i1<n1*n2; i1++)
        {
            u1[0][i1] /= norm;
            u2[0][i1] /= norm;
        }
		sf_divn(u1[0], u2[0], u3[0]);
		for(i1=0; i1<n1*n2; i1++)
		{
		//	u3 = divn(u1[0][i1], u2[0][i1]);
			s1=rad*sin(dip[0][i1]);
			c1=rad*cos(dip[0][i1]) - u3[0][i1];
			dip[0][i1] = atan2(s1, c1);
		}

	}

	for(i1=0; i1<n1*n2; i1++)
	{
		while(dip[0][i1]>SF_PI/2) dip[0][i1] -= SF_PI/2;
		while(dip[0][i1]<-SF_PI/2) dip[0][i1] += SF_PI/2;
	}
}



