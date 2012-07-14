/* omnidirectional dip estimation by PWD */

#include <rsf.h>
#include "opwd.h"

static int n1, n2, nf;
static float **u1, **u2, **u3, **u4, **u5, r;
static sf_complex **p, p0;
static bool verb, use_divn;

void odip_init(int interp, int mf, float rad,
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

	u1 = sf_floatalloc2(n1, n2);
	u2 = sf_floatalloc2(n1, n2);
	u3 = sf_floatalloc2(n1, n2);
	u4 = sf_floatalloc2(n1, n2);
	u5 = sf_floatalloc2(n1, n2);
	p = sf_complexalloc2(n1, n2);

	r=rad;
	p0 = rad*cexpf(sf_cmplx(0,SF_PI/4));
	
	opwd_init(interp, nf, r);
	if(rect[0]>1 && rect[1]>1)
	{
		sf_divn_init (2, n, nn, rect, niter, false);
		use_divn=true;
	}else 	use_divn=false;
}

void odip_close()
/*< release memory >*/
{
	free(u1[0]);
	free(u2[0]);
	free(u3[0]);
	free(p[0]);
	free(p);
	free(u1);
	free(u2);
	free(u3);
	opwd_close();
	sf_divn_close();
}

#define divn(a, b)  (a*b/(b*b+10E-15))


void odip(float **in, float **dip, int nit, float eta)
/*< omnidirectional dip estimation >*/
{
	int it, i1;
	double  norm, s1, c1;

	for(i1=0; i1<n1*n2; i1++)
		p[0][i1] = p0;

	for (it=0; it<nit; it++)
	{
		opwd(n1, n2, in, p, u1);
		opwdpd(n1, n2, in, p, u2, 0);
		opwdpd(n1, n2, in, p, u3, 1);

        if(verb)
        {
            for(i1=0, norm=0.0; i1<n1*n2; i1++)
                norm += (u1[0][i1]*u1[0][i1]);
            sf_warning("res1 %d %g", it+1, sqrtf(norm/n1/n2));
        }

		if(use_divn)
		{
	        for(i1=0, c1=0.0, s1=0.0; i1<n1*n2; i1++)
			{
    	        c1 += (u2[0][i1]*u2[0][i1]);
    	        s1 += (u3[0][i1]*u3[0][i1]);
			}
	        c1=sqrtf(c1/(n1*n2));
	        s1=sqrtf(s1/(n1*n2));
    	    for(i1=0; i1<n1*n2; i1++)
        	{
            	u1[0][i1] /= c1;
	            u2[0][i1] /= c1;
			}
			sf_divn(u1[0], u2[0], u4[0]);
    	    for(i1=0; i1<n1*n2; i1++)
        	{
            	u1[0][i1] *= c1/s1;
	            u3[0][i1] /= s1;
    	    }
			sf_divn(u1[0], u3[0], u5[0]);
		}else{
			for(i1=0; i1<n1*n2; i1++)
			{
				u4[0][i1] = divn(u1[0][i1], u2[0][i1]);
				u5[0][i1] = divn(u1[0][i1], u3[0][i1]);
			}
		}
		for(i1=0; i1<n1*n2; i1++)
		{
			p[0][i1] -= eta * sf_cmplx(u5[0][i1], u4[0][i1]);
			p[0][i1] = p[0][i1]*r/(cabsf(p[0][i1])+ 1E-15);
		}

	}
	for(i1=0; i1<n1*n2; i1++)
		dip[0][i1] = atan(tan(cargf(p[0][i1])));
}



