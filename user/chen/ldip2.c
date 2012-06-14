/* dip estimation by line-interpolating PWD */

#include <rsf.h>
#include "lphase.h"
#include "dsp.h"

static int n1, n2, nf, itp;
static float **c, **u1, **u2, **u3;

void ldip2_init(int mf, int interp,
	int m1, int m2, int niter, bool verb)
/*< initialize >*/
{
	int n;

	nf = mf;
	n1 = m1;
	n2 = m2;
	itp = interp;
	n = 2*nf+1;

	u1 = sf_floatalloc2(n1, n);
	u2 = sf_floatalloc2(n1, n);
	u3 = sf_floatalloc2(n1, n);

	c = lphase(nf, interp);
}



void ldip2_close()
/*< release memory >*/
{
	free(u1[0]);
	free(u2[0]);
	free(u3[0]);
	free(u1);
	free(u2);
	free(u3);
	free(c[0]);
	free(c);
}



#define divn(a, b)  (a*b/(b*b+0.0001))

void ldip2(float **in, float **dip, int nit)
/*< dip estimation >*/
{
	int it, i1, i2, i3, i4, m_1, n;
	double d1, d2, d3;

	n = 2*nf;

	for(i3=0; i3<=n; i3++)
		firs(-nf, nf, c[i3]+nf, in[0], 1, n1, u1[i3], 1);

	for(i2=1; i2<n2; i2++)
	{
		for(i3=0, m_1=-1; i3<=n; i3++, m_1=-m_1)
		{
			firs(-nf, nf, c[i3]+nf, in[i2], 1, n1, u2[i3], 1);
			for(i1=0; i1<n1; i1++)
			switch(itp)
			{
			case 1:
				u3[i3][i1] = (i3==0?in[i2-1][i1]:0.0) - u2[i3][i1];
				break;
			case 2:
				u3[i3][i1] = (i3==0?u1[i3][i1]:0.0) - u2[i3][i1];
				break;
			default:
				u3[i3][i1] = m_1*u1[i3][i1] + u2[i3][i1];
			}
		}
		for(i1=0; i1<n1; i1++)
		{
 			d1 = divn(-u3[0][i1], u3[1][i1]);
			for(i3=2; i3<2*nf+1; i3++)
			{
				d2 = u3[0][i1];
				d3 = 0.0;
				for(it=0; it<nit; it++)
				for(i4=1; i4<i3; i4++)
				{
					d2 += u3[i4][i1]*pow(d1, i4);
					d3 += u3[i4][i1]*pow(d1, i4-1)*i4;
					d1 += 0.5*divn(d2, d3);
				}
			}
			dip[i2][i1] = atan(d1);
		}
		memcpy(u1[0], u2[0], (2*nf+1)*n1*sizeof(float));
	}
}




