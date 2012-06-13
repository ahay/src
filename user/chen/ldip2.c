/* dip estimation by line-interpolating PWD */

#include <rsf.h>
#include "lphlag.h"
#include "lphbsp.h"
#include "lphmf.h"

static int n1, n2, nf, itp;
static float **c, **u1, **u2, **u3;

void ldip2_init(int mf, int interp,
	int m1, int m2, int niter, bool verb)
/*< initialize >*/
{
	float a;
	int i1, i2;

	nf = mf;
	n1 = m1;
	n2 = m2;
	itp = interp;

	u1 = sf_floatalloc2(n1, 2*nf+1);
	u2 = sf_floatalloc2(n1, 2*nf+1);
	u3 = sf_floatalloc2(n1, 2*nf+1);

	switch(itp)
	{
	case 1:
		c = lphlag(mf);
		break;
	case 2:
		c = lphbsp(mf);
		break;
	default:
		c = lphmf(mf);
	}
	for(i2=0; i2<2*nf+1; i2++)
	for(i1=0; i1<i2; i1++)
	{
		a = c[i2][i1];
		c[i2][i1] = c[i1][i2];
		c[i1][i2] = a;
	}
}



void ldip2_close()
/*< release memory >*/
{
	free(u1[0]);
	free(u1);
	free(u2[0]);
	free(u2);
	free(u3[0]);
	free(u3);
	free(c[0]);
	free(c);
}

static void fir(float *p, float*in, float *out)
{
	int i1, j1;
	for(i1=nf; i1<n1-nf; i1++)
	for(j1=-nf, out[i1] = 0.0; j1<=nf; j1++)
		out[i1] += p[j1+nf]*in[i1-j1];
}

#define divn(a, b)  (a*b/(b*b+0.0001))

void ldip2(float **in, float **dip, int nit)
/*< dip estimation >*/
{
	int it, i1, i2, i3, m_1, i4;
	double d1, d2, d3;

	for(i3=0; i3<2*nf+1; i3++)
		fir(c[i3], in[0], u1[i3]);
		

	for(i2=1; i2<n2; i2++)
	{
		for(i3=0, m_1=-1; i3<2*nf+1; i3++, m_1=-m_1)
		{
			fir(c[i3], in[i2], u2[i3]);
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
				for(i4=1; i4<i3; i4++)
				for(it=0; it<nit; it++)
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




