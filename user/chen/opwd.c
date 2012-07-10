/* circle-interpolating pwd*/

#include <rsf.h>
#include "lphase.h"
#include "dsp.h"

static int nf, itp;
static float *b1, *b2, *d1, *d2;
static float r, **c;

void opwd_init(int interp, int mf, float rad)
/*< initialize >*/
{
	int n;
	itp = interp;
	nf = mf;
	r = rad;
	n = 2*mf+1;
	b1 = sf_floatalloc(n);
	b2 = sf_floatalloc(n);
	d1 = sf_floatalloc(n);
	d2 = sf_floatalloc(n);

	c = lphase(mf, interp);
}

void opwd_fbank(int n1, int n2, float**in, float ****out)
/*< opwd filter bank >*/
{
	int i1, i2, j1, j2;
	float **u1, **u2;

	u1 = sf_floatalloc2(n1, n2);

	for(j2=0; j2<=2*nf; j2++)
	for(j1=0; j1<=2*nf; j1++)
	{
		for(i1=0; i1<n1; i1++)
			firs(-nf, nf, c[j2]+nf, in[0]+i1, n1, n2, u1[0]+i1, n1);
		for(i2=0; i2<n2; i2++)
			firs(-nf, nf, c[j1]+nf, u1[i2], 1, n1, out[j2][j1][i2], 1);
	}

	free(u1[0]);
	free(u1);
}

void opwd(int n1, int n2, float ****fb, float **p, float **out)
/*< apply circle interpolating PWD >*/
{
	int i1, i2, j1, j2;
	float c1, c2;

	b1[0] = 1.0;
	b2[0] = 1.0;
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		c1 = r*sin(p[i2][i1]);
		c2 = r*cos(p[i2][i1]);
		for(j1=1; j1<2*nf+1; j1++)
		{	
			b1[j1] = b1[j1-1]*c1;
			b2[j1] = b2[j1-1]*c2;
		}
		for(j2=0, out[i2][i1] = 0.0; j2<2*nf+1; j2++)
		{
			for(j1=0, c1=0.0; j1<2*nf+1; j1++)
			{
				c2 = ((j1+j2)%2==1?2.0:0.0);
				c1 += c2 * fb[j2][j1][i2][i1] * b1[j1]; 
			}
			out[i2][i1] += c1*b2[j2];
		}
	}
}


void opwdpd(int n1, int n2, float ****fb, 
	float **p, float **out, int id)
/*< partial derivative filter of circle interpolating PWD >*/
{
	int i1, i2, j1, j2;
	float c1, c2;

	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		c1 = r*sin(p[i2][i1]);
		c2 = r*cos(p[i2][i1]);
		if(id == 0)
		{
			b1[0] = 0;
			b1[1] = 1;
			for(j1=2; j1<2*nf+1; j1++)
				b1[j1] = c1*b1[j1-1]*j1/(j1-1);
			b2[0] = 1;
			for(j1=1; j1<2*nf+1; j1++)
				b2[j1] = b2[j1-1]*c2;
		}else{
			b2[0] = 0;
			b2[1] = 1;
			for(j1=2; j1<2*nf+1; j1++)
				b2[j1] = c2*b2[j1-1]*j2/(j1-1);
			b1[0] = 1;
			for(j1=1; j1<2*nf+1; j1++)
				b1[j1] = b1[j1-1]*c1;

		}
		for(j2=0, out[i2][i1] = 0.0; j2<2*nf+1; j2++)
		{
			for(j1=0, c1=0.0; j1<2*nf+1; j1++)
			{
				c2 = ((j1+j2)%2==1?2.0:0.0);
				c1 += c2 * fb[j2][j1][i2][i1] * b1[j1]; 
			}
			out[i2][i1] += c1*b2[j2];
		}
	}
}



void opwd_freq(float dip, int nk, sf_complex**out, bool iir)
/*< frequency response of circle-interpolating PWD >*/
{
	int i1, i2;
	sf_complex c1, c2, z1, z2;

	lphase_filt(nf, c, r*sin(dip), b1, 2*nf, false);
	lphase_filt(nf, c, r*cos(dip), b2, 2*nf, false);

	for(i2=-nk; i2<=nk; i2++)
	{
		c2 = fir_freq(-nf, nf, b2+nf, 0.5*i2/nk);
		for(i1=-nk; i1<=nk; i1++)
		{
			c1 = fir_freq(-nf, nf, b1+nf, 0.5*i1/nk);
			switch(itp)
			{
			case 1:
				z1 = 1.0; z2 = 1.0;
				break;
			case 2:
				z1 = fir_freq(-nf, nf, c[0]+nf, 0.5*i1/nk);
				z2 = fir_freq(-nf, nf, c[0]+nf, 0.5*i2/nk);
				break;
			default:
				z1 = conj(c1);	z2 = conj(c2);
			}
			if(iir)	out[i2+nk][i1+nk] = 1.0 - c1*c2/(z1*z2);
			else	out[i2+nk][i1+nk] = z1*z2 - c1*c2;
		}
	}

}



void opwd_close()
/*< release memory >*/
{
	free(b1);
	free(b2);
	free(d1);
	free(d2);
	free(c[0]);
	free(c);
}

