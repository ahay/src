/* line-interpolating pwd*/

#include <rsf.h>
#include "lphmf.h"
#include "lphlag.h"
#include "lphbsp.h"

static int nf, itp;
static int n1, n2;
static float **c, *b1, *b2;

void lpwd_init(int interp, int mf, int m1, int m2)
/*< initialize >*/
{

	itp = interp;
	nf = mf;
	n1 = m1;
	n2 = m2;

	switch(interp)
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
	b1 = sf_floatalloc(2*mf+1);
	b2 = sf_floatalloc(2*mf+1);
}

void lpwc(float **in, float **out, float **p)
/*< apply line interpolating PWC >*/
{
	int i1, i2, j1;
	void (*filt)(int, float**, float, float*);

	switch(itp)
	{
	case 1:
		filt=lphlag_filt;
		break;
	case 2:
		filt=lphbsp_filt;
		break;
	default:
		filt=lphmf_filt;
	}

	for(i2=0; i2<n2; i2++)
	for(i1=nf; i1<n1-nf; i1++)
	{
		filt(nf, c, p[i2][i1], b1);
		out[i2][i1] = 0.0;
		for(j1=-nf; j1<=nf; j1++)
			out[i2][i1] += b1[j1+nf]*in[i2][i1-j1];
	}
}


void lpwd(float **in, float **out, float **p, bool der)
/*< apply line interpolating PWD >*/
{
	int i1, i2, j1;
	
	switch(itp)
	{
	case 1:
		for(i2=0; i2<n2-1; i2++)
		for(i1=nf; i1<n1-nf; i1++)
		{
			if(der) lphlag_dfilt(nf, c, p[i2][i1], b1);
			else lphlag_filt(nf, c, p[i2][i1], b1);
			out[i2][i1] = in[i2][i1];
			for(j1=-nf; j1<=nf; j1++)
				out[i2][i1] -= b1[j1+nf]*in[i2+1][i1-j1];
		}
		break;
	case 2:
		for(i2=0; i2<n2-1; i2++)
		for(i1=nf; i1<n1-nf; i1++)
		{
			if(der) lphbsp_dfilt(nf, c, p[i2][i1], b1);
			else lphbsp_filt(nf, c, p[i2][i1], b1);
			out[i2][i1] = in[i2][i1];
			for(j1=-nf; j1<=nf; j1++)
				out[i2][i1] -= b1[j1+nf]*in[i2+1][i1-j1];
		}
		break;
	default:
		for(i2=0; i2<n2-1; i2++)
		for(i1=nf; i1<n1-nf; i1++)
		{
			if(der) lphmf_dfilt(nf, c, p[i2][i1], b1);
			else lphmf_filt(nf, c, p[i2][i1], b1);
			out[i2][i1] = 0.0;
			for(j1=-nf; j1<=nf; j1++)
				out[i2][i1] += b1[j1+nf]*(in[i2][i1+j1]-in[i2+1][i1-j1]);
		}
	}
}

void lpwd_freq(float dip, int nk, sf_complex**out, bool iir)
/*< frequency response of line-interpolating PWD >*/
{
	int i1, i2, j1;
	sf_complex c1, c2, z2;

	switch(itp)
	{
	case 1:
		lphlag_filt(nf, c, tan(dip), b1);
		break;
	case 2:
		lphbsp_filt(nf, c, tan(dip), b1);
		break;
	default:
		lphmf_filt(nf, c, tan(dip), b1);
	}
	for(i2=-nk; i2<=nk; i2++)
	for(i1=-nk; i1<=nk; i1++)
	switch(itp)
	{
	case 1:
	case 2:
		for(j1=-nf, c1=0.0; j1<=nf; j1++)
			c1 += b1[j1+nf]*cexpf(sf_cmplx(0.,-2*SF_PI*j1*i1/nk));
		c2 = 1.0 - cexpf(sf_cmplx(0.,2*SF_PI*i2/nk))*c1;
		out[i2+nk][i1+nk] = c2;
		break;
	default:
		c1 = sf_cmplx(0, 0);
		c2 = sf_cmplx(0, 0);
		if(iir)
		{
			for(j1=-nf; j1<=nf; j1++)
			{
				c1 += b1[j1+nf]*cexpf(sf_cmplx(0.,2*SF_PI*j1*i1/nk));
				c2 += b1[j1+nf]*cexpf(sf_cmplx(0.,-2*SF_PI*j1*i1/nk));
			}
			out[i2+nk][i1+nk] = 1-c1/c2*cexpf(sf_cmplx(0.,-2*SF_PI*i2/nk));
		}else{
			for(j1=-nf; j1<=nf; j1++)
				c1 += b1[j1+nf]*(cexpf(sf_cmplx(0.,2*SF_PI*j1*i1/nk))-
					cexpf(sf_cmplx(0., -2*SF_PI*j1*i1/nk))*
					cexpf(sf_cmplx(0., 2*SF_PI*i2/nk)));
			out[i2+nk][i1+nk] = c1;
		}
		break;
	}
}


void lpwd_close()
/*< release memory >*/
{
	free(b1);
	free(b2);
	free(c[0]);
	free(c);
}

