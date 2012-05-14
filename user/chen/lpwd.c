/* line-interpolating pwd*/

#include <rsf.h>
#include "lphase.h"

//static void *h;
static int nw, n1, n2, itp;
static float *b;
//static float **b1;

void lpwd_init(int interp, int mw, int m1, int m2)
/*< initialize >*/
{

	itp = interp;
	nw = mw;
	n1 = m1;
	n2 = m2;

	b = sf_floatalloc(2*mw+1);
//	b1 = sf_floatalloc2(2*mw+1, m1);
	lphase_init(nw, false, interp);
}

void lpwc(float **in, float **out, float **p)
/*< apply line interpolating PWC >*/
{
	int i1, i2, k;

	switch(itp)
	{
	case 1:		// lagrange fir
		for(i2=0; i2<n2; i2++)
		for(i1=nw; i1<n1-nw; i1++)
		{
			lphase_filt(p[i2][i1], b);
			out[i2][i1] = 0.0;
			for(k=-nw; k<=nw; k++)
				out[i2][i1] += b[k+nw]*in[i2][i1+k];
		}
		break;
	default:	// maxflat iir
		for(i2=0; i2<n2; i2++)
		for(i1=nw; i1<n1-nw; i1++)
		{
			lphase_filt(p[i2][i1], b);
			out[i2][i1] = 0.0;
			for(k=-nw; k<=nw; k++)
				out[i2][i1] += b[k+nw]*in[i2][i1+k];
		}
		break;
	}
}


void lpwd(float **in, float **out, float **p, bool der)
/*< apply line interpolating PWD >*/
{
	int i1, i2, j1;
	
	if(itp==1 || itp == 2)
	{
		for(i2=0; i2<n2-1; i2++)
		for(i1=nw; i1<n1-nw; i1++)
		{
			if(der) lphase_dfilt(p[i2][i1], b);
			else lphase_filt(p[i2][i1], b);
			out[i2][i1] = in[i2][i1];
			for(j1=-nw; j1<=nw; j1++)
				out[i2][i1] -= b[j1+nw]*in[i2+1][i1+j1];
		}
		return;
	}

	for(i2=0; i2<n2-1; i2++)
	for(i1=nw; i1<n1-nw; i1++)
	{
		if(der) lphase_dfilt(p[i2][i1], b);
		else lphase_filt(p[i2][i1], b);
		out[i2][i1] = 0.0;
		for(j1=-nw; j1<=nw; j1++)
			out[i2][i1] += b[j1+nw]*(in[i2][i1-j1]-in[i2+1][i1+j1]);
	}
}

void lpwd_freq(float dip, int nk, sf_complex**out, bool iir)
/*< frequency response of line-interpolating PWD >*/
{
	int i1, i2, j1;
	sf_complex c1, c2;

	lphase_filt(tan(dip), b);
	for(i2=-nk; i2<=nk; i2++)
	for(i1=-nk; i1<=nk; i1++)
	switch(itp)
	{
	case 1:
	case 2:
		for(j1=-nw, c1=0.0; j1<=nw; j1++)
			c1 += b[j1+nw]*cexpf(sf_cmplx(0.,-2*SF_PI*j1*i1/nk));
		c2 = 1.0 - cexpf(sf_cmplx(0.,2*SF_PI*i2/nk))*c1;
		out[i2+nk][i1+nk] = c2;
		break;
	default:
		c1 = sf_cmplx(0, 0);
		c2 = sf_cmplx(0, 0);
		if(iir)
		{
			for(j1=-nw; j1<=nw; j1++)
			{
				c1 += b[j1+nw]*cexpf(sf_cmplx(0.,2*SF_PI*j1*i1/nk));
				c2 += b[j1+nw]*cexpf(sf_cmplx(0.,-2*SF_PI*j1*i1/nk));
			}
			out[i2+nk][i1+nk] = 1-c1/c2*cexpf(sf_cmplx(0.,-2*SF_PI*i2/nk));
		}else{
			for(j1=-nw; j1<=nw; j1++)
				c1 += b[j1+nw]*(cexpf(sf_cmplx(0.,2*SF_PI*j1*i1/nk))-
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
	lphase_close();
	free(b);
}

