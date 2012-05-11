/* circle-interpolating pwd*/

#include <rsf.h>
#include "lphase.h"

static void *h;
static int nw, n1, n2, itp;
static float *b1, *b2;
//static float **b1;

void opwd_init(int interp, int mw, int m1, int m2)
/*< initialize >*/
{

	itp = interp;
	nw = mw;
	n1 = m1;
	n2 = m2;

	b1 = sf_floatalloc(2*mw+1);
	b2 = sf_floatalloc(2*mw+1);
	h = lphase_init(nw, false, interp);
}


void opwd(float **in, float **out, float **p, bool der)
/*< apply circle interpolating PWD >*/
{
	int i1, i2, j1, j2;
	float *d1, *d2;

	if(itp == 1)
	{
		for(i2=nw; i2<n2-nw; i2++)
		for(i1=nw; i1<n1-nw; i1++)
		{
			lphase_filt(h, sin(p[i2][i1]), b1, der);
			lphase_filt(h, cos(p[i2][i1]), b2, der);
			out[i2][i1] = 0.0;
			for(j1=-nw; j1<=nw; j1++)
			for(j2=-nw; j2<=nw; j2++)
				out[i2][i1] += b1[j1+nw]*b2[j2+nw]*in[i2+j2][i1+j1];
		}
		return;
	}

	// itp=0

	if(!der)
	{	
		for(i2=nw; i2<n2-nw; i2++)
		for(i1=nw; i1<n1-nw; i1++)
		{
			lphase_filt(h, sin(p[i2][i1]), b1, false);
			lphase_filt(h, cos(p[i2][i1]), b2, false);
			out[i2][i1] = 0.0;
			for(j1=-nw; j1<=nw; j1++)
			for(j2=-nw; j2<=nw; j2++)
				out[i2][i1] += b1[j1+nw]*b2[j2+nw]*
					(in[i2-j2][i1-j1] - in[i2+j2][i1+j1]);
		}
		return;
	}

	d1 = sf_floatalloc(2*nw+1);
	d2 = sf_floatalloc(2*nw+1);
	for(i2=nw; i2<n2-nw; i2++)
	for(i1=nw; i1<n1-nw; i1++)
	{
		lphase_filt(h, sin(p[i2][i1]), b1, false);
		lphase_filt(h, cos(p[i2][i1]), b2, false);
		lphase_filt(h, sin(p[i2][i1]), d1, true);
		lphase_filt(h, cos(p[i2][i1]), d2, true);
		out[i2][i1] = 0.0;
		for(j1=-nw; j1<=nw; j1++)
		for(j2=-nw; j2<=nw; j2++)
			out[i2][i1] += 
				(d1[j1+nw]*b2[j2+nw]*cos(p[i2][i1])-
				b1[j1+nw]*d2[j2+nw]*sin(p[i2][i1]))*
				(in[i2-j2][i1-j1] - in[i2+j2][i1+j1]);
	}

	free(d1);
	free(d2);
}

void opwd_freq(float dip, int nk, sf_complex**out, bool iir)
/*< frequency response of circle-interpolating PWD >*/
{
	int i1, i2, j1;
	sf_complex c1, c2, c3, c4;

	lphase_filt(h, sin(dip), b1, false);
	lphase_filt(h, cos(dip), b2, false);
	for(i2=-nk; i2<=nk; i2++)
	for(i1=-nk; i1<=nk; i1++)
	switch(itp)
	{
	case 1:
		for(j1=-nw, c1=0.0, c2=0.0; j1<=nw; j1++)
		{
			c1 += b1[j1+nw]*cexpf(sf_cmplx(0.,-2*SF_PI*j1*i1/nk));
			c2 += b2[j1+nw]*cexpf(sf_cmplx(0.,-2*SF_PI*j1*i2/nk));
		}
		out[i2+nk][i1+nk] = 1-c1*c2;
		break;
	default:
		for(j1=-nw, c1=0.0, c2=0.0, c3=0.0, c4=0.0; j1<=nw; j1++)
		{
			c1 += b1[j1+nw]*cexpf(sf_cmplx(0.,-2*SF_PI/nk*i1*j1));
			c2 += b2[j1+nw]*cexpf(sf_cmplx(0.,-2*SF_PI/nk*i2*j1));
			c3 += b1[j1+nw]*cexpf(sf_cmplx(0., 2*SF_PI/nk*i1*j1));
			c4 += b2[j1+nw]*cexpf(sf_cmplx(0., 2*SF_PI/nk*i2*j1));
		}
		if(iir)
			out[i2+nk][i1+nk] = 1- c1*c2/(c3*c4);
		else 
			out[i2+nk][i1+nk] = c3*c4- c1*c2;
		break;
	}
}



void opwd_close()
/*< release memory >*/
{
	lphase_close(h);
	free(b1);
	free(b2);
}

