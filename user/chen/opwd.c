/* circle-interpolating pwd*/

#include <rsf.h>
#include "lphase.h"

//static void *h;
static int nw, n1, n2, itp;
static float *b1, *b2, *b3, *b4;
static float r;
//static float **b1;

void opwd_init(int interp, int mw, int m1, int m2, float rad)
/*< initialize >*/
{

	itp = interp;
	nw = mw;
	n1 = m1;
	n2 = m2;
	r = rad;
	b1 = sf_floatalloc(2*mw+1);
	b2 = sf_floatalloc(2*mw+1);
	b3 = sf_floatalloc(2*mw+1);
	b4 = sf_floatalloc(2*mw+1);
	lphase_init(nw, false, interp);
}


void opwd(float **in, float **out, float **p, bool der)
/*< apply circle interpolating PWD >*/
{
	int i1, i2, j1, j2;
	float *d1, *d2;

	if(itp == 1 || itp == 2)
	{
		for(i2=nw; i2<n2-nw; i2++)
		for(i1=nw; i1<n1-nw; i1++)
		{
			if(der){
				lphase_dfilt(r*sin(p[i2][i1]), b1);
				lphase_dfilt(r*cos(p[i2][i1]), b2);
			}else{
				lphase_filt(r*sin(p[i2][i1]), b1);
				lphase_filt(r*cos(p[i2][i1]), b2);
			}
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
			lphase_filt(r*sin(p[i2][i1]), b1);
			lphase_filt(r*cos(p[i2][i1]), b2);
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
		lphase_filt(r*sin(p[i2][i1]), b1);
		lphase_filt(r*cos(p[i2][i1]), b2);
		lphase_dfilt(r*sin(p[i2][i1]), d1);
		lphase_dfilt(r*cos(p[i2][i1]), d2);
		out[i2][i1] = 0.0;
		for(j1=-nw; j1<=nw; j1++)
		for(j2=-nw; j2<=nw; j2++)
			out[i2][i1] += 
				(d1[j1+nw]*b2[j2+nw]*cos(p[i2][i1])-
				b1[j1+nw]*d2[j2+nw]*sin(p[i2][i1]))*r*
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

	lphase_filt(r*sin(dip), b1);
	lphase_filt(r*cos(dip), b2);
	lphase_filt(0.0, b3);
	for(i2=-nk; i2<=nk; i2++)
	for(i1=-nk; i1<=nk; i1++)
	switch(itp)
	{
	case 1:
		for(j1=-nw, c1=0.0, c2=0.0; j1<=nw; j1++)
		{
			c1 += b1[j1+nw]*cexpf(sf_cmplx(0.,-2*SF_PI*j1/nk*i1));
			c2 += b2[j1+nw]*cexpf(sf_cmplx(0.,-2*SF_PI*j1/nk*i2));
		}
		out[i2+nk][i1+nk] = 1-c1*c2;
		break;
	case 2:
		for(j1=-nw, c1=0.0, c2=0.0, c3=0.0, c4=0.0; j1<=nw; j1++)
		{
			c1 += b1[j1+nw]*cexpf(sf_cmplx(0., -2*SF_PI*j1/nk*i1));
			c2 += b2[j1+nw]*cexpf(sf_cmplx(0., -2*SF_PI*j1/nk*i2));
			c3 += b3[j1+nw]*cexpf(sf_cmplx(0., -2*SF_PI*j1/nk*i1));
			c4 += b3[j1+nw]*cexpf(sf_cmplx(0., -2*SF_PI*j1/nk*i2));
		}
		if(iir)
			out[i2+nk][i1+nk] = 1- c1*c2/(c3*c4);
		else 
			out[i2+nk][i1+nk] = c3*c4- c1*c2;
		break;
	default:
		for(j1=-nw, c1=0.0, c2=0.0, c3=0.0, c4=0.0; j1<=nw; j1++)
		{
			c1 += b1[j1+nw]*cexpf(sf_cmplx(0.,-2*SF_PI*j1/nk*i1));
			c2 += b2[j1+nw]*cexpf(sf_cmplx(0.,-2*SF_PI*j1/nk*i2));
			c3 += b1[j1+nw]*cexpf(sf_cmplx(0., 2*SF_PI*j1/nk*i1));
			c4 += b2[j1+nw]*cexpf(sf_cmplx(0., 2*SF_PI*j1/nk*i2));
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
	lphase_close();
	free(b1);
	free(b2);
}

