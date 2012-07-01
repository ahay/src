/* circle-interpolating pwd*/

#include <rsf.h>
#include "lphase.h"
#include "dsp.h"

static int nf, n1, n2, itp;
static float *b1, *b2, *d1, *d2;
static float r, **c;

void opwd_init(int interp, int mf, int m1, int m2, float rad)
/*< initialize >*/
{

	itp = interp;
	nf = mf;
	n1 = m1;
	n2 = m2;
	r = rad;
	b1 = sf_floatalloc(2*mf+1);
	b2 = sf_floatalloc(2*mf+1);
	d1 = sf_floatalloc(2*mf+1);
	d2 = sf_floatalloc(2*mf+1);

	c = lphase(mf, interp);
}

void opwd(float **in, float **out, float **p)
/*< apply circle interpolating PWD >*/
{
	int i1, i2;
	float c1, c2;

	for(i2=nf; i2<n2-nf; i2++)
	for(i1=nf; i1<n1-nf; i1++)
	{
		lphase_filt(nf, c, r*sin(p[i2][i1]), b1, 2*nf, false);
		lphase_filt(nf, c, r*cos(p[i2][i1]), b2, 2*nf, false);
		c1 = fir2(-nf, nf, b1+nf, -nf, nf, b2+nf, in[i2]+i1, 1, n1);
		switch(itp)
		{
		case 1:
			c2 = in[i2][i1];
			break;
		case 2:
			c2 = fir2(-nf, nf, c[0]+nf, -nf, nf, c[0]+nf, in[i2]+i1, 1, n1);
			break;
		default:
			c2 = fir2(-nf, nf, b1+nf, -nf, nf, b2+nf, in[i2]+i1, -1, -n1);
		}
		out[i2][i1] = c2 - c1;
	}
}

void opwdpd(float **in, float **out, float **p, int id)
/*< partial derivative filter of circle interpolating PWD >*/
{
	int i1, i2;
	float p1, p2, c1, z1;

	for(i2=nf; i2<n2-nf; i2++)
	for(i1=nf; i1<n1-nf; i1++)
	{
		p1 = r*sin(p[i2][i1]);
		p2 = r*cos(p[i2][i1]);
		if(id==0)
		{
			lphase_filt(nf, c, p1, b1, 2*nf, true);
			lphase_filt(nf, c, p2, b2, 2*nf, false);
		} else{
			lphase_filt(nf, c, p1, b1, 2*nf, false);
			lphase_filt(nf, c, p2, b2, 2*nf, true);
		}
		z1 = fir2(-nf, nf, b1+nf, -nf, nf, b2+nf, in[i2]+i1, 1, n1);
		switch(itp)
		{
		case 1:
			c1 = 0.0;
			break;
		case 2:
			c1 = 0.0;
			break;
		default:
			z1 = fir2(-nf, nf, b1+nf, -nf, nf, b2+nf, in[i2]+i1, -1, -n1);
		}
		out[i2][i1] = (c1 - z1);//*(id==0?p2:-p1);
	}
}



void opwdd(float **in, float **out, float **p)
/*< derivative filter of circle interpolating PWD >*/
{
	int i1, i2;
	float p1, p2, c1, c2, z1, z2;

	for(i2=nf; i2<n2-nf; i2++)
	for(i1=nf; i1<n1-nf; i1++)
	{
		p1 = r*sin(p[i2][i1]);
		p2 = r*cos(p[i2][i1]);
		lphase_filt(nf, c, p1, b1, 2*nf, false);
		lphase_filt(nf, c, p2, b2, 2*nf, false);
		lphase_filt(nf, c, p1, d1, 2*nf, true);
		lphase_filt(nf, c, p2, d2, 2*nf, true);
		z1 = -p1*fir2(-nf, nf, b1+nf, -nf, nf, d2+nf, in[i2]+i1, 1, n1);
		z2 =  p2*fir2(-nf, nf, d1+nf, -nf, nf, b2+nf, in[i2]+i1, 1, n1);
		switch(itp)
		{
		case 1:
			c1 = 0.0; c2 = 0.0;
			break;
		case 2:
			c1 = 0.0; c2 = 0.0;
			break;
		default:
			c1 = -p1*fir2(-nf, nf, b1+nf, -nf, nf, d2+nf, in[i2]+i1, -1, -n1);
			c2 =  p2*fir2(-nf, nf, d1+nf, -nf, nf, b2+nf, in[i2]+i1, -1, -n1);
		}
		out[i2][i1] = c1 + c2 - z1 - z2;
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
/*
	switch(itp)
	{
	case 1:
		lphlag_filt(nf, c, r*sin(dip), b1);
		lphlag_filt(nf, c, r*cos(dip), b2);
		for(i2=-nk; i2<=nk; i2++)
		for(i1=-nk; i1<=nk; i1++)
		{
			for(j1=-nf, c1=0.0, c2=0.0; j1<=nf; j1++)
			{
				c1 += b1[j1+nf]*cexpf(sf_cmplx(0., SF_PI*j1/nk*i1));
				c2 += b2[j1+nf]*cexpf(sf_cmplx(0., SF_PI*j1/nk*i2));
			}
			out[i2+nk][i1+nk] = 1-c1*c2;
		}
		break;
	case 2:
		lphbsp_filt(nf, c, r*sin(dip), b1);
		lphbsp_filt(nf, c, r*cos(dip), b2);
		for(j1=0; j1<2*nf+1; j1++)	b3[j1] = c[j1][0];
		for(i2=-nk; i2<=nk; i2++)
		for(i1=-nk; i1<=nk; i1++)
		{
			for(j1=-nf, c1=0.0, c2=0.0, c3=0.0, c4=0.0; j1<=nf; j1++)
			{
				c1 += b1[j1+nf]*cexpf(sf_cmplx(0., -SF_PI*j1/nk*i1));
				c2 += b2[j1+nf]*cexpf(sf_cmplx(0., -SF_PI*j1/nk*i2));
				c3 += b3[j1+nf]*cexpf(sf_cmplx(0., -SF_PI*j1/nk*i1));
				c4 += b3[j1+nf]*cexpf(sf_cmplx(0., -SF_PI*j1/nk*i2));
			}
			if(iir)
				out[i2+nk][i1+nk] = 1- c1*c2/(c3*c4);
			else 
				out[i2+nk][i1+nk] = c3*c4- c1*c2;
		}
		break;
	default:
		lphmf_filt(nf, c, r*sin(dip), b1);
		lphmf_filt(nf, c, r*cos(dip), b2);
		for(i2=-nk; i2<=nk; i2++)
		for(i1=-nk; i1<=nk; i1++)
		{
			for(j1=-nf, c1=0.0, c2=0.0, c3=0.0, c4=0.0; j1<=nf; j1++)
			{
				c1 += b1[j1+nf]*cexpf(sf_cmplx(0.,-SF_PI*j1/nk*i1));
				c2 += b2[j1+nf]*cexpf(sf_cmplx(0.,-SF_PI*j1/nk*i2));
				c3 += b1[j1+nf]*cexpf(sf_cmplx(0., SF_PI*j1/nk*i1));
				c4 += b2[j1+nf]*cexpf(sf_cmplx(0., SF_PI*j1/nk*i2));
			}
			if(iir)
				out[i2+nk][i1+nk] = 1- c1*c2/(c3*c4);
			else 
				out[i2+nk][i1+nk] = c3*c4- c1*c2;
		}
	}
*/
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

