/* circle-interpolating pwd*/

#include <rsf.h>
#include "lphlag.h"
#include "lphbsp.h"
#include "lphmf.h"

static int nf, n1, n2, itp;
static float *b1, *b2, *b3;
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
	b3 = sf_floatalloc(2*mf+1);
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
}


void opwd(float **in, float **out, float **p, bool der)
/*< apply circle interpolating PWD >*/
{
	int i1, i2, j1, j2;
	float *d1, *d2;

	switch(itp)
	{
	case 1:
		for(i2=nf; i2<n2-nf; i2++)
		for(i1=nf; i1<n1-nf; i1++)
		{
			if(der){
				lphlag_dfilt(nf, c, r*sin(p[i2][i1]), b1);
				lphlag_dfilt(nf, c, r*cos(p[i2][i1]), b2);
			}else{
				lphlag_filt(nf, c, r*sin(p[i2][i1]), b1);
				lphlag_filt(nf, c, r*cos(p[i2][i1]), b2);
			}
			out[i2][i1] = 0.0;
			for(j1=-nf; j1<=nf; j1++)
			for(j2=-nf; j2<=nf; j2++)
				out[i2][i1] += b1[j1+nf]*b2[j2+nf]*in[i2+j2][i1+j1];
		}
		break;
	case 2:
		for(i2=nf; i2<n2-nf; i2++)
		for(i1=nf; i1<n1-nf; i1++)
		{
			if(der){
				lphbsp_dfilt(nf, c, r*sin(p[i2][i1]), b1);
				lphbsp_dfilt(nf, c, r*cos(p[i2][i1]), b2);
			}else{
				lphbsp_filt(nf, c, r*sin(p[i2][i1]), b1);
				lphbsp_filt(nf, c, r*cos(p[i2][i1]), b2);
			}
			out[i2][i1] = 0.0;
			for(j1=-nf; j1<=nf; j1++)
			for(j2=-nf; j2<=nf; j2++)
				out[i2][i1] += b1[j1+nf]*b2[j2+nf]*in[i2+j2][i1+j1];
		}
		break;
	default:
		if(!der)
		{	
			for(i2=nf; i2<n2-nf; i2++)
			for(i1=nf; i1<n1-nf; i1++)
			{
				lphmf_filt(nf, c, r*sin(p[i2][i1]), b1);
				lphmf_filt(nf, c, r*cos(p[i2][i1]), b2);
				out[i2][i1] = 0.0;
				for(j1=-nf; j1<=nf; j1++)
				for(j2=-nf; j2<=nf; j2++)
					out[i2][i1] += b1[j1+nf]*b2[j2+nf]*
						(in[i2-j2][i1-j1] - in[i2+j2][i1+j1]);
			}
			break;
		}

		d1 = sf_floatalloc(2*nf+1);
		d2 = sf_floatalloc(2*nf+1);
		for(i2=nf; i2<n2-nf; i2++)
		for(i1=nf; i1<n1-nf; i1++)
		{
			lphmf_filt(nf, c, r*sin(p[i2][i1]), b1);
			lphmf_filt(nf, c, r*cos(p[i2][i1]), b2);
			lphmf_dfilt(nf, c, r*sin(p[i2][i1]), d1);
			lphmf_dfilt(nf, c, r*cos(p[i2][i1]), d2);
			out[i2][i1] = 0.0;
			for(j1=-nf; j1<=nf; j1++)
			for(j2=-nf; j2<=nf; j2++)
				out[i2][i1] += 
					(d1[j1+nf]*b2[j2+nf]*cos(p[i2][i1])-
					b1[j1+nf]*d2[j2+nf]*sin(p[i2][i1]))*r*
					(in[i2-j2][i1-j1] - in[i2+j2][i1+j1]);
		}

		free(d1);
		free(d2);
	}
}

void opwd_freq(float dip, int nk, sf_complex**out, bool iir)
/*< frequency response of circle-interpolating PWD >*/
{
	int i1, i2, j1;
	sf_complex c1, c2, c3, c4;

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
}



void opwd_close()
/*< release memory >*/
{
	free(b1);
	free(b2);
	free(b3);
	free(c[0]);
	free(c);
}

