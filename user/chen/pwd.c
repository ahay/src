/* plane-wave destruction by line/circle interpolating */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <rsf.h>
#include "pcmf.h"

static void *h;
static int nw, n1, n2;
static double **b;

void pwd_init(int mw, int m1, int m2)
/*< initialize >*/
{
	int i1;
	h = pcmf_init(mw);
	nw = mw;
	n1 = m1;
	n2 = m2;
	b = sf_alloc((2*mw+1), sizeof(double*));
	b[0] = sf_alloc((2*mw+1)*(2*nw+1), sizeof(double));
	for(i1=1; i1<2*nw+1; i1++)
		b[i1] =b[i1-1] + 2*nw+1;
}

void lpwd(float **in, float **out, float **p)
/*< apply line interpolating PWD >*/
{
	int i1, i2, k;

	memset(out[0], 0, n1*n2*sizeof(float));
	for(i2=0; i2<n2-1; i2++)
	{
		for(i1=nw; i1<n1-nw; i1++)
		{
			pcmf_filt_1d(h, p[i2][i1], b[0]);
			for(k=-nw; k<=nw; k++)
				out[i2][i1] += b[0][nw-k]*(in[i2][i1+k]-in[i2+1][i1-k]);
		}
	}
}

void lpwd_der(float **in, float **out, float **p)
/*< apply derivative filter of line interpolating PWD >*/
{
	int i1, i2, k;

	memset(out[0], 0, n1*n2*sizeof(float));
	for(i2=0; i2<n2-1; i2++)
	{
		for(i1=nw; i1<n1-nw; i1++)
		{
			pcmf_der_1d(h, p[i2][i1], b[0]);
			for(k=-nw; k<=nw; k++)
				out[i2][i1] += b[0][nw-k]*(in[i2][i1+k]-in[i2+1][i1-k]);
		}
	}
}


void opwd(float **in, float **out, float **p)
/*< apply circle interpolating PWD >*/
{
	int i1, i2, j1, j2;

	memset(out[0], 0, n1*n2*sizeof(float));
	for(i2=nw; i2<n2-nw; i2++)
	{
		for(i1=nw; i1<n1-nw; i1++)
		{
			pcmf_filt_2d(h, p[i2][i1], b);
			for(j1=-nw; j1<=nw; j1++)
			for(j2=-nw; j2<=nw; j2++)
				out[i2][i1] += b[nw-j1][nw-j2]*
					(in[i2+j2][i1+j1]-in[i2-j2][i1-j1]);
		}
	}
}

void opwd_der(float **in, float **out, float **p)
/*< apply the derivative operator of circle interpolating PWD >*/
{
	int i1, i2, j1, j2;

	memset(out[0], 0, n1*n2*sizeof(float));
	for(i2=nw; i2<n2-nw; i2++)
	{
		for(i1=nw; i1<n1-nw; i1++)
		{
			pcmf_der_2d(h, p[i2][i1], b);
			for(j1=-nw; j1<=nw; j1++)
			for(j2=-nw; j2<=nw; j2++)
				out[i2][i1] += b[nw-j1][nw-j2]*
					(in[i2+j2][i1+j1]-in[i2-j2][i1-j1]);
		}
	}
}

void lpwd_freq(float p, int nk, sf_complex**out, bool iir)
/*< frequency response of line-interpolating PWD >*/
{
	int i1, i2, j1;
	sf_complex c1, c2;

	pcmf_filt_1d(h, p, b[0]);
	for(i2=-nk; i2<=nk; i2++)
	for(i1=-nk; i1<=nk; i1++)
	if(iir)
	{
		for(j1=-nw, c1=0.0, c2=0.0; j1<=nw; j1++)
		{
		    c1 += b[0][j1+nw]*cexpf(sf_cmplx(0.,2*SF_PI*j1*i1/nk));
		    c2 += b[0][j1+nw]*cexpf(sf_cmplx(0.,-2*SF_PI*j1*i1/nk));
		}
		out[i2+nk][i1+nk] = 1-c1/c2*cexpf(sf_cmplx(0.,2*SF_PI*i2/nk));
	}
	else
	{
		for(j1=-nw, c1=0.0; j1<=nw; j1++)
		    c1 += b[0][j1+nw]*(cexpf(sf_cmplx(0.,-2*SF_PI*j1*i1/nk))-
				       cexpf(sf_cmplx(0.,2*SF_PI*j1*i1/nk))*
				       cexpf(sf_cmplx(0.,2*SF_PI*i2/nk)));
		out[i2+nk][i1+nk] = c1;
	}
}


void opwd_freq(float dip, int nk, sf_complex**out, bool iir)
/*< frequency response of circle-interpolating PWD >*/
{
	int i1, i2, j1, j2;
	sf_complex c1,c2;

	pcmf_filt_2d(h, dip, b);
	for(i2=-nk; i2<=nk; i2++)
	for(i1=-nk; i1<=nk; i1++)
	if(iir)
	{
		c1=0.0;
		c2=0.0;
		for(j1=-nw; j1<=nw; j1++)
		for(j2=-nw; j2<=nw; j2++)
		{
		    c1 += b[j1+nw][j2+nw]*cexpf(sf_cmplx(0.,-2*SF_PI/nk*(i2*j2+i1*j1)));
		    c2 += b[j1+nw][j2+nw]*cexpf(sf_cmplx(0.,+2*SF_PI/nk*(i2*j2+i1*j1)));
		}
		out[i2+nk][i1+nk] = 1- c1/c2;
	}else{
		c1=0.0;
		for(j1=-nw; j1<=nw; j1++)
		for(j2=-nw; j2<=nw; j2++)
			c1 += b[j1+nw][j2+nw]*
			    (cexpf(sf_cmplx(0.,-2*SF_PI/nk*(i2*j2+i1*j1)))-
			     cexpf(sf_cmplx(0.,+2*SF_PI/nk*(i2*j2+i1*j1))));
		out[i2+nk][i1+nk] = c1;
	}
}


void lpwd_phase(float p, int n, float*out)
/*< frequency response of line-interpolating PWD >*/
{
	int i1, j1;
	sf_complex c1, c2;

	pcmf_filt_1d(h, p, b[0]);
	for(i1=0; i1<n; i1++)
	{
		for(j1=-nw, c1=0.0, c2=0.0; j1<=nw; j1++)
		{
		    c1 += b[0][j1+nw]*cexpf(sf_cmplx(0.,+2*SF_PI*j1*i1/n));
		    c2 += b[0][j1+nw]*cexpf(sf_cmplx(0.,-2*SF_PI*j1*i1/n));
		}
		out[i1] = cargf(c1/c2);
	}
}


void pwd_close()
/*< release memory >*/
{
	free(b[0]);
	free(b);
	pcmf_close(h);
}


