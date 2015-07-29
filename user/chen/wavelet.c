

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

#include<rsf.h>

void sf_wvlt_frck(int nf,
		  sf_complex *buf,	/* in frequency; out = ricker wavelet */
		  float *par)		/* 0 ref freq, 1 phase */
/*< frequency domain ricker wavelet >*/
{
	int ifr;
	double freq;
	for(ifr=0;ifr<nf;ifr++)
	{
		freq=crealf(buf[ifr])/par[0];
		freq=freq*freq;
		buf[ifr]=sf_cmplx(4.0*sqrt(SF_PI)*freq*exp(-freq),0.0f);
	}
}

void sf_wvlt_rck(int nt,
		 float *buf,		/* in time index; out = ricker wavelet */
		 float *par)		/* 0 ref freq, 1 phase */
/*< time domain ricker wavelet >*/
{
	int it;
	double d1, w0;

	w0 = 2.0*SF_PI*par[0];
	for(it=0; it<nt; it++)
	{
		d1 = w0 * buf[it];
		d1 = d1*d1;
		buf[it] = (1.0 - 0.5*d1)*exp(-0.25*d1);
	}
}


void sf_wvlt_harmonic(int nt,
		      float *buf,		/* in time index; out = ricker wavelet */
		      float *par)		/* 0 ref freq, 1 phase */
/*< harmonics: sin signal >*/
{
	int it;
	double d1, w0;

	w0 = 2.0*SF_PI*par[0];
	for(it=0; it<nt; it++)
	{
		d1 = w0 * buf[it];
		buf[it] = sin(d1+par[1]);
	}
}


void sf_wvlt_sinc(int nt,
		  float *buf,		/* in time index; out = sinc wavelet */
		  float *par)		/* 0- fw; */
/*< truncated sinc wavelet >*/
{
	int it;
	float x;

	x = SF_PI*par[0];
	for(it=0;it<nt;it++)
	{
		if(buf[it] == 0.0) buf[it] = 1.0;
		else	buf[it] = sin(x*buf[it])/(SF_PI*buf[it]);
	}
}


