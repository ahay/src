/* digital filters */

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

float fir(int ma, int na, float *a,
	  float *x, int d1)
/*< FIR filter
  na
  y(n) = sum a[k]*x[n-k] 
  k=ma
  >*/
{
    int k;
    double out=0.0;
    for(k=ma; k<=na; k++)
	out += a[k]*x[-k*d1];
    return out;
}


float fir2(int m1, int n1, float *a1, 
	   int m2, int n2, float *a2,
	   float *x, int d1, int d2)
/*< FIR filter
  n2    n1
  y[m][n] = sum   sum a2[k2]*a1[k1]*x[m-k2][n-k1] 
  k2=m2 k1=m1
  >*/
{
    int k1, k2;
    double out=0.0;
    for(k2=m2; k2<=n2; k2++)
	for(k1=m1; k1<=n1; k1++)
	    out += a2[k2]*a1[k1]*x[-k1*d1-k2*d2];
    return out;
}

float iir(int na, float *a, int nb, float *b,
	  float *x, int d1, float *y, int d2)
/*< IIR causal filter
  na                nb
  y(n) = ( sum a[k]*x[n-k] - sum b[k]*y[n-k] )/b[0]
  k=ma              k=mb
  >*/
{
    int k;
    double out=0.0;
    for(k=0; k<=na; k++)
	out += a[k]*x[-k*d1];
    for(k=1; k<=nb; k++)
	out -= b[k]*y[-k*d2];
    return out/b[0];
}


void firs(int ma, int na, float *a,
	  float *x, int d1, int n1, 
	  float *y, int d2)
/*< FIR filter
  na
  y(n) = sum a[k]*x[n-k] 
  k=ma
  >*/
{
    int k, n, min, max, sign;
    double out;
    for(n=0; n<n1; n++)
    {
	if(d1>0){
	    min = SF_MAX(ma, n-n1+1);
	    max = SF_MIN(na, n);
	    sign = 1;
	}else{
	    min = SF_MAX(ma, -n);
	    max = SF_MIN(n1-n-1, na);
	    sign = -1;
	}
	out = 0.0;
	for(k=min; k<=max; k++)
	    out += a[k]*x[(n-k*sign)*abs(d1)];
	y[n*d2] = out;
    }
}


void iirs(int na, float *a, int nb, float *b,
	  float *x, int d1, int n1, 
	  float *y, int d2)
/*< IIR filter
  na                nb
  y(n) = ( sum a[k]*x[n-k] - sum b[k]*y[n-k] )/b[0]
  k=0               k=1
  >*/
{
    int k, n;
    double out;
    for(n=0, out=0.0; n<n1; n++)
    {
	for(k=0; (k<=na && k<=n); k++)
	    out += a[k]*x[(n-k)*d1];
	for(k=1; (k<=nb && k<=n); k++)
	    out -= b[k]*y[(n-k)*d2];
	y[n*d2] = out/b[0];
    }
}


void allpass(int na, float *a,
	     float *x, int d1, int n1, 
	     float *y, int d2)
/*< allpass filter
  na                na
  y(n) = ( sum a[k]*x[n-k] - sum a[na-k]*y[n-k] )/a[na]
  k=0               k=1
  >*/
{
    int k, n;
    double out;
    for(n=0, out=0.0; n<n1; n++)
    {
	for(k=0; (k<=na && k<=n); k++)
	    out += a[k]*x[(n-k)*d1];
	for(k=1; (k<=na && k<=n); k++)
	    out -= a[na-k]*y[(n-k)*d2];
	y[n*d2] = out/a[na];
    }
}

sf_complex fir_freq(int m, int n, float *h, float f)
/*< FIR frequency response 
  n
  H(e^{jw}) = sum h[k] e^{-jwk}
  k=m
  >*/
{
    int k;
    sf_complex out=sf_cmplx(0.0,0.0);
    for(k=m; k<=n; k++) {
#ifdef SF_HAS_COMPLEX_H
	out += h[k]*cexpf(sf_cmplx(0, -2*SF_PI*f*k));
#else
	out = sf_cadd(out,sf_cmul(sf_cmplx(h[k],0.0),cexpf(sf_cmplx(0, -2*SF_PI*f*k))));
#endif
    }
    return out;
}


sf_complex iir_freq(int ma, int na, float *a, 
		    int mb, int nb, float *b, float f)
/*< IIR frequency response 
  na
  sum a[k] e^{-jwk}
  k=ma
  H(e^{jw}) = ---------------------
  nb
  sum b[k] e^{-jwk}
  k=mb
  >*/
{
    int k;
    sf_complex z1, z2;
    for(k=ma, z1=sf_cmplx(0.0,0.0); k<=na; k++) {
#ifdef SF_HAS_COMPLEX_H
	z1 += a[k]*cexpf(sf_cmplx(0, -2*SF_PI*f*k));
#else
	z1 = sf_cadd(z1,sf_cmul(sf_cmplx(a[k],0.0),cexpf(sf_cmplx(0, -2*SF_PI*f*k))));
#endif
    }
    for(k=mb, z2=sf_cmplx(0.0,0.0); k<=nb; k++) {
#ifdef SF_HAS_COMPLEX_H
	z2 += b[k]*cexpf(sf_cmplx(0, -2*SF_PI*f*k));
#else
	z2 = sf_cadd(z2,sf_cmul(sf_cmplx(b[k],0.0),cexpf(sf_cmplx(0, -2*SF_PI*f*k))));
#endif
    }
#ifdef SF_HAS_COMPLEX_H
    return (z1*conjf(z2))/(z2*conjf(z2)+0.00001);
#else
    return sf_cdiv(sf_cmul(z1,conjf(z2)),sf_cadd(sf_cmul(z2,conjf(z2)),sf_cmplx(0.00001,0.0)));
#endif
}



sf_complex allpass_freq(int n, float *a, float f)
/*< frequency response of allpass filter
  n
  sum a[k] e^{-jwk}
  k=-n
  H(e^{jw}) = ---------------------
  n
  sum a[k] e^{jwk}
  k=-n
  >*/
{
    int k;
    sf_complex z0;
    for(k=-n, z0=sf_cmplx(0.0,0.0); k<=n; k++) {
#ifdef SF_HAS_COMPLEX_H
	z0 += a[k]*cexpf(sf_cmplx(0, -2*SF_PI*f*k));
#else
	z0 = sf_cadd(z0,sf_cmul(sf_cmplx(a[k],0.0),cexpf(sf_cmplx(0, -2*SF_PI*f*k))));
#endif
    }
#ifdef SF_HAS_COMPLEX_H
    return (z0*z0)/(z0*conjf(z0)+0.00001);
#else
    return sf_cdiv(sf_cmul(z0,z0),sf_cadd(sf_cmul(z0,conjf(z0)),sf_cmplx(0.00001,0.0)));
#endif
}


