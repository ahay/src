/* large dip calculation via non-stationary regularization  */
/*
  Copyright (C) 2012 The University of Texas at Austin
  
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
/*#include "dipl.h"*/

/*from poly.c*/
/* polynomial operators */
/* p(x) = c[0] + c[1]*x + c[2]*x^2 + c[n]*x^n */

static double prod_n_m(double *r, int n, int m)
{
    if( n<m || n<=0) return 0.0;
    if( m == 0 ) return 1.0;
    if( m == 1 ) 
	return (r[0]+prod_n_m(r+1, n-1, 1));
    return (r[0]*prod_n_m(r+1, n-1, m-1) +
	    prod_n_m(r+1, n-1, m));
}

void poly(int n, float *r)
/*< polynomial root to coefficients 
  [1.0, 2.0, 1.0] => [2.0, 3.0, 1.0]
  >*/
{
    double *p;
    int i1;

    p = sf_alloc(n, sizeof(double));
    for(i1=0; i1<n; i1++)	p[i1] = -r[i1];

    for(i1=0; i1<=n; i1++)
    {
	r[i1]=r[n]*prod_n_m(p, n, n-i1);
    }

    free(p);
}

sf_complex poly_dval(int n, float *c, sf_complex x)
/*< polynomial derivative value >*/
{
    int i;
    sf_complex val;

    val = sf_cmplx(0.0,0.0);
    for(i=n; i>=1; i--)
    {
#ifdef SF_HAS_COMPLEX_H
	val = val*x + c[i]*i;
#else
	val = sf_cadd(sf_cmul(val,x),sf_cmplx(c[i]*i,0.0));
#endif
    }
    return val;
}



sf_complex poly_val(int n, float *c, sf_complex x )
/*< polynomial value >*/
{
    int i;
    sf_complex val;

    val = sf_cmplx(0.0,0.0);
    for(i=n; i>=0; i--)
    {
#ifdef SF_HAS_COMPLEX_H
	val = val*x + c[i];
#else
	val = sf_cadd(sf_cmul(val,x),sf_cmplx(c[i],0.0));
#endif
    }
    return val;
}

sf_complex plyr_val(int n, float *r, sf_complex x)
/*< polynomial value >*/
{
    int i;
    sf_complex val;

    val = sf_cmplx(r[n],0.0);
    for(i=0; i<n; i++)
    {
#ifdef SF_HAS_COMPLEX_H
	val = val*(x-r[i]);
#else
	val = sf_cmul(val,sf_cadd(x,sf_cmplx(-r[i],0.0)));
#endif
    }
    return val;
}


sf_complex plyr_dval(int n, float *r, sf_complex x)
/*< polynomial derivative value >*/
{
    int i, k;
    sf_complex v1, v2;

    v1 = sf_cmplx(1.0,0.0);
    k = 0;
    for(i=0; i<n; i++)
    {
#ifdef SF_HAS_COMPLEX_H
	v2 = x - r[i];
#else
	v2 = sf_cadd(x,sf_cmplx(-r[i],0.0));
#endif
	if(crealf(v2) == 0.0f && cimagf(v2) == 0.0f) k++;
#ifdef SF_HAS_COMPLEX_H
	else	v1 *= v2;
#else
	else    v1 = sf_cmul(v1,v2);
#endif
    }
    if(k >= 2) return sf_cmplx(0.0,0.0);
#ifdef SF_HAS_COMPLEX_H
    if(k == 1) return v1*r[n];
#else
    if(k == 1) return sf_crmul(v1,r[n]);
#endif

    v2 = sf_cmplx(0.0,0.0);
    for(i=0; i<n; i++)
    {
#ifdef SF_HAS_COMPLEX_H
	v2 += v1/(x-r[i]);
#else
	v2 = sf_cadd(v2,sf_cdiv(v1,sf_cadd(x,sf_cmplx(-r[i],0.0))));
#endif
    }

#ifdef SF_HAS_COMPLEX_H
    return v2*r[n];
#else
    return sf_crmul(v2,r[n]);
#endif
}


sf_complex cplyr_val(int n, sf_complex *r, sf_complex x)
/*< polynomial value >*/
{
    int i;
    sf_complex val;

    val = r[n];
    for(i=0; i<n; i++)
    {
#ifdef SF_HAS_COMPLEX_H
	val = val*(x-r[i]);
#else
	val = sf_cmul(val,sf_cadd(x,sf_crmul(r[i],-1.0)));
#endif
    }
    return val;
}

sf_complex cplyr_dval(int n, sf_complex *r, sf_complex x)
/*< polynomial derivative value >*/
{
    int i, k;
    sf_complex v1, v2;

    v1 = sf_cmplx(1.0,0.0);
    k = 0;
    for(i=0; i<n; i++)
    {
#ifdef SF_HAS_COMPLEX_H
	v2 = x - r[i];
#else
	v2 = sf_cadd(x,sf_crmul(r[i],-1.0));
#endif
	if(crealf(v2) == 0.0f && cimagf(v2) == 0.0f) k++;
#ifdef SF_HAS_COMPLEX_H
	else	v1 *= v2;
#else
	else	v1 = sf_cmul(v1,v2);
#endif
    }
    if(k >= 2) return sf_cmplx(0.0,0.0);
#ifdef SF_HAS_COMPLEX_H
    if(k == 1) return v1*r[n];
#else
    if(k == 1) return sf_cmul(v1,r[n]);
#endif

    v2 = sf_cmplx(0.0,0.0);
    for(i=0; i<n; i++)
    {
#ifdef SF_HAS_COMPLEX_H
	v2 += v1/(x-r[i]);
#else
	v2 = sf_cadd(v2,sf_cdiv(v1,sf_cadd(x,sf_crmul(r[i],-1.0))));
#endif
    }

#ifdef SF_HAS_COMPLEX_H
    return v2*r[n];
#else
    return sf_cmul(v2,r[n]);
#endif
}

/*from lphpoly.c*/
float **lphlag(int m, int n)
/*< Linear phase filter by noncausal Lagrange approximation
[    1     p    p^2  |
   ------------------+---
    0.0   0.5   0.5  |  Z^{m}
    1.0   0.0  -1.0	 |
    0.0  -0.5   0.5  |  Z^{-n} ]
 >*/
{
	float **p, c;
	int i, k, n1, i1, i2;

	n1 = m+n;
	p = sf_floatalloc2(n1+1, n1+1);

	for(k=-m; k<=n; k++)
	{
		p[k+m][n1] = 1.0;
		for(i=-m; i<k; i++)
		{
			p[k+m][n1] /= (i-k);
			p[k+m][i+m] = -i;
		}
		for(i=k+1; i<=n; i++)
		{
			p[k+m][n1] /= (i-k);
			p[k+m][i+m-1] = -i;
		}
		poly(n1, p[k+m]);
	}

	for(i2=0; i2<=n1; i2++)
	for(i1=0; i1<i2; i1++)
	{
		c = p[i2][i1];
		p[i2][i1] = p[i1][i2];
		p[i1][i2] = c;
	}

	return p;
}

static double factorial(int s, int e)
{
	double c=1.0;
	int i;
	for(i=s; i<=e; i++) c*=i;
	return c;	
}

float **lphbsp(int m, int n)
/*< Linear phase filter by noncausal B-Spline approximation
[    1     p    p^2  |
   ------------------+---
   0.125   1.0   1.0  |  Z
   0.75    0.0  -2.0	 |
   0.125  -1.0   1.0  |  Z^{-1} ]
 >*/
{
	float **p;
	double c, c0, c1, c2, c3;
	int i, j, k, n1;

	n1 = m+n;
	p = sf_floatalloc2(n1+1, n1+1);

	c0 = factorial(1, n1+1);
	for(j=0; j<=n1; j++)
	{
		c1 = factorial(1, j) * factorial(1, n1-j);
		for(k=-m; k<=n; k++)
		{
			c = 0;
			for(i=0; i<k+1+0.5*n1; i++)
			{
				c2 = (i%2 == 1)? -1.0 : 1.0;
				c2 *= factorial(1, i) * factorial(1, n1+1-i);
				c3 = pow((n1+1.0)/2+k-i, n1-j);
				c += c0*c3/(c1*c2);
			}
			p[j][k+m] = c;	
		}
	}
	return p;
}


float **lphmf(int m, int n)
/*< Linear phase filter by noncausal maxflat approximation 
[    1     p    p^2  |
   ------------------+---
    0.1667   0.5   0.3333  |  Z
    0.6667   0.0  -0.6667  |
    0.1667  -0.5   0.3333  |  Z^{-1} ]
 >*/
{
	float **p;
	int i, k, n1, i1, i2;
	double c1, c2;

	n1 = m+n;
	p = sf_floatalloc2(n1+1, n1+1);

	c1 = factorial(1, n1);
	c1 = c1/factorial(n1+1, 2*n1);
	for(k=-m; k<=n; k++)
	{
		c2 = ((m+k)%2 == 1)? -1.0 : 1.0;
		c2 *= factorial(1, m+k) * factorial(1, n-k);
		p[k+m][n1] =c1/c2;
		for(i=0; i<m+k; i++)
			p[k+m][i] = (2*m-i);
		for(i=0; i<n-k; i++)
			p[k+m][i+m+k] = (i-2*n);
		poly(n1, p[k+m]);
	}
	for(i2=0; i2<=n1; i2++)
	for(i1=0; i1<i2; i1++)
	{
		c1 = p[i2][i1];
		p[i2][i1] = p[i1][i2];
		p[i1][i2] = c1;
	}

	for(i1=n1; i1>=1; i1--)
	{
		for(i=n1; i>=i1; i--)	
		for(k=-m; k<=n; k++)
		p[i][k+m] *= 2.0;
	}

	return p;
}

float **lphpoly(int m, int n, char* interp)
/*< linear phase filter bank >*/
{
	if(strcmp(interp,"lagrange")==0)	return lphlag(m, n);
	else if(strcmp(interp,"bspline")==0)	return lphbsp(m, n);
	else if(strcmp(interp,"maxflat")==0)	return lphmf(m, n);
	else return NULL;
}


void lphpoly_coef(int nf, float **c, float delay, 
	float *out, int np, bool der)
/*< filter >*/
{
	int k, i;
	double b;

	for(k=0; k<=nf; k++)
	{
		b = der?0:c[0][k];
		for(i=1; i<=np; i++)
			b += c[i][k] * 
				(der?(powf(delay, i-1)*i):powf(delay, i));
		out[k] = b;
	}
}


/*from dsp.c*/
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



/*from opwd.c*/
static int nf1, nf2;
static float *b1, *b2;
static float r, **c;

void opwd_init(int mf1, int mf2, char *interp, float rad)
/*< initialize >*/
{

	nf1 = mf1;
	nf2 = mf2;
	r = rad;
	b1 = sf_floatalloc(mf1+mf2+1);
	b2 = sf_floatalloc(mf1+mf2+1);

	c = lphpoly(mf1, mf2, interp);
}

void opwd(int n1, int n2, float **in, sf_complex **p, float **out)
/*< apply circle interpolating PWD >*/
{
	int i1, i2;
	float c1, c2;

	for(i2=nf1; i2<n2-nf2; i2++)
	for(i1=nf1; i1<n1-nf2; i1++)
	{
		lphpoly_coef(nf1+nf2, c, r*cimagf(p[i2][i1]), b1, nf1+nf2, false);
		lphpoly_coef(nf1+nf2, c, r*crealf(p[i2][i1]), b2, nf1+nf2, false);
		c1 = fir2(-nf1, nf2, b1+nf1, -nf1, nf2, b2+nf1, in[i2]+i1, 1, n1);
		c2 = fir2(-nf1, nf2, b1+nf1, -nf1, nf2, b2+nf1, in[i2]+i1, -1, -n1);
		out[i2][i1] = c2 - c1;
	}
}

void opwdpd(int n1, int n2, float **in, 
	sf_complex **p, float **out, int id)
/*< partial derivative filter of circle interpolating PWD >*/
{
	int i1, i2;
	float p1, p2, c1, z1;

	for(i2=nf1; i2<n2-nf2; i2++)
	for(i1=nf1; i1<n1-nf2; i1++)
	{
		p1 = r*cimagf(p[i2][i1]);
		p2 = r*crealf(p[i2][i1]);
		if(id==0)
		{
			lphpoly_coef(nf1+nf2, c, p1, b1, nf1+nf2, true);
			lphpoly_coef(nf1+nf2, c, p2, b2, nf1+nf2, false);
		} else{
			lphpoly_coef(nf1+nf2, c, p1, b1, nf1+nf2, false);
			lphpoly_coef(nf1+nf2, c, p2, b2, nf1+nf2, true);
		}
		z1 = fir2(-nf1, nf2, b1+nf1, -nf1, nf2, b2+nf1, in[i2]+i1, 1, n1);
		c1 = fir2(-nf1, nf2, b1+nf1, -nf1, nf2, b2+nf1, in[i2]+i1, -1, -n1);
		out[i2][i1] = (c1 - z1);/*(id==0?p2:-p1);*/
	}
}



void opwd_freq(sf_complex dip, int nk, sf_complex**out, bool iir)
/*< frequency response of circle-interpolating PWD >*/
{
	int i1, i2;
	sf_complex c1, c2, z1, z2;

	lphpoly_coef(nf1+nf2, c, r*cimagf(dip), b1, nf1+nf2, false);
	lphpoly_coef(nf1+nf2, c, r*crealf(dip), b2, nf1+nf2, false);

	for(i2=-nk; i2<=nk; i2++)
	{
		c2 = fir_freq(-nf1, nf2, b2+nf1, 0.5*i2/nk);
		for(i1=-nk; i1<=nk; i1++)
		{
			c1 = fir_freq(-nf1, nf2, b1+nf1, 0.5*i1/nk);
			z1 = conjf(c1);	z2 = conjf(c2);
#ifdef SF_HAS_COMPLEX_H
			if(iir)	out[i2+nk][i1+nk] = 1.0 - c1*c2/(z1*z2);
			else	out[i2+nk][i1+nk] = z1*z2 - c1*c2;
#else
			if(iir)	out[i2+nk][i1+nk] = sf_cadd(sf_cmplx(1.0,0.0),sf_cneg(sf_cdiv(sf_cmul(c1,c2),sf_cmul(z1,z2))));
			else	out[i2+nk][i1+nk] = sf_cadd(sf_cmul(z1,z2),sf_cneg(sf_cmul(c1,c2)));
#endif
		}
	}

}



void opwd_close()
/*< release memory >*/
{
	free(b1);
	free(b2);
	free(c[0]);
	free(c);
}









static int n1, n2, nf1, nf2;
static float **u1, **u2, **u3, **u4, **u5, r, p0, eps;
static sf_complex **p;
static bool verb, use_divn;

void odip_init(char* interp, int mf1, int mf2, float rad,
	       int m1, int m2, int *rect, int niter, float dip0, bool vb)
/*< initialize >*/
{
    int n, nn[2];
    nf1 = mf1;
    nf2 = mf2;
    n1 = m1;
    n2 = m2;
    verb = vb;

    u1 = sf_floatalloc2(n1, n2);
    u2 = sf_floatalloc2(n1, n2);
    u3 = sf_floatalloc2(n1, n2);
    u4 = sf_floatalloc2(n1, n2);
    u5 = sf_floatalloc2(n1, n2);
    p = sf_complexalloc2(n1, n2);

    r=rad;
    p0 = dip0;
	
    opwd_init(nf1, nf2, interp, r);
    if(rect[0]>0 && rect[1]>0)
    {
	n = n1*n2;
	nn[0] = n1;
	nn[1] = n2;
	sf_divn_init (2, n, nn, rect, niter, false);
	use_divn=true;
    }else 	use_divn=false;
}

void odip_close()
/*< release memory >*/
{
    free(u1[0]);
    free(u2[0]);
    free(u3[0]);
    free(u4[0]);
    free(u5[0]);
    free(u1);
    free(u2);
    free(u3);
    free(u4);
    free(u5);
    free(p[0]);
    free(p);
    opwd_close();
    if(use_divn)	sf_divn_close();
}

#define divn(a, b)  (a*b/(b*b+10E-15))

static void odip_verb(int it, int n, float *u)
{
    int i1;
    double norm;

    if(verb)
    {
	for(i1=0, norm=0.0; i1<n; i1++)
	    norm += (u[i1]*u[i1]);
	sf_warning("res1 %d %g", it, sqrtf(norm/n));
    }
}

void odip(float **in, float **dip, int nit, float eta)
/*< omnidirectional dip estimation >*/
{
    int it, i1;
    double  c1;
    sf_complex dip0;

#ifdef SF_HAS_COMPLEX_H
    dip0 = r*cexpf(sf_cmplx(0.0, p0));
#else
    dip0 = sf_crmul(cexpf(sf_cmplx(0.0, p0)),r);
#endif

    for(i1=0; i1<n1*n2; i1++)
    {
	dip[0][i1] = p0;
	p[0][i1] = dip0;
    }

    for (it=0; it<nit; it++)
    {
	opwd(n1, n2, in, p, u1);
	opwdpd(n1, n2, in, p, u2, 0);
	opwdpd(n1, n2, in, p, u3, 1);

	odip_verb(it+1, n1*n2, u1[0]);
	for(i1=0, c1=0.0; i1<n1*n2; i1++)
	{
	    u4[0][i1] = u2[0][i1]*crealf(p[0][i1]) - u3[0][i1]*cimagf(p[0][i1]);
	    c1 += (u4[0][i1]*u4[0][i1]);
	}
	c1=sqrtf(c1/(n1*n2));
	for(i1=0; i1<n1*n2; i1++)
	{
	    u1[0][i1] /= c1;
	    u4[0][i1] /= c1;
	}

	if(use_divn)
	{
	    sf_divn(u1[0], u4[0], u5[0]);
	}else{
	    for(i1=0; i1<n1*n2; i1++)
	    {
		u5[0][i1] = divn(u1[0][i1], u4[0][i1]);
	    }
	}
	for(i1=0; i1<n1*n2; i1++)
	{
	    dip[0][i1] -= eta * u5[0][i1];
#ifdef SF_HAS_COMPLEX_H
	    p[0][i1] = r*cexpf(sf_cmplx(0, dip[0][i1]));
#else
	    p[0][i1] = sf_crmul(cexpf(sf_cmplx(0, dip[0][i1])),r);
#endif
	}

    }
}

void oslope(float **in, float **dip, int nit, float eta)
/*< omnidirectional slope estimation >*/
{
    int it, i1;
    double  s1, c1;
    sf_complex dip0;

#ifdef SF_HAS_COMPLEX_H
    dip0 = r*cexpf(sf_cmplx(0, p0));
#else
    dip0 = sf_crmul(cexpf(sf_cmplx(0, p0)),r);
#endif
    for(i1=0; i1<n1*n2; i1++)
	p[0][i1] = dip0;

    for (it=0; it<nit; it++)
    {
	opwd(n1, n2, in, p, u1);
	opwdpd(n1, n2, in, p, u2, 0);
	opwdpd(n1, n2, in, p, u3, 1);

	odip_verb(it+1, n1*n2, u1[0]);

	if(use_divn)
	{
	    for(i1=0, c1=0.0, s1=0.0; i1<n1*n2; i1++)
	    {
		c1 += (u2[0][i1]*u2[0][i1]);
		s1 += (u3[0][i1]*u3[0][i1]);
	    }
	    c1=sqrtf(c1/(n1*n2));
	    s1=sqrtf(s1/(n1*n2));
	    for(i1=0; i1<n1*n2; i1++)
	    {
		u1[0][i1] /= c1;
		u2[0][i1] /= c1;
	    }
	    sf_divn(u1[0], u2[0], u4[0]);
	    for(i1=0; i1<n1*n2; i1++)
	    {
		u1[0][i1] *= c1/s1;
		u3[0][i1] /= s1;
	    }
	    sf_divn(u1[0], u3[0], u5[0]);
	}else{
	    for(i1=0; i1<n1*n2; i1++)
	    {
		u4[0][i1] = divn(u1[0][i1], u2[0][i1]);
		u5[0][i1] = divn(u1[0][i1], u3[0][i1]);
	    }
	}
	for(i1=0; i1<n1*n2; i1++)
	{
#ifdef SF_HAS_COMPLEX_H
	    p[0][i1] -= eta * sf_cmplx(u5[0][i1], u4[0][i1]);
	    p[0][i1] = p[0][i1]*r/(cabsf(p[0][i1])+ 1E-15);
#else
	    p[0][i1] = sf_cadd(p[0][i1],sf_crmul(sf_cmplx(u5[0][i1], u4[0][i1]),-eta));
	    p[0][i1] = sf_crmul(p[0][i1],r/(cabsf(p[0][i1])+ 1E-15));
#endif
	}

    }
    for(i1=0; i1<n1*n2; i1++)
	dip[0][i1] = atan(tan(cargf(p[0][i1])));
}

