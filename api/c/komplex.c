/* Complex number operations */
/*
  Copyright (C) 2006 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <math.h>

#include "komplex.h"
#include "error.h"
#include "_defs.h"

#include "kiss_fft.h"
#include "c99.h"
/*^*/

#ifndef SF_HAS_COMPLEX_H
/*^*/

#define crealf  sf_crealf
#define creal   sf_creal
#define cimagf  sf_cimagf
#define cimag   sf_cimag
#define conjf   sf_conjf
#define cabsf   sf_cabsf
#define cabs    sf_cabsd
#define cargf   sf_cargf
#define carg    sf_carg
#define ccosf   sf_ccosf
#define csinf   sf_csinf
#define ctanf   sf_ctanf
#define cacosf  sf_cacosf
#define casinf  sf_casinf
#define catanf  sf_catanf
#define ccoshf  sf_ccoshf
#define csinhf  sf_csinhf
#define ctanhf  sf_ctanhf
#define cacoshf sf_cacoshf
#define casinhf sf_casinhf
#define catanhf sf_catanhf
#define cexpf   sf_cexpf
#define clogf   sf_clogf
#define csqrtf  sf_csqrtf
#define cpowf   sf_cpowf
/*^*/

double sf_creal(sf_double_complex c)
/*< real part >*/
{
    return c.r;
}

double sf_cimag(sf_double_complex c)
/*< imaginary part >*/
{
    return c.i;
}

sf_double_complex sf_dcneg(sf_double_complex a)
/*< unary minus >*/
{
    a.r = -a.r;
    a.i = -a.i;
    return a;
}

sf_double_complex sf_dcadd(sf_double_complex a, sf_double_complex b)
/*< complex addition >*/
{
    a.r += b.r;
    a.i += b.i;
    return a;
}

sf_double_complex sf_dcsub(sf_double_complex a, sf_double_complex b)
/*< complex subtraction >*/
{
    a.r -= b.r;
    a.i -= b.i;
    return a;
}

sf_double_complex sf_dcmul(sf_double_complex a, sf_double_complex b)
/*< complex multiplication >*/
{
    sf_double_complex c;
    c.r = a.r*b.r-a.i*b.i;
    c.i = a.i*b.r+a.r*b.i;
    return c;
}

kiss_fft_cpx sf_dccmul(sf_double_complex a, kiss_fft_cpx b)
/*< complex multiplication >*/
{
    kiss_fft_cpx c;
    c.r = a.r*b.r-a.i*b.i;
    c.i = a.i*b.r+a.r*b.i;
    return c;
}

sf_double_complex sf_dcdmul(sf_double_complex a, kiss_fft_cpx b)
/*< complex multiplication >*/
{
    sf_double_complex c;
    c.r = a.r*b.r-a.i*b.i;
    c.i = a.i*b.r+a.r*b.i;
    return c;
}

sf_double_complex sf_dcrmul(sf_double_complex a, double b)
/*< complex by real multiplication >*/
{
    a.r *= b;
    a.i *= b;
    return a;
}

sf_double_complex sf_dcdiv(sf_double_complex a, sf_double_complex b)
/*< complex division >*/
{
    sf_double_complex c;
    double r,den;
    if (fabsf(b.r)>=fabsf(b.i)) {
	r = b.i/b.r;
	den = b.r+r*b.i;
	c.r = (a.r+r*a.i)/den;
	c.i = (a.i-r*a.r)/den;
    } else {
	r = b.r/b.i;
	den = b.i+r*b.r;
	c.r = (a.r*r+a.i)/den;
	c.i = (a.i*r-a.r)/den;
    }
    return c;
}

double sf_carg(sf_double_complex z)
/*< replacement for cargf >*/
{
    extern double atan2(double,double);
    return atan2(z.i,z.r);
}

double sf_cabsd(sf_double_complex z)
/*< replacement for cabs >*/
{
    extern double hypot(double,double);
    return hypot(z.r,z.i);
}

#endif
/*^*/

#if !defined(__cplusplus)
/*^*/

float sf_cabs(sf_complex c)
/*< complex absolute value >*/
{
    return hypotf(crealf(c),cimagf(c));
}\

#endif
/*^*/

float sf_crealf(kiss_fft_cpx c)
/*< real part >*/
{
    return c.r;
}

float sf_cimagf(kiss_fft_cpx c)
/*< imaginary part >*/
{
    return c.i;
}

void cprint (sf_complex c)
/*< print a complex number (for debugging purposes) >*/
{
    sf_warning("%g+%gi",crealf(c),cimagf(c));
}

kiss_fft_cpx sf_cadd(kiss_fft_cpx a, kiss_fft_cpx b)
/*< complex addition >*/
{
    a.r += b.r;
    a.i += b.i;
    return a;
}

kiss_fft_cpx sf_csub(kiss_fft_cpx a, kiss_fft_cpx b)
/*< complex subtraction >*/
{
    a.r -= b.r;
    a.i -= b.i;
    return a;
}

kiss_fft_cpx sf_csqrtf (kiss_fft_cpx c)
/*< complex square root >*/
{

#if !defined(hpux) && !defined(__hpux)
    extern float copysignf(float x, float y);
#endif

    float d, r, s;
    kiss_fft_cpx v;

    if (c.i == 0) {
      if (c.r < 0) {
	  v.r = 0.;
	  v.i = copysignf (sqrtf (-c.r), c.i);
      } else {
	  v.r =  fabsf (sqrtf (c.r));
	  v.i =  copysignf (0, c.i);
      }
    } else if (c.r == 0) {
	r = sqrtf (0.5 * fabsf (c.i));
	v.r = r;
	v.i = copysignf (r, c.i);
    } else {
	d = hypotf (c.r, c.i);
	/* Use the identity   2  Re res  Im res = Im x
	   to avoid cancellation error in  d +/- Re x.  */
	if (c.r > 0) {
	    r = sqrtf (0.5f * d + 0.5f * c.r);
	    s = (0.5f * c.i) / r;
        } else {
	    s = sqrtf (0.5f * d - 0.5f * c.r);
	    r = fabsf ((0.5f * c.i) / s);
        }
	v.r = r;
	v.i = copysignf (s, c.i);
    }
    return v;
}

kiss_fft_cpx sf_cdiv(kiss_fft_cpx a, kiss_fft_cpx b)
/*< complex division >*/
{
    kiss_fft_cpx c;
    float r,den;
    if (fabsf(b.r)>=fabsf(b.i)) {
	r = b.i/b.r;
	den = b.r+r*b.i;
	c.r = (a.r+r*a.i)/den;
	c.i = (a.i-r*a.r)/den;
    } else {
	r = b.r/b.i;
	den = b.i+r*b.r;
	c.r = (a.r*r+a.i)/den;
	c.i = (a.i*r-a.r)/den;
    }
    return c;
}

kiss_fft_cpx sf_cmul(kiss_fft_cpx a, kiss_fft_cpx b)
/*< complex multiplication >*/
{
    kiss_fft_cpx c;
    c.r = a.r*b.r-a.i*b.i;
    c.i = a.i*b.r+a.r*b.i;
    return c;
}

kiss_fft_cpx sf_crmul(kiss_fft_cpx a, float b)
/*< complex by real multiplication >*/
{
    a.r *= b;
    a.i *= b;
    return a;
}

kiss_fft_cpx sf_cneg(kiss_fft_cpx a)
/*< unary minus >*/
{
    a.r = -a.r;
    a.i = -a.i;
    return a;
}


kiss_fft_cpx sf_conjf(kiss_fft_cpx z)
/*< complex conjugate >*/
{
    z.i = -z.i;
    return z;
}

float sf_cabsf(kiss_fft_cpx z)
/*< replacement for cabsf >*/
{
    return hypotf(z.r,z.i);
}


float sf_cargf(kiss_fft_cpx z)
/*< replacement for cargf >*/
{
    return atan2f(z.i,z.r);
}

kiss_fft_cpx sf_ctanhf(kiss_fft_cpx z)
/*< complex hyperbolic tangent >*/
{
    float x, y, d;

    x = z.r;
    y = z.i;

    d = coshf(2*x) + cosf(2*y);
    z.r = sinhf(2*x)/ d;
    z.i = sinf (2*y)/ d;

    return z;
}

kiss_fft_cpx sf_ccosf(kiss_fft_cpx z)
/*< complex cosine >*/
{
    float x, y;

    x = z.r;
    y = z.i;

    z.r = coshf(y)*cosf(x);
    z.i = -sinhf(y)*sinf(x);

    return z;
}

kiss_fft_cpx sf_ccoshf(kiss_fft_cpx z)
/*< complex hyperbolic cosine >*/
{
    float x, y;

    x = z.r;
    y = z.i;

    z.r = coshf(x)*cosf(y);
    z.i = sinhf(x)*sinf(y);

    return z;
}


kiss_fft_cpx sf_csinf(kiss_fft_cpx z)
/*< complex sine >*/
{
    float x, y;

    x = z.r;
    y = z.i;

    z.r = coshf(y)*sinf(x);
    z.i = sinhf(y)*cosf(x);

    return z;
}

kiss_fft_cpx sf_csinhf(kiss_fft_cpx z)
/*< complex hyperbolic sine >*/
{
    float x, y;

    x = z.r;
    y = z.i;

    z.r = sinhf(x)*cosf(y);
    z.i = coshf(x)*sinf(y);

    return z;
}

kiss_fft_cpx sf_clogf(kiss_fft_cpx z)
/*< complex natural logarithm >*/
{
    float x, y;

    x = z.r;
    y = z.i;

    z.r = logf(hypotf(x,y));
    z.i = atan2f(y,x);

    return z;
}

kiss_fft_cpx sf_cexpf(kiss_fft_cpx z)
/*< complex exponential >*/
{
    float x, y;

    x = expf(z.r);
    y = z.i;

    z.r = x*cosf(y);
    z.i = x*sinf(y);

    return z;
}

kiss_fft_cpx sf_ctanf(kiss_fft_cpx z)
/*< complex tangent >*/
{
    return sf_cdiv(sf_csinf(z),sf_ccosf(z));
}

kiss_fft_cpx sf_casinf(kiss_fft_cpx z)
/*< complex hyperbolic arcsine >*/
{
    float x, y;
    kiss_fft_cpx z2;

    x = z.r;
    y = z.i;

    if (0.0 == y) {
	z2.r = asinf(x);
	z2.i = 0.0;
    } else { 
	z.r = 1.0 - (x - y) * (x + y);
	z.i = -2.0 * x * y;
	z = sf_csqrtf(z);

	z.r -= y;
	z.i += x;
	z = sf_clogf(z);

	z2.r =  z.i;
	z2.i = -z.r;
    }
  
    return z2;
}

kiss_fft_cpx sf_cacosf(kiss_fft_cpx z)
/*< complex hyperbolic arccosine >*/
{
    z = sf_casinf(z);
    z.r = SF_PI/2 - z.r;
    z.i = - z.i;
    return z;
}

kiss_fft_cpx sf_catanf(kiss_fft_cpx z)
/*< complex arctangent >*/
{
    kiss_fft_cpx z2;

    z2.r = -z.r;
    z2.i = 1.0-z.i;
    z.i += 1.0;
    
    z2 = sf_clogf(sf_cdiv(z,z2));
    z.r = -0.5f*z2.i;
    z.i = 0.5f*z2.r; 
    /* signs? */

    return z;
}     

kiss_fft_cpx sf_catanhf(kiss_fft_cpx z)
/*< complex hyperbolic arctangent >*/
{
    kiss_fft_cpx z2;

    z2.r = -z.i;
    z2.i =  z.r;
    z2 = sf_catanf(z2);
    z.r =  z2.i;
    z.i = -z2.r; 
    /* signs? */

    return z;
}     

kiss_fft_cpx sf_casinhf(kiss_fft_cpx z)
/*< complex hyperbolic sine >*/
{
    kiss_fft_cpx z2;

    z2.r = -z.i;
    z2.i =  z.r;
    z2 = sf_casinf(z2);
    z.r =  z2.i;
    z.i = -z2.r; 
    /* signs? */

    return z;
}     

kiss_fft_cpx sf_cacoshf(kiss_fft_cpx z)
/*< complex hyperbolic cosine >*/
{
    kiss_fft_cpx z2;

    z2 = sf_casinf(z);
    z.r = z2.i;
    z.i = SF_PI/2.0-z2.r;
    /* signs? */

    return z;
}

kiss_fft_cpx sf_cpowf(kiss_fft_cpx a, kiss_fft_cpx b)
/*< complex power >*/
{
    float i, r, rho, theta;
    kiss_fft_cpx c;

    r = sf_cabsf(a);
    i = sf_cargf (a);

    if (r == 0.0) {
	c.r = 0.0;
	c.i = 0.0;
    } else {
	theta = i * b.r;
 
	if (b.i == 0.0) {
	    rho = powf (r,b.r);
	} else {
	    r = logf(r);

	    theta += r * b.i;
	    rho = expf(r * b.r - i * b.i);
	}

	c.r = rho * cosf (theta);
	c.i = rho * sinf (theta);
    }

    return c;
}

