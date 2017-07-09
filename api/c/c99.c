/* Supplying compatibility with the C99 standard */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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
/*^*/

#include "c99.h"
#include "error.h"
#include "_defs.h"

#ifndef _sf_c99_h

#if !defined (__cplusplus) && !defined(NO_COMPLEX) && defined(__STDC__) && ((__STDC_VERSION__ >= 199901L) || defined(__ICC))
/*^*/

#define SF_HAS_COMPLEX_H
/*^*/

/* The following from C99 - must define for C90 */
#include <complex.h>
#undef I
#ifdef sun
#define I (1.0fi)
#else
#define I _Complex_I
#endif
/*^*/

typedef float complex sf_complex;
typedef double complex sf_double_complex;
/*^*/

#else
/*^*/

#include "kiss_fft.h"
typedef kiss_fft_cpx sf_complex;
typedef struct {
    double r, i;
} sf_double_complex;
/*^*/

#endif
/*^*/

#endif

#ifdef SF_HAS_COMPLEX_H
/*^*/

float complex sf_cmplx(float re, float im)
/*< complex number >*/
{
    float complex c;
#ifdef __GNUC__
    __real__ c = re;
    __imag__ c = im;
#else
    c = re+im*_Complex_I;
#endif
    return c;
}

double complex sf_dcmplx(double re, double im)
/*< complex number >*/
{
    double complex c;
#ifdef __GNUC__
    __real__ c = re;
    __imag__ c = im;
#else
    c = re+im*_Complex_I;
#endif
    return c;
}

#else
/*^*/

kiss_fft_cpx sf_cmplx(float re, float im)
/*< complex number >*/
{
    kiss_fft_cpx c;
    c.r = re;
    c.i = im;
    return c;
}

sf_double_complex sf_dcmplx(double re, double im)
/*< complex number >*/
{
    sf_double_complex c;
    c.r = re;
    c.i = im;
    return c;
}

#endif
/*^*/

#if !defined(__cplusplus) && !defined(SF_HAS_COMPLEX_H)
/*^*/

#if !defined(hpux) && !defined(__hpux)
/*^*/

float copysignf(float x, float y)
/*< float copysign >*/
{ return (float) copysign(x,y);}

#endif
/*^*/

float sqrtf(float x) 
/*< float sqrt >*/
{ return (float) sqrt(x);}

float logf(float x)  
/*< float log >*/
{ return (float) log(x);}

float log10f(float x) 
/*< float log10 >*/
{ return (float) log10(x);}

float expf(float x) 
/*< float exp >*/
{ return (float) exp(x);}

float erff(float x) 
/*< float erf >*/
{ return (float) erf(x);}

float erfcf(float x) 
/*< float erfc >*/
{ return (float) erfc(x);}

#if !defined(hpux) && !defined(__hpux)
/*^*/

float fabsf(float x) 
/*< float fabs >*/
{ return (float) fabs(x);}

#endif
/*^*/

float fmaxf(float x, float y)
/*< float fmax >*/
{ return SF_MAX(x,y);}

float fminf(float x, float y)
/*< float fmin >*/
{ return SF_MIN(x,y);}

float floorf(float x)
/*< float floor >*/
{ return (float) floor(x);}

float ceilf(float x) 
/*< float ceil >*/
{ return (float) ceil(x);}

float roundf(float x) 
/*< round to nearest integer >*/
{ return ((x < 0.0)? ceilf(x-0.5): floorf(x+0.5));}

float fmodf(float x, float y) 
/*< float fmod >*/
{ return (float) fmod(x,y);}

float cosf(float x) 
/*< float cos >*/
{ return (float) cos(x);}

float sinf(float x) 
/*< float sin >*/
{ return (float) sin(x);}

float tanf(float x) 
/*< float tan >*/
{ return (float) tan(x);}

float acosf(float x) 
/*< float acos >*/
{ return (float) acos(x);}

float asinf(float x) 
/*< float asin >*/
{ return (float) asin(x);}

float atanf(float x) 
/*< float atan >*/
{ return (float) atan(x);}

float atan2f(float x, float y) 
/*< float atan2 >*/
{ return (float) atan2(x,y);}

float log2f(float x) 
/*< float log2 >*/
{ return (float) log(x)/log(2.0);}

float coshf(float x) 
/*< float cosh >*/
{ return (float) cosh(x);}

float sinhf(float x) 
/*< float sinh >*/
{ return (float) sinh(x);}

float tanhf(float x) 
/*< float tanh >*/
{ return (float) tanh(x);}

float acoshf(float x) 
/*< float acosh >*/
{ extern double acosh(double x);
return (float) acosh(x);}

float asinhf(float x) 
/*< float asinh >*/
{ extern double asinh(double x);
return (float) asinh(x);}

float atanhf(float x) 
/*< float atanh >*/
{ extern double atanh(double x);
 return (float) atanh(x);}

float powf(float x, float y) 
/*< float pow >*/
{ return (float) pow(x,y);}

float hypotf(float x, float y) 
/*< float hypot >*/
{ extern double hypot(double x, double y);
 return (float) hypot(x,y);}

long lrint(double num)
/*< round to integer >*/
{ return (long)(num < 0.0 ? (num - 0.5) : (num + 0.5)); }

long long llround(double num)
/*< round to integer >*/
{ return (long long)(num < 0.0 ? (num - 0.5) : (num + 0.5)); }

#if defined(hpux) || defined(__hpux)
/*^*/

static char	*digits = "0123456789abcdefghijklmnopqrstuvwxyz";

long long
strtoll(const char *ptr, const char **endptr, int base)
/*< strtoll replacement >*/
{
	const char	*cp;
	long long	 ret;
	char		*dig;
	int		 d;

	for (ret = 0, cp = ptr ; *cp && (dig = strchr(digits, *cp)) != NULL  && (d = (int)(dig - digits)) < base ; cp++) {
		ret = (ret * base) + d;
	}
	if (endptr != NULL) {
		*endptr = cp;
	}
	return ret;
}

unsigned long long
strtoull(const char *ptr, const char **endptr, int base)
/*< strtoull replacement >*/
{
	const char	*cp;
	unsigned long long	 ret;
	char		*dig;
	int		 d;

	for (ret = 0, cp = ptr ; *cp && (dig = strchr(digits, *cp)) != NULL  && (d = (unsigned int)(dig - digits)) < base ; cp++) {
		ret = (ret * base) + d;
	}
	if (endptr != NULL) {
		*endptr = cp;
	}
	return ret;
}

#endif
/*^*/

#ifndef _sf_c99_h

#ifdef sun
extern int finite(double x);
#define isfinite(x) finite(x)
int isinf(double x) { return !finite(x) && x==x; }
#endif
/*^*/

#ifndef _sf_c99_h

#endif
/*^*/

