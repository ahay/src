#ifndef _sf_c99_h
#define _sf_c99_h

#ifndef __cplusplus

#include <math.h>

#if !defined(__APPLE__) && !defined(__CYGWIN__) && defined(__STDC__) && (__STDC_VERSION__ >= 199901L)

/* The following from C99 - must define for C90 */
#include <stdbool.h>       /* define bool, true, false */
#include <complex.h>

#else

typedef enum {false, true} bool;

/* What do we do with complex? */
#define complex  
#define I 0.0
#define csqrtf sqrtf
#define clogf logf
#define clog log
#define cexpf expf
#define cexp exp
#define cabsf fabsf
#define crealf sf_crealf
#define cimagf sf_cimagf
#define conjf sf_conjf
#define cargf sf_cargf
#define creal sf_creal
#define cimag sf_cimag
#define conj sf_conj

#endif

void cprint (float complex c);
/*< print a complex number (for debugging purposes) >*/

#if defined(__APPLE__) || defined(__CYGWIN__) || defined (__STDC__) || (__STDC_VERSION__ < 199901L)

float sf_crealf(/*@unused@*/ float complex c) ;
/*< real part of a complex number >*/

float sf_cimagf(/*@unused@*/ float complex c) ;
/*< imaginary part of a complex number >*/

double sf_creal(/*@unused@*/ double complex c) ;
/*< real part of a complex number >*/

double sf_cimag(/*@unused@*/ double complex c) ;
/*< imaginary part of a complex number >*/


float sf_cargf(/*@unused@*/ float complex c);
/*< phase of a complex number >*/

float complex sf_conjf(/*@unused@*/ float complex c);
/*< conjugate of a complex number >*/

double complex sf_conj(/*@unused@*/ double complex c) ;
/*< conjugate of a complex number >*/

float sqrtf(float x) ;
/*< float sqrt >*/

float logf(float x)  ;
/*< float log >*/

float log10f(float x) ;
/*< float log10 >*/

float expf(float x) ;
/*< float exp >*/

float fabsf(float x) ;
/*< float fabs >*/

#ifndef __APPLE__

float floorf(float x);
/*< float floor >*/

float ceilf(float x) ;
/*< float ceil >*/

float fmodf(float x, float y) ;
/*< float fmod >*/

#endif

float cosf(float x) ;
/*< float cos >*/

float sinf(float x) ;
/*< float sin >*/

float acosf(float x) ;
/*< float acos >*/

float asinf(float x) ;
/*< float asin >*/

float atanf(float x) ;
/*< float atan >*/

float atan2f(float x, float y) ;
/*< float atan2 >*/

float powf(float x, float y) ;
/*< float pow >*/

float hypotf(float x, float y) ;
/*< float hypot >*/

#endif

#endif /* c++ */

#endif
