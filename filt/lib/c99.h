#ifndef _sf_c99_h
#define _sf_c99_h

#include <math.h>

#ifndef __cplusplus

#if !defined(__APPLE__) && defined(__STDC__) && (__STDC_VERSION__ >= 199901L)

/* The following from C99 - must define for C90 */
#include <stdbool.h>       /* define bool, true, false */
#include <complex.h>

#else

typedef enum {false, true} bool;

/* What do we do with float functions? */
float sqrtf(float); 
float logf(float);
float log10f(float);
float expf(float);
float fabsf(float);

#ifndef __APPLE__
float floorf(float);
float ceilf(float);
float fmodf(float,float);
#endif

float cosf(float);
float sinf(float);
float tanf(float);
float acosf(float);
float asinf(float);
float atanf(float);
float atan2f(float,float);
float powf(float,float);
float hypotf(float,float);


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

float sf_crealf(float complex c); 
float sf_cimagf(float complex c);
float sf_cargf(float complex c);
float complex sf_conjf(float complex c);
double sf_creal(double complex c); 
double sf_cimag(double complex c);
double complex sf_conj(double complex c);

#endif

void cprint (float complex c);

#endif /* c++ */

#endif

/* 	$Id: c99.h,v 1.8 2004/06/15 16:27:42 fomels Exp $	 */
