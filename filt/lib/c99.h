#ifndef _sf_c99_h
#define _sf_c99_h

#if defined(__STDC__) && (__STDC_VERSION__ >= 199901L)

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
float floorf(float);
float ceilf(float);
float cosf(float);
float sinf(float);
float tanf(float);
float acosf(float);
float asinf(float);
float atanf(float);
float atan2f(float,float);
float powf(float,float);
float hypotf(float,float);
float fmodf(float,float);

/* What do we do with complex? */
#define complex  
#define I 0.0
#define csqrtf sqrtf
#define clogf logf
#define cexpf expf
#define cabsf fabsf
#define crealf sf_crealf
#define cimagf sf_cimagf
#define conjf sf_conjf
#define cargf sf_cargf

float sf_crealf(float complex c); 
float sf_cimagf(float complex c);
float sf_argf(float complex c);
float complex sf_conjf(float complex c);

#endif

void cprint (float complex c);

#endif

