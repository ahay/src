#ifndef _sf_c99_h
#define _sf_c99_h

/* Is it 1994 or 1999?? */
#if defined(__STDC__) && (__STDC_VERSION__ >= 199401L)

/* The following from C99 - must define for C90 */
#include <stdbool.h>       /* define bool, true, false */
#include <complex.h>

#else

typedef enum {false, true} bool;

/* What do we do with float functions? */
/* #define sqrtf sqrt 
#define logf log
#define expf exp
#define fabsf fabs
*/

/* What do we do with complex? */
#define complex  
#define I 0.0
#define csqrtf sqrtf
#define clogf logf
#define cexpf expf
#define cabsf fabsf
#define crealf sf_crealf
#define cimagf sf_cimagf

float sf_crealf(float complex c); 
float sf_cimagf(float complex c);

#endif

#endif
