#ifndef _sf_bool_h
#define _sf_bool_h

/* Standard defines of operators, usable on C90 up */
#if defined(__STDC__) && (__STDC_VERSION__ >= 199901L)
   /* The following from C99 - must define for C90 */
    #include <stdbool.h>       /* define bool, true, false */
    #include <complex.h>
#else
    typedef enum {false=0, true} Bool;
    #define false false
    #define true  true
    #define bool Bool 

/* What do we do with complex? */
    #define complex  
    #define I 0.0
    #define sqrtf sqrt
    #define csqrtf sqrtf
    #define logf log
    #define expf exp
    #define fabsf fabs
    #define clogf logf
    #define cexpf expf
    #define cabsf fabsf
    #define crealf sf_crealf
    #define cimagf sf_cimagf
#endif
    float sf_crealf(float complex c); 
    float sf_cimagf(float complex c); 
#endif
