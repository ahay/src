#include <math.h>

#include "c99.h"
#include "error.h"

float sf_crealf(/*@unused@*/ float complex c) 
{ 
    sf_warning("No support for complex types!!!\n"
	       "Please use a C99-compliant compiler");
    return c; 
}

float sf_cimagf(/*@unused@*/ float complex c) { 
    sf_warning("No support for complex types!!!\n"
	       "Please use a C99-compliant compiler");
    return 0.; 
} 

float sf_cargf(/*@unused@*/ float complex c) 
{ 
    sf_warning("No support for complex types!!!\n"
	       "Please use a C99-compliant compiler");
    return c; 
}

float complex sf_conjf(/*@unused@*/ float complex c) 
{ 
    sf_warning("No support for complex types!!!\n"
	       "Please use a C99-compliant compiler");
    return c; 
}

float sqrtf(float x) { return (float) sqrt(x);}
float logf(float x) { return (float) log(x);}
float log10f(float x) { return (float) log10(x);}
float expf(float x) { return (float) exp(x);}
float fabsf(float x) { return (float) fabs(x);}
float floorf(float x) { return (float) floor(x);}
float ceilf(float x) { return (float) ceil(x);}
float cosf(float x) { return (float) cos(x);}
float sinf(float x) { return (float) sin(x);}
float tanf(float x) { return (float) tan(x);}
float acosf(float x) { return (float) acos(x);}
float asinf(float x) { return (float) asin(x);}
float atanf(float x) { return (float) atan(x);}
float atan2f(float x, float y) { return (float) atan2(x,y);}
float powf(float x, float y) { return (float) pow(x,y);}
float hypotf(float x, float y) { return (float) hypot(x,y);}
float fmodf(float x, float y) { return (float) fmod(x,y);}
