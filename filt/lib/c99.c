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
