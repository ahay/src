#include <math.h>

#include <rsf.h>

#include "symes.h"

float symes_vel(void* par, float* x)
{
    float v;
    v = 1. + 0.2*sinf(3.*SF_PI*x[1])*sinf(0.5*SF_PI*x[0]);
    v = 1./(v*v);
    return v;
}

void symes_vgrad(void* par, float* x, float* grad) {
    int i;
    grad[1] = 0.6*SF_PI*cosf(3.*SF_PI*x[1])*sinf(0.5*SF_PI*x[0]);
    grad[0] = 0.1*SF_PI*sinf(3.*SF_PI*x[1])*cosf(0.5*SF_PI*x[0]);
    for (i=0; i < 2; i++) {
	grad[i] *= -symes_vel(par, x)/
	    (1. + 0.2*sinf(3.*SF_PI*x[1])*sinf(0.5*SF_PI*x[0]));
    }
}

/* 	$Id$	 */
