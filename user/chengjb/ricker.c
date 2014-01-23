/*************************************************************************
 * Ricker wavelet 
 * f0: peak frequency
 * t0: time lag
 * A: amplitude
 *
 * Copyright: Tongji University (Jiubing Cheng)
 *  2000.1.2
 * *************************************************************************/
#include <rsf.h>

#include "_cjb.h"

float Ricker(float t, float f0, float t0, float A) 
/*< ricker wavelet:
 * f0: peak frequency
 * t0: time lag
 * A: amplitude
 * ************************>*/
{
        float x=pow(SF_PI*f0*(t-t0),2);
        return -A*exp(-x)*(1-2*x);
}
