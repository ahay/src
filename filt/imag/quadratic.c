#include <math.h>

#include "quadratic.h"

/* solves a x^2 + 2 b x + c == 0 for smallest positive x */
float quadratic_solve (float a, float b, float c) 
{
    float d;
    
    if (fabsf(c) < EPS && 
	((b > EPS && a  < -EPS) || (b < - EPS && a > EPS))) 
	return (-2.*b/a);
    
    d = b*b - a*c;
    if (d < 0.) return HUG;
    
    d = sqrtf(d)+fabsf(b);
    
    if (b*c <= 0.) {
	if (d > EPS && fabsf(c) > 0.) return (fabsf(c)/d);
    } else if (b*a <= 0.) {
	if (fabsf(a) > EPS && d > 0.) return (d/fabsf(a));
    } 
    
    return HUG;
}

/* 	$Id: quadratic.c,v 1.2 2003/09/30 14:30:53 fomels Exp $	 */
