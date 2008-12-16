#include <rsf.h>
/*^*/

#include "slinint.h"

float a, b;

void slinint_init(float t)
/*< initialize with shift >*/
{
    a = 1./(1.-t);
    b = t;
}

void slinint_pre(int n, float *x)
/*< prefilter x[n] for interpolation >*/
{
    int i;
    float y0, y1;

    y1 = 0.;
    for (i=0; i < n; i++) { 
	y0 = a*(x[i] - b*y1);
	x[i] = y1 = y0;
    }
}

bool slinint_int(int n, float *c, float x, float *y)
/*< interpolate c[n] to find y(x), return true on success >*/
{
    int i;

    x -= b;
    i = floorf(x);

    if (i < 0 || i > n-1) return false;

    x -= i;

    if (x == 0.) {
	*y = c[i];
	return true;
    }

    if (i == n-1) return false;

    *y = (1.-x)*c[i]+x*c[i+1];
    return true;
}




    
