#include <stdlib.h>
#include <math.h>

#include "randn.h"

static const float s = 0.449871, tt = -0.386595;
static const float a = 0.19600, b = 0.25472;
static const float r1 = 0.27597, r2 = 0.27846;
  
float random_one (void)
{
    return ((float) rand())/RAND_MAX;
}

float randn_one (void)
{
    float u, v, x, y, q;
    
    do {
	u = random_one ();
	v = random_one ();
	v = 1.7156 * (v - 0.5);
	/*     Evaluate the quadratic form */
	x = u - s;
	y = fabsf(v) - tt;
	q = x*x + y*(a*y - b*x);
	/*     Accept P if inside inner ellipse */
	if (q < r1) break;
	/*     Reject P if outside outer ellipse */
	if (q > r2) continue;
	/*     Reject P if outside acceptance region */
    } while (v*v >= -4.0*logf(u)*u*u);
    return (v/u);
}

void randn (int nr, float *r)
{
    int i;

    for (i = 0; i < nr; i++) {
	r[i] = randn_one ();
    }
}

void random0 (int nr, float *r)
{
    int i;

    for (i = 0; i < nr; i++) {
	r[i] = random_one ();
    }
}
