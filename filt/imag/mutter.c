#include <math.h>

#include <rsf.h>

#include "mutter.h"

static float dt, t0;
static int nt;              

void mutter_init (int n1, float o1, float d1)
{
    nt = n1;
    t0 = o1;
    dt = d1;
}

void mutter (float tp, float slope0, float slopep, float x, float *data)
{
    int it;
    float wt, t;

    x = fabsf(x);

    for (it=0; it < nt; it++) {
	t = t0+it*dt;
	wt = t - x * slope0;
	if (wt < 0.) {
	    data[it] = 0.;
	} else {
	    if (t <= tp + x * slopep) {
		wt = sinf(0.5 * SF_PI * 
			  (t-x*slope0)/(tp+x*(slopep-slope0)));
		data[it] *= (wt*wt);
	    } 
	}
    }
}
