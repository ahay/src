#include <math.h>

#include <rsf.h>

#include "mutter.h"

static float dt, t0;
static int nt;
static bool abs0, inner;              

void mutter_init (int n1, float o1, float d1, bool abs1, bool inner1)
{
    nt = n1;
    t0 = o1;
    dt = d1;
    abs0 = abs1;
    inner = inner1;
}

void mutter (float tp, float slope0, float slopep, float x, float *data)
{
    int it;
    float wt, t;

    if (abs0) x = fabsf(x);

    for (it=0; it < nt; it++) {
	t = t0+it*dt;
	wt = t - x * slope0;
	if ((inner && wt > 0.) || (!inner && wt < 0.)) {
	    data[it] = 0.;
	} else {
	    wt = t - tp - x * slopep;
	    if ((inner && wt >=0.) || (!inner && wt <= 0.)) {
		wt = sinf(0.5 * SF_PI * 
			  (t-x*slope0)/(tp+x*(slopep-slope0)));
		data[it] *= (wt*wt);
	    } 
	}
    }
}

/* 	$Id: mutter.c,v 1.4 2004/06/15 16:27:42 fomels Exp $	 */

