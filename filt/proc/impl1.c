#include <math.h>

#include <rsf.h>

#include "impl1.h"
#include "edge.h"
#include "tridiagonal.h"
#include "quantile.h"

static float t, *w, *d1, *w1;
static int nstep, n, nclip;
static bool up;
static tris slv;

void impl1_init (float r, int n1, float tau, float pclip, bool up_in)
{
    float q;

    q = 0.25/SF_PI;
    t = r*r*q;

    nstep = t/tau;
    if (nstep > 1) {
	t /= nstep;
    } else {
	nstep = 1;
    }

    n = n1;

    w = sf_floatalloc(n);
    d1 = sf_floatalloc(n);
    w1 = sf_floatalloc(n);

    slv = tridiagonal_init (n);

    nclip = (int) n*pclip*0.01;
    if (nclip < 1) {
	nclip = 1;
    } else if (nclip > n) {
	nclip = n;
    }
    nclip--;

    sf_warning("%s: nstep=%d tau=%g nclip=%d",__FILE__,nstep,t,nclip);
}

void impl1_close (void)
{
    free(w);
    free(d1);
    free(w1);

    tridiagonal_close (slv);
}

void impl1_apply (float *x)
{
    int istep, i;
    float a, xsum, wsum;

    for (istep=0; istep < nstep; istep++) {
	grad2(n,x,w);

	for (i=0; i < n; i++) {
	    w1[i] = w[i];
	}

	a = quantile(nclip,n,w1);
	
	if (a==0.) sf_error("%s: clip at nclip=%d is zero, use a higher pclip",
			    __FILE__,nclip);

	wsum = xsum = 0.;
	for (i=0; i < n; i++) {
	    w[i] = sqrtf(1.+w[i]/a);
	    wsum += w[i];
	    xsum += x[i]*x[i];
	}

	sf_warning("step %d of %d, xsum=%g, wsum=%g, a=%g",  
		   istep, nstep, xsum, wsum, a);

	for (i=0; i < n; i++) {
	    if (up) w[i] = 1./w[i];
	    x[i] *= w[i];
	}

	for (i=0; i < n; i++) {
	    d1[i] = w[i];
	    if (up) {
		w1[i] = -t*d1[i];
	    } else {
		w1[i] = -t/d1[i];
	    }
	}

	d1[0] -= w1[0];
	for (i=1; i < n-1; i++) {
	    d1[i] -= w1[i] + w1[i-1];
	}
	d1[n-1] -= w1[n-2];

	tridiagonal_define (slv, d1, w1); 
	tridiagonal_solve (slv, x);
    }
}







