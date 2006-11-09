/* Anisotropic diffusion, 3-D */

#include <math.h>

#include <rsf.h>
/*^*/

#include "impl3.h"

static float t1, t2, t3, ***y, ***z, ***w, ***ww, ***t, *dist;
static float *tmp1, *tmp2, *tmp3, *d1, *d2, *d3, *w1, *w2, *w3; 
static int nstep, n1, n2, n3, n, nclip, ns, nsnap;
static bool up, verb;
static sf_tris slv1, slv2, slv3;
static sf_file snap;

void impl3_init (float r1, float r2, float r3    /* radius */, 
		 int n1_in, int n2_in, int n3_in /* data size */, 
		 float tau            /* duration */, 
		 float pclip          /* percentage clip */, 
		 bool up_in           /* weighting case */,
		 bool verb_in         /* verbosity flag */,
		 float *dist_in       /* optional distance function */,
		 int nsnap_in         /* number of snapshots */,
		 sf_file snap_in      /* snapshot file */)
/*< initialize >*/
{
    int i;

    t1 = (r1*r1-1.)/12.;
    t2 = (r2*r2-1.)/12.;
    t3 = (r3*r3-1.)/12.;

    nstep = SF_MAX(SF_MAX(t1,t2),t3)/tau;
    if (nstep > 1) {
	t1 /= nstep;
	t2 /= nstep;
	t3 /= nstep;
    } else {
	nstep = 1;
    }
    
    n1 = n1_in;
    n2 = n2_in;
    n3 = n3_in;
    up = up_in;
    verb = verb_in;
    dist = dist_in;
    nsnap = nsnap_in;
    snap = snap_in;
    n = n1*n2*n3;

    ns = SF_MAX(nstep/nsnap,1);
    
    y = sf_floatalloc3(n1,n2,n3);
    z = sf_floatalloc3(n1,n2,n3);
    w = sf_floatalloc3(n1,n2,n3);
    t = sf_floatalloc3(n1,n2,n3);
    tmp1 = sf_floatalloc(n1);
    tmp2 = sf_floatalloc(n2);
    tmp3 = sf_floatalloc(n3);
    ww = sf_floatalloc3(n1,n2,n3);

    d1 = sf_floatalloc(n1);
    w1 = sf_floatalloc(n1);
    d2 = sf_floatalloc(n2);
    w2 = sf_floatalloc(n2);
    d3 = sf_floatalloc(n3);
    w3 = sf_floatalloc(n3);

    slv1 = sf_tridiagonal_init (n1);
    slv2 = sf_tridiagonal_init (n2);
    slv3 = sf_tridiagonal_init (n3);

    nclip = (int) n*pclip*0.01;
    if (nclip < 1) {
	nclip = 1;
    } else if (nclip > n) {
	nclip = n;
    }
    nclip--;

    sf_warning("%s: nstep=%d tau=(%g,%g,%g) nclip=%d",
	       __FILE__,nstep,t1,t2,t3,nclip);

    for (i=0; i < n; i++) {
	w[0][0][i] = 1.;
	ww[0][0][i] = 1.;
    }
}

void impl3_close (void)
/*< free allocated storage >*/
{
    free(**y);
    free(*y);
    free(y);
    free(**z);
    free(*z);
    free(z);
    free(**w);
    free(*w);
    free(w);
    free(tmp1);
    free(tmp2);
    free(tmp3);
    free(**ww);
    free(*ww);
    free(ww);
    free(d1);
    free(w1);
    free(d2);
    free(w2);
    free(d3);
    free(w3);

    sf_tridiagonal_close (slv1);
    sf_tridiagonal_close (slv2);
    sf_tridiagonal_close (slv3);
}

void impl3_set(float *** x)
/*< compute weighting function >*/
{
    int i;
    float a, xsum, wsum, *tmp;

    sf_sobel32(n1,n2,n3,x,w);
    tmp = ww[0][0];
    
    for (i=0; i < n; i++) {
	tmp[i] = w[0][0][i];
    }

    a = sf_quantile(nclip,n,tmp);
    if (a==0.) sf_error("%s: clip at nclip=%d is zero, use a higher pclip",
			__FILE__,nclip);

    for (i=0; i < n; i++) {
	w[0][0][i] = sqrtf(1.+w[0][0][i]/a);
	if (NULL != dist) w[0][0][i] *= dist[i];	
	tmp[i] = 1./w[0][0][i];
	if (up) w[0][0][i] = tmp[i];
    }

    if (verb) {
	wsum = xsum = 0.;
	for (i=0; i < n; i++) {
	    wsum += w[0][0][i];
	    xsum += x[0][0][i]*x[0][0][i];
	}

	sf_warning("xsum=%g, wsum=%g, a=%g", xsum, wsum, a);
    }
}

void impl3_apply (float ***x, bool set, bool adj)
/*< apply diffusion >*/
{
    int istep, i1, i2, i3, i, is;

    is=0;
    for (istep=0; istep < nstep; istep++) {
	if (NULL != snap && 0==istep%ns && is < nsnap) {
	    sf_floatwrite(x[0][0],n,snap);
	    is++;
	}

	if (set) impl3_set(x);

	for (i=0; i < n; i++) {
	    if (!adj) x[0][0][i] *= w[0][0][i];
	    y[0][0][i] = x[0][0][i];
	    z[0][0][i] = x[0][0][i];
	}

	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    w1[i1] = -3.*t1*ww[i3][i2][i1];
		    tmp1[i1] = x[i3][i2][i1];
		    d1[i1] = w[i3][i2][i1];
		}
		d1[0] -= w1[0];
		for (i1=1; i1 < n1-1; i1++) {
		    d1[i1] -= w1[i1] + w1[i1-1];
		}
		d1[n1-1] -= w1[n1-2];
		sf_tridiagonal_define (slv1, d1, w1);
		sf_tridiagonal_solve (slv1, tmp1);
		for (i1=0; i1 < n1; i1++) {
		    x[i3][i2][i1] = tmp1[i1];
		}
	    }
	}
	
	for (i3=0; i3 < n3; i3++) {
	    for (i1=0; i1 < n1; i1++) {
		for (i2=0; i2 < n2; i2++) {
		    w2[i2] = -3.*t2*ww[i3][i2][i1];
		    tmp2[i2] = y[i3][i2][i1];
		    d2[i2] = w[i3][i2][i1];
		}
		d2[0] -= w2[0];
		for (i2=1; i2 < n2-1; i2++) {
		    d2[i2] -= w2[i2] + w2[i2-1];
		}
		d2[n2-1] -= w2[n2-2];
		sf_tridiagonal_define (slv2, d2, w2);
		sf_tridiagonal_solve (slv2, tmp2);
		for (i2=0; i2 < n2; i2++) {
		    y[i3][i2][i1] = tmp2[i2];
		}
	    }
	}

	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		for (i3=0; i3 < n3; i3++) {
		    w3[i3] = -3.*t3*ww[i3][i2][i1];
		    tmp3[i3] = z[i3][i2][i1];
		    d3[i3] = w[i3][i2][i1];
		}
		d3[0] -= w3[0];
		for (i3=1; i3 < n3-1; i3++) {
		    d3[i3] -= w3[i3] + w3[i3-1];
		}
		d3[n3-1] -= w3[n3-2];
		sf_tridiagonal_define (slv3, d3, w3);
		sf_tridiagonal_solve (slv3, tmp3);
		for (i3=0; i3 < n3; i3++) {
		    z[i3][i2][i1] = tmp3[i3];
		}
	    }
	}

	for (i=0; i < n; i++) {
	    x[0][0][i] = (x[0][0][i]+y[0][0][i]+z[0][0][i])/3;
	    if (adj) x[0][0][i] *= w[0][0][i];
	}
    }
}

/* 	$Id: impl3.c 2348 2006-11-02 01:53:52Z sfomel $	 */

