/* Anisotropic diffusion, 2-D */

#include <math.h>

#include <rsf.h>
/*^*/

#include "impl2.h"
#include "edge.h"

static float t1, t2, **y, **w, **ww, *tmp, *d1, *d2, *w1, *w2, **t, *dist;
static int nstep, n1, n2, n, nclip, ns, nsnap;
static bool up, verb;
static sf_tris slv1, slv2;
static sf_file snap;

void impl2_init (float r1, float r2   /* radius */, 
		 int n1_in, int n2_in /* data size */, 
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

    nstep = SF_MAX(t1,t2)/tau;
    if (nstep > 1) {
	t1 /= nstep;
	t2 /= nstep;
    } else {
	nstep = 1;
    }

    n1 = n1_in;
    n2 = n2_in;
    up = up_in;
    verb = verb_in;
    dist = dist_in;
    nsnap = nsnap_in;
    snap = snap_in;
    n = n1*n2;

    ns = SF_MAX(nstep/nsnap,1);
    
    y = sf_floatalloc2(n1,n2);
    w = sf_floatalloc2(n1,n2);
    t = sf_floatalloc2(n1,n2);
    tmp = sf_floatalloc(n);
    ww = sf_floatalloc2(n1,n2);

    d1 = sf_floatalloc(n1);
    w1 = sf_floatalloc(n1);
    d2 = sf_floatalloc(n2);
    w2 = sf_floatalloc(n2);

    slv1 = sf_tridiagonal_init (n1);
    slv2 = sf_tridiagonal_init (n2);

    nclip = (int) n*pclip*0.01;
    if (nclip < 1) {
	nclip = 1;
    } else if (nclip > n) {
	nclip = n;
    }
    nclip--;

    sf_warning("%s: nstep=%d tau=(%g,%g) nclip=%d",
	       __FILE__,nstep,t1, t2,nclip);

    for (i=0; i < n; i++) {
	w[0][i] = 1.;
	ww[0][i] = 1.;
    }
}

void impl2_close (void)
/*< free allocated storage >*/
{
    free(*y);
    free(y);
    free(*w);
    free(w);
    free(tmp);
    free(*ww);
    free(ww);
    free(d1);
    free(w1);
    free(d2);
    free(w2);

    sf_tridiagonal_close (slv1);
    sf_tridiagonal_close (slv2);
}

void impl2_set(float ** x)
/*< compute weighting function >*/
{
    int i;
    float a, xsum, wsum;

    sobel2(n1,n2,x,w);

    for (i=0; i < n; i++) {
	tmp[i] = w[0][i];
    }

    a = sf_quantile(nclip,n,tmp);
    if (a==0.) sf_error("%s: clip at nclip=%d is zero, use a higher pclip",
			__FILE__,nclip);

    for (i=0; i < n; i++) {
	w[0][i] = sqrtf(1.+w[0][i]/a);
	if (NULL != dist) w[0][i] *= dist[i];	
	ww[0][i] = 1./w[0][i];
	if (up) w[0][i] = ww[0][i];
    }

    if (verb) {
	wsum = xsum = 0.;
	for (i=0; i < n; i++) {
	    wsum += w[0][i];
	    xsum += x[0][i]*x[0][i];
	}

	sf_warning("xsum=%g, wsum=%g, a=%g", xsum, wsum, a);
    }
}

void impl2_apply (float **x, bool set, bool adj)
/*< apply diffusion >*/
{
    int istep, i1, i2, i, is;

    is=0;
    for (istep=0; istep < nstep; istep++) {
	if (NULL != snap && 0==istep%ns && is < nsnap) {
	    sf_floatwrite(x[0],n,snap);
	    is++;
	}

	if (set) impl2_set(x);

	for (i=0; i < n; i++) {
	    if (!adj) x[0][i] *= w[0][i];
	    y[0][i] = x[0][i];
	}

	for (i1=0; i1 < n1; i1++) {
	    for (i2=0; i2 < n2; i2++) {
		w2[i2] = -t2*ww[i2][i1];
		tmp[i2] = x[i2][i1];
		d2[i2] = w[i2][i1];
	    }
	    d2[0] -= w2[0];
	    for (i2=1; i2 < n2-1; i2++) {
		d2[i2] -= w2[i2] + w2[i2-1];
	    }
	    d2[n2-1] -= w2[n2-2];
	    sf_tridiagonal_define (slv2, d2, w2);
	    sf_tridiagonal_solve (slv2, tmp);
	    for (i2=0; i2 < n2; i2++) {
		x[i2][i1] = tmp[i2];
	    }
	}
	for (i=0; i < n; i++) {
	    x[0][i] *= w[0][i];
	}

	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		w1[i1] = -t1*ww[i2][i1];
		d1[i1] = w[i2][i1];
	    }
	    d1[0] -= w1[0];
	    for (i1=1; i1 < n1-1; i1++) {
		d1[i1] -= w1[i1] + w1[i1-1];
	    }
	    d1[n1-1] -= w1[n1-2];
	    sf_tridiagonal_define (slv1, d1, w1); 
	    sf_tridiagonal_solve (slv1, x[i2]);
	    sf_tridiagonal_solve (slv1, y[i2]);
	}
	for (i=0; i < n; i++) {
	    y[0][i] *= w[0][i];
	}

	for (i1=0; i1 < n1; i1++) {
	    for (i2=0; i2 < n2; i2++) {
		w2[i2] = -t2*ww[i2][i1];
		tmp[i2] = y[i2][i1];
		d2[i2] = w[i2][i1];
	    }
	    d2[0] -= w2[0];
	    for (i2=1; i2 < n2-1; i2++) {
		d2[i2] -= w2[i2] + w2[i2-1];
	    }
	    d2[n2-1] -= w2[n2-2];
	    sf_tridiagonal_define (slv2, d2, w2);
	    sf_tridiagonal_solve (slv2, tmp);
	    for (i2=0; i2 < n2; i2++) {
		y[i2][i1] = tmp[i2];
	    }
	}
	for (i=0; i < n; i++) {
	    x[0][i] = 0.5*(x[0][i]+y[0][i]);
	    if (adj) x[0][i] *= w[0][i];
	}
    }
}

void impl2_lop (bool adj, bool add, int nx, int ny, float* x, float* y)
/*< linear operator >*/
{
    int i1, i2;

    sf_adjnull(adj,add,nx,ny,x,y);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    t[i2][i1] = adj? y[i1+i2*n1]: x[i1+i2*n1];
	}
    }

    impl2_apply(t,false,adj);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    if (adj) {
		x[i1+i2*n1] += t[i2][i1];
	    } else {
		y[i1+i2*n1] += t[i2][i1];
	    }
	}
    }
}

/* 	$Id$	 */

