#include <rsf.h>

#include "radon.h"
#include "cadj.h"

static float dp, p0;
static int nx, np;
static float complex *c0, *dc, czero;

void radon_init (int nx_in, int np_in, float dp_in, float p0_in) 
{
    nx = nx_in; np = np_in;
    dp = dp_in; p0 = p0_in;
    c0 = sf_complexalloc(nx);
    dc = sf_complexalloc(nx);
    czero = 0.;
}

void radon_close () 
{
    free (c0);
    free (dc);
}

void radon_set (float w, float* xx)
{
    int ix;

    for (ix=0; ix < nx; ix++) {
	dc[ix] = cexpf(w*dp*xx[ix]*I);
	c0[ix] = cexpf(w*p0*xx[ix]*I); 
    }
}

void radon_toep (float complex *qq, float eps)
{
    int ix, ip;
    float complex c, q;
    
    qq[0] = eps*eps;
    for (ip=1; ip < np; ip++) {
	qq[ip] = czero;
    }

    for (ix=0; ix < nx; ix++) {
	c = conjf(dc[ix]);
	q = 1.;
	for (ip=0; ip < np-1; ip++) {
	    qq[ip] += q;
	    q *= c;
	}
	qq[np-1] += q;
    }
}

void radon_lop (bool adj, bool add, int nm, int nd, 
		float complex *mm, float complex *dd)
{
    int ix, ip;
    float complex c, d;

    if (nm != np || nd != nx) 
	sf_error("%s: mismatched data sizes",__FILE__);

    adjnull(adj, add, nm, nd, mm, dd);
    
    for (ix=0; ix < nx; ix++) {
	if (adj == true) {
	    c = dc[ix];
	    d = dd[ix]*c0[ix];
	    for (ip=0; ip < np-1; ip++) {
		mm[ip] += d;
		d *= c; 
	    }
	    mm[np-1] += d;
	} else {
	    c = conjf(dc[ix]);
	    d = mm[np-1];
	    for (ip=np-2; ip >= 0; ip--) {
		d = d*c + mm[ip];
	    }
	    dd[ix] += d*conjf(c0[ix]);
	}
    }
}

/* 	$Id: radon.c,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */

