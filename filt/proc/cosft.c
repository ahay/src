#include <rsf.h>

#include "cosft.h"

static int nt, nw, n;
static float *p /* , dt */;
static float complex *pp;

void cosft_init(int n1 /*, float o1, float d1 */) {
    nt = sf_npfar(2*(n1-1));
    nw = nt/2+1;
    sf_warning("nw=%d",nw);
    p  = sf_floatalloc (nt);
    pp = sf_complexalloc(nw);
    n = n1;
/*    dt = 2.*SF_PI*o1/(nt*d1); */
}

void cosft_close(void) {
    free (p);
    free (pp);
}

void cosft_frw (float *q, int o1, int d1) {
    int i;

    for (i=0; i < n; i++) {
	p[i] = q[o1+i*d1];
    }
    for (i=n; i <= nt/2; i++) {
	p[i]=0.0;
    }
    for (i=nt/2+1; i < nt; i++) {
	p[i] = p[nt-i];
    }
    
    sf_pfarc(1,nt,p,pp);

    for (i=0; i < n; i++) {
	q[o1+i*d1] = crealf(pp[i]);
    }
}

void cosft_inv (float *q, int o1, int d1) {
    int i;

    for (i=0; i < n; i++) {
	pp[i] = q[o1+i*d1];
    }
    for (i=n; i < nw; i++) {
	pp[i] = 0.;
    }

/*

    if (0. != dt) {
	for (i=0; i < n; i++) {
	    pp[i] *= cexpf(I*i*dt);
	}
    }
  
*/
  
    sf_pfacr(-1,nt,pp,p);

    for (i=0; i < n; i++) {
	q[o1+i*d1] = p[i]/nt;
    }
}

/* 	$Id: cosft.c,v 1.6 2003/12/04 05:13:21 fomels Exp $	 */
