#include <rsf.h>

#include "cosft.h"

static int nt, nw;
static float *p /* , dt */;
static float complex *pp;
static kiss_fftr_cfg forw, invs;

void cosft_init(int nw /*, float o1, float d1 */) {
    nt = 2*(nw-1);
    sf_warning("nw=%d",nw);
    p  = sf_floatalloc (nt);
    pp = sf_complexalloc(nw);
    forw = kiss_fftr_alloc(nt,0,NULL,NULL);
    invs = kiss_fftr_alloc(nt,1,NULL,NULL);
}

void cosft_close(void) {
    free (p);
    free (pp);
    free (forw);
    free (invs);
}

void cosft_frw (float *q, int o1, int d1) {
    int i;

    for (i=0; i < nw; i++) {
	p[i] = q[o1+i*d1];
    }
    for (i=nw; i < nt; i++) {
	p[i] = p[nt-i];
    }
    
    kiss_fftr(forw, p, (kiss_fft_cpx *) pp);
    
    for (i=0; i < nw; i++) {
	q[o1+i*d1] = crealf(pp[i]);
    }
}

void cosft_inv (float *q, int o1, int d1) {
    int i;

    for (i=0; i < nw; i++) {
	pp[i] = q[o1+i*d1];
    }

/*

    if (0. != dt) {
	for (i=0; i < n; i++) {
	    pp[i] *= cexpf(I*i*dt);
	}
    }
  
*/
  
    kiss_fftri(invs,(const kiss_fft_cpx *) pp, p);
    
    for (i=0; i < nw; i++) {
	q[o1+i*d1] = p[i]/nt;
    }
}

/* 	$Id$	 */
