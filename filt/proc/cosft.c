#include <rsf.h>

#include "cosft.h"

static int nt, nw, n1;
static float *p, dt;
static complex float *pp;

void cosft_init(int n1, float o1, float d1) {
    nt = sf_npfar(2*(n1-1));
    nw = nt/2+1;
    p  = sf_floatalloc (nt);
    pp = sf_complexalloc(nw);
    n1 = n1;
    dt = 2.*SF_PI*o1/(nt*d1);
}

void cosft_close(void) {
    free (p);
    free (pp);
}

void cosft_frw (float *q, int o1, int d1) {
    int i1;

    for (i1=0; i1 < n1; i1++) {
	p[i1] = q[o1+i1*d1];
    }
    for (i1=n1; i1 <= nt/2; i1++) {
	p[i1]=0.0;
    }
    for (i1=nt/2+1; i1 < nt; i1++) {
	p[i1] = p[nt-i1];
    }
    
    sf_pfarc(1,nt,p,pp);
	    
    if (0. != dt) {
	for (i1=0; i1 < nw; i1++) {
	    pp[i1] *= cexpf(I*i1*dt);
	}
    }

    for (i1=0; i1 < n1; i1++) {
	q[o1+i1*d1] = crealf(pp[n1]);
    }
}

void cosft_inv (float *q, int o1, int d1) {
    int i1;

    for (i1=0; i1 < n1; i1++) {
	pp[i1] = q[o1+i1*d1];
    }
    for (i1=n1; i1 < nw; i1++) {
	pp[i1] = 0.;
    }
    
    sf_pfacr(-1,nt,pp,p);

    for (i1=0; i1 < n1; i1++) {
	q[o1+i1*d1] = p[i1]/nt;
    }
}
