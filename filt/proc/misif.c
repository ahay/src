#include <rsf.h>

#include "misif.h"

#include "peftc.h"
#include "tcai2.h"

void misif1 (int niter, int na, int nx, float *xx, float *aa, bool *mm) {
    float *bb, *dat, *x;
    int id, nd;

    nd = nx+na;

    dat = sf_floatalloc(nd);
    x =  sf_floatalloc(nd);
    bb = x+nx;

    for (id=0; id < nd; id++) {
	dat[id] = 0.;
    }

    for (id=0; id < nx; id++) {
	x[id] = mm[id]? xx[id]: 0.;
    }
    for (id=nx; id < nd; id++) {
	x[id] = mm[id]? 1.: 0.;
    }

    peftc_init (na, nx, bb, x);
    tcai2_init (na, nx, bb);

    sf_solver (peftc_lop, sf_cgstep, nd, nd, x, dat, niter, "x0", x, 
               "nloper", tcai2_lop, "known", mm, "end");
    sf_cgstep_close ();

    for (id=0; id < nx; id++) {
	xx[id] = x[id];
    }
    for (id=nx; id < nd; id++) {
	aa[id-nx] = x[id];
    }

    free(x);
    free(dat);
}
