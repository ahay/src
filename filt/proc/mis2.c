#include <rsf.h>

#include "mis2.h"
#include "mask.h"
#include "helicon.h" 
#include "polydiv.h" 
#include "cgstep.h"
#include "bigsolver.h"

void mis2(int niter, int nx, float *xx, filter aa, 
	  const bool *known, bool doprec) 
{
    int ix;
    float *dd;

    if (doprec) {                          /*  preconditioned */
	mask_init(known);
	polydiv_init(nx, aa);
	solver_prec(mask_lop, cgstep, polydiv_lop, nx, nx, nx, xx, xx, 
		    niter, 0., "end");
	polydiv_close();
    } else {                               /*  regularized */
	dd = sf_floatalloc(nx);
	for (ix=0; ix < nx; ix++) {
	    dd[ix]=0.;
	}

	helicon_init(aa);
	solver (helicon_lop, cgstep, nx, nx, xx, dd, niter, 
		"known", known, "end");
	free(dd);
    }
    cgstep_close();
}



