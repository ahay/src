#include <rsf.h>

#include "mis2.h"
#include "mask.h"
#include "helicon.h" 
#include "polydiv.h" 

void mis2(int niter, int nx, float *xx, filter aa, 
	  const bool *known, bool doprec) 
{
    int ix;
    float *dd;

    if (doprec) {                          /*  preconditioned */
	mask_init(known);
	polydiv_init(nx, aa);
	sf_solver_prec(mask_lop, sf_cgstep, polydiv_lop, nx, nx, nx, xx, xx, 
		       niter, 0., "end");
	polydiv_close();
    } else {                               /*  regularized */
	dd = sf_floatalloc(nx);
	for (ix=0; ix < nx; ix++) {
	    dd[ix]=0.;
	}

	helicon_init(aa);
	sf_solver (helicon_lop, sf_cgstep, nx, nx, xx, dd, niter, 
		   "known", known, "x0", xx, "end");
	free(dd);
    }
    sf_cgstep_close();
}

/* 	$Id: mis2.c,v 1.4 2004/04/06 02:03:03 fomels Exp $	 */

