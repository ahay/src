/* Missing data interpolation with helical filters */

#include <rsf.h>
#include "helix.h"
/*^*/

#include "mis2.h"
#include "helicon.h" 
#include "polydiv.h" 

void mis2(int niter         /* number of iterations */, 
	  int nx            /* model size */, 
	  float *xx         /* model */, 
	  filter aa         /* helix filter */, 
	  const bool *known /* mask for known data */,
	  float eps         /* regularization parameter */,
	  bool doprec       /* to apply preconditioning */) 
/*< interpolate >*/
{
    int ix;
    float *dd;

    if (doprec) {                          /*  preconditioned */
	sf_mask_init(known);
	polydiv_init(nx, aa);
	sf_solver_prec(sf_mask_lop, sf_cgstep, polydiv_lop, nx, nx, nx, 
		       xx, xx, niter, eps, "end");
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

/* 	$Id$	 */

