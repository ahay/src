#include "mspef.h"
#include "mshconest.h"

void find_pef(int nd, float* dd, msfilter aa, int niter) 
{
    float *ee;
    int is, id, ns;

    ns = aa->ns;
    ee = sf_floatalloc(nd*ns);
    for (is=0; is < ns; is++) {
	for (id=0; id < nd; id++) {
	    ee[id+is*nd] = dd[id];
	}
    }

    mshconest_init( dd, aa);
    sf_solver(mshconest_lop, sf_cgstep, aa->nh, nd*ns, aa->flt, ee, niter, 
	      "x0", aa->flt, "end");
    sf_cgstep_close();

    free(ee);
}

/* 	$Id: mspef.c,v 1.1 2004/06/11 10:51:34 fomels Exp $	 */

