#include <rsf.h>

#include "msmis2.h"
#include "mask.h"
#include "mshelicon.h" 

void msmis2(int niter, int nx, int ns, float *xx, msfilter aa, 
	    const bool *known) 
{
    int ix, nxs;
    float *dd;

    nxs = ns*nx;

    dd = sf_floatalloc(nxs);
    for (ix=0; ix < nxs; ix++) {
	dd[ix]=0.;
    }
    
    mshelicon_init(aa);
    sf_solver (mshelicon_lop, sf_cgstep, nx, nxs, xx, dd, niter, 
	       "known", known, "x0", xx, "end");
    
    sf_cgstep_close();
    free(dd);
}

/* 	$Id: msmis2.c,v 1.1 2004/06/11 10:51:34 fomels Exp $	 */

