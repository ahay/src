#include "pef.h"
#include "hconest.h"

void find_pef(int nd, float* dd, filter aa, int niter) 
{
    hconest_init( dd, aa);
    sf_solver(hconest_lop, sf_cgstep, aa->nh, nd, aa->flt, dd, niter, 
	   "x0", aa->flt, "end");
    sf_cgstep_close();
}

/* 	$Id: pef.c,v 1.4 2003/10/21 15:09:08 fomels Exp $	 */

