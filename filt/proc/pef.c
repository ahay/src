#include "pef.h"
#include "hconest.h"
#include "cgstep.h"
#include "bigsolver.h"

void find_pef(int nd, float* dd, filter aa, int niter) 
{
    hconest_init( dd, aa);
    solver(hconest_lop, cgstep, aa->nh, nd, aa->flt, dd, niter, 
	   "x0", aa->flt, "end");
    cgstep_close();
}

/* 	$Id: pef.c,v 1.3 2003/10/01 22:45:56 fomels Exp $	 */

