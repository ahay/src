#include "pef.h"
#include "hconest.h"
#include "cgstep.h"
#include "bigsolver.h"

void find_pef(int nd, float* dd, filter aa, int niter) 
{
    int i;

    for (i=0; i < nd; i++) dd[i] = -dd[i];

    hconest_init( dd, aa);
    solver(hconest_lop, cgstep, aa->nh, nd, aa->flt, dd, niter, 
	   "x0", aa->flt, "end");
    cgstep_close();
    
    for (i=0; i < nd; i++) dd[i] = -dd[i];	
}
