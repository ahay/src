#include <rsf.h>

#include "deburst.h"

#include "copy.h"
#include "icai1.h"

void deburst (int n, int niter, sf_weight wght, float eps, 
	      const float *dd, float *hh) 
{
    float aa[] = {-1.,2.,-1.}; /* laplacian filter */

    icai1_init (3, aa, 1); 
    sf_solver_reg(copy_lop, sf_cgstep, icai1_lop, n, n, n, hh, dd, niter, eps,
		  "wght", wght, "nfreq", 1, "nmem", 0, "end");
    sf_cgstep_close();
}
