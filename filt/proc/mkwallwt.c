#include <rsf.h>

#include "mkwallwt.h"
#include "patch.h"

/* mkwallwt
   --------
   Make wall weighting from window weighting.
   dim - dimensionality
   npatch - number of patches
   nwall - data size
   nwind - patch size
   windwt[product(nwind)] - window weighting (input)
   wallwt[product(nwall)] - wall weighting (output)
*/
void mkwallwt(int dim, int* npatch, int* nwall, int* nwind, 
	      float* windwt, float* wallwt)
{
    int i, j, ip, np, n, nw;

    np = 1; 
    n = 1;
    nw = 1;

    for (j=0; j < dim; j++) {
	np *= npatch[j];
	n *= nwall[j];
	nw *= nwind[j];
    }

    for (i = 0; i < n; i++) {
	wallwt[i] = 0.;
    }

    patch_init(dim, npatch, nwall, nwind);
  
    for (ip=0; ip < np; ip++) {
	patch_lop(true, true, n, nw, wallwt, windwt);
	patch_close ();
    }

    for (i = 0; i < n; i++) {
	if ( wallwt[i] != 0.) wallwt[i] = 1. / wallwt[i];
    }
}


