#include <stdio.h>

#include <rsf.h>

#include "ocmkwallwt.h"
#include "ocpatch.h"
#include "oc.h"

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
void ocmkwallwt(bool inv, int dim, int* npatch, int* nwall, int* nwind, 
		float* windwt, FILE* wallwt)
{
    int j, iw, ip, np, nw;
    size_t n;
    float *tmp;
    
    np = 1; 
    n = sizeof(float);
    nw = 1;

    for (j=0; j < dim; j++) {
	np *= npatch[j];
	n *= nwall[j];
	nw *= nwind[j];
    }

    tmp = sf_floatalloc(nw);

    oc_zero (n, wallwt);
    ocpatch_init(dim, nw, np, npatch, nwall, nwind);
  
    for (ip=0; ip < np; ip++) {
	ocpatch_lop (ip, false, wallwt, tmp);
	for (iw=0; iw < nw; iw++) {
	    tmp[iw] += windwt[iw];
	}
	ocpatch_lop(ip, true, wallwt, tmp);
    }

    if (inv) oc_invert(n,wallwt);
}


