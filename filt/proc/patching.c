#include <rsf.h>

#include "patching.h"
#include "patch.h"
#include "mkwallwt.h"

/* patching
   --------
   Apply a linear operator in patches.
   modl -input. 
   data - output,
   dim - dimensionality,
   npatch - number of patches
   nwall - data size
   nwind - window size,
   windwt - window weight */
void patching(sf_operator oper, float* modl, float* data, 
	      int dim, int* npatch, int* nwall, int* nwind, float* windwt)
{
    float *winmodl, *windata, *wallwt;
    int i, j, iw, ip, np, n, nw;

    np = 1; 
    n = 1;
    nw = 1;

    for (j=0; j < dim; j++) {
	np *= npatch[j];
	n *= nwall[j];
	nw *= nwind[j];
    }
  
    winmodl = sf_floatalloc(nw);
    windata = sf_floatalloc(nw);
    wallwt  = sf_floatalloc(n);

    for (i=0; i < n; i++) {
	data[i] = 0.;
    }

    patch_init(dim, npatch, nwall, nwind);
    for (ip = 0; ip < np; ip++) {
	/* modl -> winmodl */
	patch_lop(false, false, n, nw, modl, winmodl);
	/* winmodl -> windata */
	oper(false, false, nw, nw, winmodl, windata);
	/* apply window weighting */
	for (iw=0; iw < nw; iw++) {
	    windata[iw] *= windwt[iw];
	}
	/* data <- windata */
	patch_lop(true, true, n, nw, data, windata);
	patch_close();
    }

    /* windwt -> wallwt */
    mkwallwt(dim, npatch, nwall, nwind, windwt, wallwt);

    /* apply wall weighting */
    for (i=0; i < n; i++) {
	data[i] *= wallwt[i];
    }

    free (winmodl);
    free (windata);
    free (wallwt);
}

