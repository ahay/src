#include <stdio.h>

#include <rsf.h>

#include "ocparcel.h"
#include "ocpatch.h"

static int np, nw;
static float *tmp;

void ocparcel_init(int dim, int *npatch, int *nwall, int *nwind)
{
    int i;

    np = 1;
    nw = 1;  
    for (i=0; i < dim; i++) {
	np *= npatch[i]; /* compute number of patches */
	nw *= nwind[i];  /* compute window size */
    }
    
    ocpatch_init (dim, nw, np, npatch, nwall, nwind);
    tmp = sf_floatalloc(nw);
}

void ocparcel_close(void)
{
    free (tmp);
}

void ocparcel_lop(bool adj, int n, int mw, FILE* wall, float* wind)
{
    int ip, iw;

    if (mw != np*nw) sf_error("%s: wrong dimensions",__FILE__);
  
    for (ip=0; ip < np; ip++) {
	if (adj) {
	    ocpatch_lop (ip, false, wall, tmp);
	    for (iw=0; iw < nw; iw++) {
		tmp[iw] += wind[iw+ip*nw];
	    }
	    ocpatch_lop (ip, true, wall, tmp);
	} else {
	    ocpatch_lop (ip, false, wall, wind+ip*nw);
	}
    }
}
