#include <rsf.h>

#include "parcel.h"
#include "patch.h"

static int np, nw;

/* parcel_init
   -----------
   Initialize parcel parameters.
   dim - dimensionality
   npatch - number of patches
   nwall - data size
   nwind - patch size
*/
void parcel_init(int dim, int *npatch, int *nwall, int *nwind)
{
    int i;

    patch_init (dim, npatch, nwall, nwind);
    np = 1;
    nw = 1;  
    for (i=0; i < dim; i++) {
	np *= npatch[i]; /* compute number of patches */
	nw *= nwind[i];  /* compute window size */
    }
}

void parcel_lop(bool adj, bool add, int n, int mw, 
		float* wall, float* wind)
{
    int ip;

    if (mw != np*nw) sf_error("%s: wrong dimensions",__FILE__);
  
    sf_adjnull (adj, add, n, mw, wall, wind);
    for (ip=0; ip < np; ip++) {
	patch_lop (adj, true, n, nw, wall, wind+ip*nw);
	patch_close ();
    }
}
