#include <rsf.h>

#include "patch.h"

static int dim, ipatch;
static int *npatch, *nwind, *nwall, ii[SF_MAX_DIM], jj[SF_MAX_DIM]; 

/* patch_init
   ----------
   Initialize patch parameters.
   dim - dimensionality
   npatch - number of patches
   nwall - data size
   nwind - patch size
*/
void patch_init(int dim_in, int* npatch_in, int* nwall_in, int* nwind_in)
{
    dim = dim_in;
    npatch = npatch_in; 
    nwall = nwall_in; 
    nwind = nwind_in; 
    ipatch = 0;
}

/* patch_lop
   ---------
   Patch operation. 
   if (!adj) extract a patch (wind from wall)
   if (adj)  insert  a patch (wind to wall)
*/
void patch_lop (bool adj, bool add, 
		int nx, int ny, float* wall, float* wind)
{
    int i, j, shift;
 
    sf_adjnull (adj, add, nx, ny, wall, wind);
    sf_line2cart(dim, npatch, ipatch, jj);  
    for(i = 0; i < dim; i++) {
	if(npatch[i] == 1) {
	    jj[i] = 0;
	} else if (jj[i] == npatch[i]-1) {
	    jj[i] = nwall[i] - nwind[i];
	} else {	    
	    jj[i] = jj[i]*(nwall[i] - nwind[i])/(npatch[i] - 1.0);
	}
    }

    /* shift till the patch start */
    shift = sf_cart2line(dim, nwall, jj); 
    for(i = 0; i < ny; i++) {
	sf_line2cart(dim, nwind, i, ii); 
	j = sf_cart2line(dim, nwall, ii) + shift;   
	if (adj) wall[j] += wind[i];
	else     wind[i] += wall[j];
    }
}

/* patch_close
   -----------
   Move to the next patch.
*/
void patch_close(void)
{
    ipatch++; /* increase patch counter */
}
