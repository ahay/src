#include "hconest.h"
#include "adjnull.h"
#include "helix.h"

static float *x;
static filter aa;

void hconest_init(float *x_in, filter aa_in)
{
    x = x_in;
    aa = aa_in;
}

void hconest_lop(bool adj, bool add, int na, int ny, float *a, float *y)
{
    int  ia, ix, iy;
    
    if (na != aa->nh) sf_error("%s: Wrong data dimensions",__FILE__);

    adjnull(adj, add, na, ny, a, y);

    for (ia = 0; ia < na; ia++) {
	for (iy = aa->lag[ia]; iy < ny; iy++) {  
	    if(aa->mis[iy]) continue;
  
            ix = iy - aa->lag[ia];

	    if( adj) a[ia] -=  y[iy] * x[ix];
	    else     y[iy] -=  a[ia] * x[ix];
	}
    }
}

/* 	$Id: hconest.c,v 1.3 2003/10/01 22:45:56 fomels Exp $	 */
