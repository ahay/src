#include "hconest.h"
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
    for (ia = 0; ia < na; ia++) {
	for (iy = aa->lag[ia]; iy < ny; iy++) {  
	    if(aa->mis[iy]) continue;
  
            ix = iy - aa->lag[ia];

	    if( adj) a[ia] +=  y[iy] * x[ix];
	    else     y[iy] +=  a[ia] * x[ix];
	}
    }
}
