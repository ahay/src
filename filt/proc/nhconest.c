#include <rsf.h>

#include "nhconest.h"
#include "nhelix.h"

static float *x;
static nfilter aa;
static int nhmax;

void nhconest_init(float *x_in, nfilter aa_in, int nhmax_in)
{
    x = x_in;
    aa = aa_in;
    nhmax = nhmax_in;
}

void nhconest_lop(bool adj, bool add, int naa, int ny, float *a, float *y)
{
    int ia, na, ix, iy, ip, *lag;

    if (naa != nhmax*aa->np) sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull(adj,add,naa,ny,a,y);

    for (iy=0; iy < ny; iy++) {
	if (aa->mis[iy]) continue;

        ip = aa->pch[iy]; 
	lag = aa->hlx[ip]->lag;
	na = aa->hlx[ip]->nh;

        for (ia=0; ia < na; ia++) {
	    ix = iy - lag[ia]; 
	    if (ix < 0) continue;
  
	    if (adj) {
		a[ia+nhmax*ip] += y[iy] * x[ix];
	    } else {
		y[iy] += a[ia+nhmax*ip] * x[ix];
	    }
        }
    }
}
