#include "nhelicon.h"
#include "nhelix.h"
#include "copy.h"

static nfilter aa;

void nhelicon_init(nfilter aa_in)
{
    aa = aa_in;
}

void nhelicon_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy)
{
    int iy, ix, ia, ip, na, *lag;
    float *flt;

    copy_lop(adj,add,nx,ny,xx,yy);

    for (iy=0; iy < ny; iy++) {    
        if (NULL != aa->mis && aa->mis[iy]) continue;
        ip = aa->pch[iy]; 
	lag = aa->hlx[ip]->lag; 
	flt = aa->hlx[ip]->flt;
	na = aa->hlx[ip]->nh;
        for (ia=0; ia < na; ia++) { 
            ix = iy - lag[ia]; 
	    if(ix < 0) continue;
            if (adj) {
		xx[ix] += yy[iy] * flt[ia];
	    } else {
		yy[iy] += xx[ix] * flt[ia];
            }
	}
    }
}
