#include <rsf.h>

#include "npolydiv.h"
#include "nhelix.h"

static int nd;
static nfilter aa;
static float *tt;

void npolydiv_init (int nd_in, nfilter aa_in)
{
    nd = nd_in;
    aa = aa_in;
    tt = sf_floatalloc(nd);
}

void npolydiv_close(void)
{
    free (tt);
}


void npolydiv_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy)
{
    int id, ia, na, ix, iy, ip, *lag;
    float *flt;

    sf_adjnull(adj,add,nx,ny,xx,yy);

    for (id=0; id < nd; id++) {
	tt[id] = adj? yy[id]: xx[id];
    }

    if (adj) {
        for (iy=nd-1; iy >= 0; iy--) { 
	    ip = aa->pch[iy];
	    lag = aa->hlx[ip]->lag; 
	    flt = aa->hlx[ip]->flt;
	    na = aa->hlx[ip]->nh;
	    for (ia=0; ia < na; ia++) {
		ix = iy - lag[ia];     
		if (ix < 0)  continue;
		tt[ix] -=  flt[ia] * tt[iy];
	    } 
	}
	for (id=0; id < nd; id++) {
	    xx[id] += tt[id];
	}
    } else { 
        for (iy=0; iy < nd; iy++) { 
	    ip = aa->pch[iy];
	    lag = aa->hlx[ip]->lag; 
	    flt = aa->hlx[ip]->flt;
	    na = aa->hlx[ip]->nh;
	    for (ia=0; ia < na; ia++) {
		ix = iy - lag[ia]; 
		if (ix < 0)  continue;
		tt[iy] -=  flt[ia] * tt[ix];
	    } 
	}
	for (id=0; id < nd; id++) {
	    yy[id] += tt[id];
        }
    }
}

