#include <rsf.h>

#include "hypotenusei.h"

static int nt, nx;
static float x0,dx, t0,dt;

void imospray_init (float slow, float y0, float dy, float z0, float dz, 
		    int nz, int ny)
{
    x0 = y0*slow;
    dx = dy*slow;
    t0 = z0;
    dt = dz;
    nt = nz;
    nx = ny;

    hypotenusei_init (nt);
}

void imospray_lop(bool adj, bool add, int n1, int n2, 
		   float *stack, float *gather)
{
    int ix;
    float x;

    if (n1 != nt || n2 != nt*nx) sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull(adj,add,n1,n2,stack,gather);

    for (ix=0; ix < nx; ix++) {
	x = x0 + dx*ix;

        hypotenusei_set (t0, dt, x);
	hypotenusei_lop (adj, true, nt, nt, stack, gather+ix*nt);
    }
}

void imospray_close(void)
{
    hypotenusei_close ();
}
