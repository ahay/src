#include <rsf.h>

#include "mis1.h"

#include "tcai1.h"

static int nx, ny;
static float *zero;

void mis1_init(int n1, int na, float *aa)
{
    int iy;

    nx=n1;
    ny=nx+na-1;

    zero = sf_floatalloc(ny);
    for (iy=0; iy < ny; iy++) {
	zero[iy]=0.;
    }

    tcai1_init(na,aa);
}

void mis1_close(void)
{
    free(zero);
}

void mis1(int niter, float *xx, const bool *known) 
{
    sf_solver (tcai1_lop, sf_cgstep, nx, ny, xx, zero, niter, 
	       "x0", xx, "known", known, "end");
    sf_cgstep_close();
}

/* 	$Id: mis1.c,v 1.1 2004/04/01 15:41:54 fomels Exp $	 */

