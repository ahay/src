#include <math.h>

#include <rsf.h>

#include "int1.h"

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif

static int nd, nf, m1, *nx;
static bool *mask, allocated=false;
static float **w1;

void  int1_init (float* coord, float o1, float d1, int n1, 
		 interpolator interp, int nf_in, int nd_in)
{
    int id, i1; 
    float rx;

    nf = nf_in;
    nd = nd_in;
    m1 = n1;

    if (!allocated) {
	nx = sf_intalloc(nd);
	mask = sf_boolalloc(nd);
	w1 = sf_floatalloc2(nf,nd);
	allocated = true;
    }

    for (id = 0; id < nd; id++) {
	rx = (coord[id] - o1)/d1;
	i1 = (int) floorf(rx + 1. - 0.5*nf);
	rx -= floorf(rx);
   
	if (i1 > - nf && i1 < n1) {
	    mask[id] = false; 
	    interp (rx, nf, w1[id]);
	    nx[id] = i1;
	} else {
	    mask[id] = true;
	}
    }
}

void  int1_lop (bool adj, bool add, int nm, int ny, float* x, float* ord)
{ 
    int id, i0, i, im;
    
    if (ny != nd) sf_error("%s: wrong data size: %d != %d",__FILE__,ny,nd);

    sf_adjnull (adj,add,nm,nd,x,ord);

    for (id=0; id < nd; id++) {
	if (mask[id]) continue;
	i0 = nx[id];

	for (i = MAX(0,-i0); i < MIN(nf,m1-i0); i++) { 
	    im = i+i0;
	    if( adj) { 
		x[im] += ord[id] * w1[id][i];
	    } else {
		ord[id] += x[im] * w1[id][i];
	    }
	}
    }
}

void int1_close (void)
{
    if (allocated) {
	free (nx);
	free (mask);
	free (*w1);
	free (w1);
    }
}

/* 	$Id: int1.c,v 1.6 2004/04/12 15:40:43 fomels Exp $	 */
