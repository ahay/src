#include <math.h>

#include <rsf.h>

#include "int1.h"
#include "adjnull.h"

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif

static int nd, nf, m1, *nx;
static bool *mask;
static float **w1;

void  int1_init (float* coord, float o1, float d1, int n1, 
		 interpolator interp, int nf_in, int nd_in)
{
    int id, i1; 
    float x1, rx;

    nf = nf_in;
    nd = nd_in;
    m1 = n1;

    nx = sf_intalloc(nd);
    mask = sf_boolalloc(nd);
    w1 = sf_floatalloc2(nf,nd);

    for (id = 0; id < nd; id++) {
	rx = (coord[id] - o1)/d1;
	i1 = (int) floor(rx + 1. - 0.5*nf);
	x1 = rx - floor(rx);
   
	if (i1 > - nf && i1 < n1) {
	    mask[id] = false; 
	    interp (x1, nf, w1[id]);
	    nx[id] = i1;
	} else {
	    mask[id] = true;
	}
    }
}

void  int1_lop (bool adj, bool add, int nm, int nd, float* x, float* ord)
{ 
    int id, i0, i, im;
    
    adjnull (adj,add,nm,nd,x,ord);

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
    free (nx);
    free (mask);
    free (*w1);
    free (w1);
}
