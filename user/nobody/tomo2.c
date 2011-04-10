#include <math.h>

#include <rsf.h>

#include "tomo2.h"

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif

static int nd, nf, m1, m2, ***nxy, nr, *rl;
static bool **mask;
static float ***w1, ***w2;

void tomo2_init (float*** rays, int *raylen, int nrays, 
		 float o1, float o2, float d1, float d2,
		 int n1, int n2, 
		 interpolator interp, int nf_in)
{
    int ir, id, i1, i2, nd; 
    float x1, x2, rx;

    nf = nf_in;
    m1 = n1;
    m2 = n2;
    nr = nrays;
    rl = raylen;

    nxy = (int***) sf_alloc(nr,sizeof(int**));
    mask = (bool**) sf_alloc(nr,sizeof(bool*));
    w1 = (float***) sf_alloc(nr,sizeof(float**));
    w2 = (float***) sf_alloc(nr,sizeof(float**));

    for (ir = 0; ir < nr; ir++) {
	nd = rl[ir];
	nxy[ir] = sf_intalloc2(2,nd);
	mask[ir] = sf_boolalloc(nd);
	w1[ir] = sf_floatalloc2(nf,nd);
	w2[ir] = sf_floatalloc2(nf,nd);
	for (id = 0; id < nd; id++) {
	    rx = (rays[ir][id][0] - o1)/d1;
	    i1 = (int) floor(rx + 1. - 0.5*nf);
	    x1 = rx - floor(rx);
	
	    rx = (rays[ir][id][1] - o2)/d2;
	    i2 = (int) floor(rx + 1. - 0.5*nf);
	    x2 = rx - floor(rx);
   
	    if (i1 > - nf && i1 < n1 &&
		i2 > - nf && i2 < n2) {
		mask[ir][id] = false; 
		interp (x1, nf, w1[ir][id]);
		interp (x2, nf, w2[ir][id]);
		nxy[ir][id][0] = i1;
		nxy[ir][id][1] = i2;
	    } else {
		mask[ir][id] = true;
	    }
	}
    }
}

void  tomo2_lop (bool adj, bool add, int nm, int ny, float* x, float* ord)
{ 
    int ir, id, i0, j0, i, j, im;
    float w;

    if (ny != nr) sf_error("%s: wrong dimensions: %d != %d",__FILE__,ny,nr);

    sf_adjnull (adj,add,nm,nd,x,ord);

    for (ir = 0; ir < nr; ir++) {
	nd = rl[ir];
	for (id=0; id < nd; id++) {
	    if (mask[ir][id]) continue;
	    i0 = nxy[ir][id][0]; 
	    j0 = nxy[ir][id][1]; 
	    for (j = MAX(0,-j0); j < MIN(nf,m2-j0); j++) {
		w = w2[ir][id][j];
		for (i = MAX(0,-i0); i < MIN(nf,m1-i0); i++) { 
		    im = (i+i0) + (j+j0)*m1;
		    if( adj) { 
			x[im] += ord[ir] * w * w1[ir][id][i];
		    } else {
			ord[ir] += x[im] * w * w1[ir][id][i];
		    }
		}
	    }
	}
    }
}

void tomo2_close (void)
{
    int ir;

    for (ir = 0; ir < nr; ir++) {
	free (**nxy);
	free (*nxy);
	free (*mask);
	free (**w1);
	free (**w2);
	free (*w1);
	free (*w2);
    }
    free (nxy);
    free (mask);
    free (w1);
    free (w2);
}

/* 	$Id$	 */
