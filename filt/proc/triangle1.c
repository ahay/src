#include "adjnull.h"
#include "triangle.h"
#include "triangle1.h"

#include <rsf.h>

static int nd;
static triangle tr;
static float *tmp;

void triangle1_init (int nbox, int ndat)
{
    nd = ndat;
    tr = triangle_init (nbox,ndat);
    tmp = sf_floatalloc (ndat);
}

void triangle1_lop (bool adj, bool add, int nx, int ny, float* x, float* y)
{
    int i;

    if (nx != ny || nx != nd) sf_error("%s: Wrong data dimensions",__FILE__);

    adjnull (adj,add,nx,ny,x,y);
    if (adj) {
	for (i=0; i < nd; i++) {
	    tmp[i] = y[i];
	}
	smooth (tr, 0, 1, false, tmp);
	for (i=0; i < nd; i++) {
	    x[i] += tmp[i];
	}
    } else {
	for (i=0; i < nd; i++) {
	    tmp[i] = x[i];
	}
	smooth (tr, 0, 1, false, tmp);
	for (i=0; i < nd; i++) {
	    y[i] += tmp[i];
	}
    }
}

void triangle1_close(void)
{
    free (tmp);
    triangle_close (tr);
}
