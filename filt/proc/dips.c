#include <rsf.h>

#include "dips.h"

#include "callp2.h"

static int nd, n1, n2, n12, skip;
static callpass2 *ap;
static float *x, **tmp1, **tmp2;

void dips_init(int nd1, int nw, int nj, int nx, int ny, float** x1)
{
    int id;

    x = x1[0];
    nd = nd1;
    ap = (callpass2 *) sf_alloc(nd,sizeof(callpass2));
    for (id=0; id < nd; id++) {
	ap[id] = callpass2_init(nw, nj, nx, ny);
    }
    tmp1 = sf_floatalloc2(nx,ny);
    tmp2 = sf_floatalloc2(nx,ny);

    n1 = nx;
    n2 = ny;
    n12 = nx*ny;

    skip = nd*nw*nj;
}

void dips_close(void)
{
    int id;

    for (id=0; id < nd; id++) {
	free (ap[id]);
    }
    free (ap);
    free (tmp1[0]);
    free (tmp1);
    free (tmp2[0]);
    free (tmp2);
}

void dips(const float *d, float *b, float **aa)
{
    int id, jd, i1, i2, i;
    float **tmp;

    for (i=0; i < n12; i++) {
	tmp2[0][i] = x[i];
    }

    for (id=0; id < nd; id++) {
	tmp = tmp1; tmp1 = tmp2; tmp2 = tmp;

	callpass21_set (ap[id], d[id]);
	callpass21 (false, ap[id], tmp1, tmp2);
    }

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    i = i2*n1+i1;
	    if (i2 < n2-nd && i1 >= skip && i1 < n1-skip) {
		b[i] = tmp2[i2][i1];
	    } else {
		b[i] = 0;
	    }
	}
    }

    for (id=0; id < nd; id++) {
	for (i=0; i < n12; i++) {
	    tmp2[0][i] = x[i];
	}

	for (jd=0; jd < nd; jd++) {
	    tmp = tmp1; tmp1 = tmp2; tmp2 = tmp;
	    callpass21 (jd == id, ap[jd], tmp1, tmp2);
	}
	
	
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		i = i2*n1+i1;
		if (i2 < n2-nd && i1 >= skip && i1 < n1-skip) {
		    aa[i][id] = tmp2[i2][i1];
		} else {
		    aa[i][id] = 0.;
		}
	    }
	}
    }
}


