#include <rsf.h>

#include "laplac2.h"

static int n1, n2;

void laplac2_init(int m1, int m2)
{
    n1 = m1;
    n2 = m2;
}

void laplac2_lop(bool adj, bool add, int np, int nr, float *p, float *r)
{
    int i1,i2,j;

    sf_adjnull(adj,add,np,nr,p,r);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    j = i1+i2*n1;
	    if (i1 > 0) {
		if (adj) {
		    p[j-1] -= r[j];
		    p[j]   += r[j];
		} else {
		    r[j] += p[j] - p[j-1];
		}
	    }
	    if (i1 < n1-1) {
		if (adj) {
		    p[j+1] -= r[j];
		    p[j]   += r[j];
		} else {
		    r[j] += p[j] - p[j+1];
		}
	    }
	    if (i2 > 0) {
		if (adj) {
		    p[j-n1] -= r[j];
		    p[j]    += r[j];
		} else {
		    r[j] += p[j] - p[j-n1];
		}
	    }
	    if (i2 < n2-1) {
		if (adj) {
		    p[j+n1] -= r[j];
		    p[j]    += r[j];
		} else {
		    r[j] += p[j] - p[j+n1];
		}
	    }
	}
    }
}

/* 	$Id: laplac2.c,v 1.1 2004/02/14 06:52:41 fomels Exp $	 */
