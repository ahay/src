#include <rsf.h>

#include "causint.h"

/* Causal integration */
void causint_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy)
{
    int i;       
    float t;

    sf_adjnull (adj, add, nx, ny, xx, yy);

    t = 0.;
    if ( adj) {
	for (i=nx-1; i >= 0; i--) {
	    t += yy[i];
	    xx[i] += t;
	}
    } else {
	for (i=0; i < nx; i++) {
	    t += xx[i];
	    yy[i] += t;
	}
    }
}

/* 	$Id: causint.c,v 1.3 2003/10/21 15:09:08 fomels Exp $	 */

