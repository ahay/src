#include <rsf.h>

#include "cadj.h"

void adjnull(bool adj, bool add, int nx, int ny, 
	     float complex* x, float complex* y)
{
    int i;
    
    if(add) return;
    if(adj) {
	for (i=0; i < nx; i++) {
	    x[i] = 0.;
	}
    } else {
	for (i=0; i < ny; i++) {
	    y[i] = 0.;
	}
    }
}

/* 	$Id: cadj.c,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */

