#include <rsf.h>

#include "igrad1.h"

void  igrad1_lop(bool adj, bool add, int nx, int ny, float *xx, float *yy)
{
    int i;

    sf_adjnull(adj,add,nx,ny,xx,yy);

    for (i=0; i < nx-1; i++) {
        if (adj) {
	    xx[i+1] += yy[i];
	    xx[i]   -= yy[i];
	} else {
	    yy[i] += xx[i+1] - xx[i];
	}
    }
}

/* 	$Id: igrad1.c,v 1.1 2004/02/14 06:52:41 fomels Exp $	 */


