#include <math.h>

#include <rsf.h>

#include "compress.h"
#include "helix.h"

filter compress( filter aa, float eps) {
    filter bb;
    bool* keep;
    int i, k;

    keep = sf_boolalloc (aa->nh);
    k = 0;
    for(i=0; i < aa->nh; i++) {
	if ( fabs( aa->flt[i]) > eps) {
	    keep[i] = true;
	    k++;
	} else {
	    keep[i] = false;
	}
    }
    bb = allocatehelix( k);
    k = 0;
    for(i=0; i < aa->nh; i++) {
	if (keep[i]) {
	    bb->flt[k] = aa->flt[i];
	    bb->lag[k] = aa->lag[i];
	    k++;
	}
    }
    deallocatehelix( aa);
    free (keep);

    return bb;
}

/* 	$Id: compress.c,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */
