#include <rsf.h>

#include "boxfilter.h"
#include "helix.h"

void box (int dim, const int *nd, const int *center, const int *na, 
	  const filter aa, int nc, float* cube) 
{
    int ii[SF_MAX_DIM];
    int j, lag0a, lag0d, id, ia;

    for (ia=0; ia < nc; ia++) {
	cube[ia] = 0.;
    }
    lag0a = sf_cart2line(dim, na, center);  /* locate the 1.0 in na_cube. */
    cube[lag0a] = 1.;                       /* place it. */
    lag0d = sf_cart2line(dim, nd, center);  /* locate the 1.0 in nd_cube. */
    for (j=0; j < aa->nh; j++) { /* inspect the entire helix */
	id = aa->lag[j] + lag0d;
	sf_line2cart(dim, nd, id, ii);	/* ii = cartesian indices  */
	ia = sf_cart2line(dim, na, ii);	/* ia = linear index in aa */
	cube[ia] = aa->flt[j];		/* copy the filter coefficient */
    }
}

/* 	$Id: boxfilter.c,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */

