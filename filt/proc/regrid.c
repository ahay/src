#include <rsf.h>

#include "regrid.h"

/* Regrid
   ------
   Convert a helix filter from one data size to another. */
void regrid( int dim, const int* nold, const int* nnew, filter aa) {
    int i, h0, h1, h, ii[SF_MAX_DIM];

    for (i=0; i < dim; i++) {
	ii[i] = nold[i]/2-1;
    }
  
    h0 = sf_cart2line( dim, nold, ii); /* lag of near middle point on nold */
    h1 = sf_cart2line( dim, nnew, ii); /* lag                      on nnew */
    for (i=0; i < aa->nh; i++) { /* forall given filter coefficients */
	h = aa->lag[i] + h0;
	sf_line2cart( dim, nold, h, ii);
	aa->lag[i] = sf_cart2line( dim, nnew, ii) - h1;
    }
}

/* 	$Id: regrid.c,v 1.1 2004/06/11 10:47:16 fomels Exp $	 */

