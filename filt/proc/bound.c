#include <stdlib.h>
#include <stdio.h>

#include "bound.h"
#include "helicon.h"
#include "regrid.h"

/* mark helix filter outputs where input is off data. */
void bound (int dim, const int *nold, const int *nd, const int *na, 
	    const filter aa) 
{
    int nb[SF_MAX_DIM], ii[SF_MAX_DIM];
    float *xx, *yy;
    int iy, my, ib, mb, i;
    
    my = mb = 1;
    for (i=0; i < dim; i++) {
	nb[i] = nd[i] + 2*na[i]; /* nb is a bigger space to pad into. */   
	mb *= nb[i];
	my *= nd[i];
    }

    /* two large spaces, equal size */
    xx = sf_floatalloc(mb);
    yy = sf_floatalloc(mb);

    for (ib=0; ib < mb; ib++) {
	xx[ib] = 0.;  /* zeros */
	/* surround the zeros with many ones */
	sf_line2cart(dim, nb, ib, ii);           /* ii( ib) */
	for (i=0; i < dim; i++) {
	    if(ii[i]+1 <= na[i] ||  ii[i]+1 > nb[i]-na[i]) {
		xx[ib] = 1.; 
		break;
	    }
	}	
    }
    helicon_init( aa);		     /* give aa to helicon.lop */
    regrid(dim, nold, nb, aa);  
    for (i=0; i < aa->nh; i++) {
	aa->flt[i] = 1.;		/* put all 1's in filter */
    }
    helicon_lop(false, false, mb, mb, xx, yy);	/* apply filter */
    regrid(dim, nb, nd, aa);  
    for (i=0; i < aa->nh; i++) {
	aa->flt[i] = 0.;		/* remake filter for orig data */
    }

    aa->mis = sf_boolalloc(my); /* attach missing designation to y_filter */

    for (iy = 0; iy < my; iy++) {  /* map from unpadded to padded space */
	sf_line2cart(dim, nd, iy, ii);
	for (i=0; i < dim; i++) {
	    ii[i] += na[i];
	}
	ib = sf_cart2line(dim, nb, ii);
	
	aa->mis[iy] = ( yy[ib] > 0.);  /* true where inputs missing */
    }
    
    free (xx);
    free (yy);
} 

/* 	$Id: bound.c,v 1.1 2004/06/11 10:47:16 fomels Exp $	 */
