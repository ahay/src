#include <stdlib.h>

#ifdef __GNUC__
#define alloca __builtin_alloca
#else /* not GNU C  */
#if (!defined (__STDC__) && defined (sparc)) || defined (__sparc__) || defined (__sparc) || defined (__sgi)
#include <alloca.h>
#endif
#endif

#include <rsf.h>

#include "createhelix.h"
#include "helix.h"

filter createhelix(int ndim, int* nd, int* center, int* gap, int* na)
{
    filter aa;
    int ii[SF_MAX_DIM], na123, ia, nh, lag0a,lag0d, *lag, i;
    bool skip;

    for (na123 = 1, i=0; i < ndim; i++) {
	na123 *= na[i];
    }

    lag = (int*) alloca(na123*sizeof(int));

    lag0a = sf_cart2line (ndim, na, center); /* index pointing to the "1.0" */

    nh= 0;   
    for (ia = 1+lag0a; ia < na123; ia++) { /* ia  is fortran linear index. */
	sf_line2cart(ndim, na, ia, ii);
	
	skip = false;
	for (i=0; i < ndim; i++) {
	    if (ii[i] < gap[i]) {
		skip = true;
		break;
	    }
	}
	if (skip) continue;
	
	lag[nh] = sf_cart2line(ndim, nd, ii);
	nh++;                        /* got another live one */
    }
    lag0d = sf_cart2line(ndim,  nd, center); /* center shift for nd cube */

    aa = allocatehelix(nh);		 /* nh becomes size of filter */

    for (ia=0; ia < nh; ia++) {
	aa->lag[ia] = lag[ia] - lag0d; 
	aa->flt[ia] = 0.;
    }

    return aa;
}

/* 	$Id: createhelix.c,v 1.5 2003/11/17 19:42:01 fomels Exp $	 */
