#include <stdlib.h>

#ifdef __GNUC__
#ifndef alloca
#define alloca __builtin_alloca
#endif
#else /* not GNU C  */
#if (!defined (__STDC__) && defined (sparc)) || defined (__sparc__) || defined (__sparc) || defined (__sgi)
#include <alloca.h>
#endif
#endif

#include <rsf.h>

#include "createmshelix.h"
#include "createhelix.h"
#include "bound.h"

msfilter createmshelix(int ndim, int* nd, int* center, int* gap, 
		       int ns, int *jump, int* na)
{
    msfilter msaa;
    filter aa;
    int is, ih, nh, id, n123, nb[SF_MAX_DIM];

    aa = createhelix(ndim, nd, center, gap, na);
    nh = aa->nh;

    msaa = msallocate(nh, ns);
    for (is=0; is < ns; is++) {
	for (ih=0; ih < nh; ih++) {
	    msaa->lag[is][ih] = aa->lag[ih]*jump[is]; /* expanded scale lags */
	}
    }
    deallocatehelix(aa);

    n123=1;
    for (id=0; id < ndim; id++) {
	n123 *= nd[id];
    }
    
    msaa->mis = sf_boolalloc2(n123,ns);
    aa = msaa->one;

    for (is=0; is < ns; is++) { /* for all scales */
	onescale( is, msaa); /* extract a filter */  
	aa->mis = NULL;
	for (id=0; id < ndim; id++) {
	    nb[id] = na[id]*jump[is];
	}
	bound(ndim, nd, nd, nb, aa); /* set up its boundaries */
	for (id=0; id < n123; id++) {
	    msaa->mis[is][id] = aa->mis[id];  /* save them */
	}
	free (aa->mis);
    }
 
    return msaa;
}

/* 	$Id: createmshelix.c,v 1.1 2004/06/11 10:51:33 fomels Exp $	 */
