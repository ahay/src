#include <rsf.h>

#include "mshelix.h"
#include "helix.h"

msfilter msallocate(int nh, int ns) {
    msfilter aa;

    aa = (msfilter) sf_alloc(1,sizeof(*aa));
    aa->nh = nh;
    aa->ns = ns;
    aa->flt = sf_floatalloc(nh);
    aa->lag = sf_intalloc2(nh,ns);
    aa->mis = NULL;
    aa->one = (filter) sf_alloc(1,sizeof(*aa->one));
    aa->one->flt = aa->flt;
    aa->one->nh = nh;

    return aa;
}

void msdeallocate( msfilter aa) {
    free( aa->flt);
    free( aa->lag[0]);
    free( aa->lag);
    if (NULL != aa->mis) {
	free( aa->mis[0]);
	free( aa->mis);
    }
    free( aa->one);
    free( aa);
}

void onescale(int i, msfilter aa) {
    aa->one->lag = aa->lag[i];
    if (NULL != aa->mis) aa->one->mis = aa->mis[i];
}

/* 	$Id: mshelix.c,v 1.1 2004/06/11 10:51:34 fomels Exp $	 */
