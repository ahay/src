#include <rsf.h>

#include "helix.h"

filter allocatehelix( int nh) {
    filter aa;

    aa = (filter) sf_alloc(1,sizeof(*aa));
    aa->nh = nh;
    aa->flt = sf_floatalloc(nh);
    aa->lag = sf_intalloc(nh);
    aa->mis = NULL;
    
    return aa;
}

void deallocatehelix( filter aa) {
    free( aa->flt);
    free( aa->lag);
    if (NULL != aa->mis) free( aa->mis);
    free( aa);
}

/* 	$Id: helix.c,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */
