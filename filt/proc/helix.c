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

