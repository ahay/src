#ifndef _helix_h
#define _helix_h

#include <rsf.h>

typedef struct helixfilter {
    int nh;
    float* flt;
    int* lag;
    bool* mis;
} *filter;

filter allocatehelix( int nh);
void deallocatehelix( filter aa);

#endif

/* 	$Id: helix.h,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */
