#ifndef _mshelix_h
#define _mshelix_h

#include <rsf.h>

#include "helix.h"

typedef struct mshelixfilter {
    int nh, ns;
    float* flt;
    int** lag;
    bool** mis;
    filter one;
} *msfilter;

msfilter msallocate(int nh, int ns);
void msdeallocate(msfilter msaa);
void onescale(int i, msfilter msaa); 

#endif

/* 	$Id: mshelix.h,v 1.1 2004/06/11 10:51:34 fomels Exp $	 */
