#ifndef _nhelix_h
#define _nhelix_h

#include <rsf.h>

#include "helix.h"

typedef struct nhelixfilter {
    int np;    
    filter* hlx;
    bool* mis;
    int *pch;
} *nfilter;

nfilter nallocate(int np, int nd, int *nh, int *pch);
void ndeallocate( nfilter aa);

#endif

/* 	$Id: nhelix.h,v 1.1 2004/06/18 01:06:45 fomels Exp $	 */
