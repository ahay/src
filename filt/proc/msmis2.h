#ifndef _msmis2_h
#define _msmis2_h

#include <rsf.h>

#include "mshelix.h"

void msmis2(int niter, int nx, int ns, float *xx, msfilter aa, 
	    const bool *known);

#endif

/* 	$Id: msmis2.h,v 1.1 2004/06/11 10:51:34 fomels Exp $	 */
