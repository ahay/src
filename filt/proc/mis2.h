#ifndef _mis2_h
#define _mis2_h

#include <rsf.h>

#include "helix.h"

void mis2(int niter, int nx, float *xx, filter aa, 
	  const bool *known, bool doprec);

#endif

/* 	$Id: mis2.h,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */
