#ifndef _split1_h
#define _split1_h

#include <rsf.h>

void split1 (bool verb, bool inv, float eps,  
	     int nt, float dt, 
	     int nz, float dz, 
	     int nx, float dx,
	     float **vt, float *v,
	     float complex **cp, float **q);

#endif

/* 	$Id: split1.h,v 1.2 2003/09/30 14:30:53 fomels Exp $	 */
