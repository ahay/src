#ifndef _radon_h
#define _radon_h

#include <rsf.h>

void radon_init (int nx_in, int np_in, float dp_in, float p0_in);
void radon_close ();
void radon_set (float w, float* xx);
void radon_toep (float complex *qq, float eps);
void radon_lop (bool adj, bool add, int nm, int nd, 
		float complex *mm, float complex *dd);

#endif

/* 	$Id: radon.h,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */
