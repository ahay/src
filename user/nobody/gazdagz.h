#ifndef _gazdagz_h
#define _gazdagz_h

#include <rsf.h>

void gazdagz_init (float eps, int nt, float dt, 
		  int nz, float dz, float *vt, bool depth, bool midpoint); 
void gazdagz_close ();
void gazdagz (bool inv, float k2, float complex *p, float complex *q);

#endif

/* 	$Id$	 */
