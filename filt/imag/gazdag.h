#ifndef _gazdag_h
#define _gazdag_h

#include <rsf.h>

void gazdag_init (float eps, int nt, float dt, 
		  int nz, float dz, float *vt, bool depth); /* , float *gt); */
void gazdag_close ();
void gazdag (bool inv, float k2, float complex *p, float complex *q);

#endif

/* 	$Id: gazdag.h,v 1.4 2004/01/15 02:36:32 fomels Exp $	 */
