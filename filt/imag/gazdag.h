#ifndef _gazdag_h
#define _gazdag_h

#include <rsf.h>

void gazdag_init (float eps, int nt, float dt, 
		  int nz, float dz, float *vt, bool depth); /* , float *gt); */
void gazdag_close ();
void gazdag (bool inv, float k2, float *p, float *q);

#endif

/* 	$Id$	 */
