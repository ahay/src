#ifndef _gazdag_h
#define _gazdag_h

void gazdag_init (float eps, int nt, float dt, 
		  int nz, float dz, float *vt, float *gt);
void gazdag_close ();
void gazdag (bool inv, float k2, float complex *p, float complex *q);

#endif

/* 	$Id: gazdag.h,v 1.3 2003/11/22 22:42:28 fomels Exp $	 */
