#ifndef _dsr_h
#define _dsr_h

void dsr_init (float eps1, int nt, float dt, 
	       int nz1, float dz1, float *vt1, bool depth1);
void dsr_close ();
void dsr (bool inv, float kx, float kh, float *p, float *q);

#endif

/* 	$Id$	 */
