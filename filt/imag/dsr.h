#ifndef _dsr_h
#define _dsr_h

void dsr (int inv, float eps, float kx, float kh, 
	  int nw, float dw, float fw, 
	  int nz, float dz, 
	  float *vt, float complex *p, float complex *q);

#endif

/* 	$Id: dsr.h,v 1.2 2003/09/30 14:30:52 fomels Exp $	 */
