#ifndef _dsr_h
#define _dsr_h

void dsr (int inv, float eps, float kx, float kh, 
	  int nw, float dw, float fw, 
	  int nz, float dz, 
	  float *vt, float complex *p, float complex *q);

#endif
