#ifndef _offruffp_h
#define _offruffp_h

#include <rsf.h>

void offruffp_init (float h0, int nh, float dh,
		    int nx, float dx, float w, int num);
float complex offruffp_c1 (void);
float complex offruffp_c2 (void);
void offruffp_lop (bool adj, bool add, int n1, int n2, 
		   complex float *x, complex float *y);
void hderp_lop (bool adj, bool add, int n1, int n2, 
		complex float *x, complex float *y);

#endif
