#ifndef _impl2_h
#define _impl2_h

#include <rsf.h>

void impl2_init (float r1, float r2, int n1, int n2, 
		 float tau, float pclip, bool up);
void impl2_close (void);
void impl2_apply (float **x);


#endif
