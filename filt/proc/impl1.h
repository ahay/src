#ifndef _impl1_h
#define _impl1_h

#include <rsf.h>

void impl1_init (float r, int n1, float tau, float pclip, bool up);
void impl1_close (void);
void impl1_apply (float *x);


#endif
