#ifndef _expsignoi2_h
#define _expsignoi2_h

#include <rsf.h>

void expsignoi2_init (int m1,int m2, float eps, float **aa);
void expsignoi2_close(void);
void expsignoi2_lop (bool adj, bool add, int ns, int nd, 
		     float *sig, float *dat);

#endif
