#ifndef _explanesignoi_h
#define _explanesignoi_h

#include <rsf.h>

void explanesignoi_init (int m1,int m2, float eps, float **aa,
			 int nw, int nj1, int nj2, float **nn, float **ss);
void explanesignoi_close(void);
void explanesignoi_lop (bool adj, bool add, int ns, int nd, 
			float *sig, float *dat);

#endif
