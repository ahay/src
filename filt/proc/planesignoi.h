#ifndef _planesignoi_h
#define _planesignoi_h

#include <rsf.h>

void planesignoi_init (int nw, int nj1, int nj2, int nx, int ny, 
		       float **nn, float **ss, float eps1);
void planesignoi_lop (bool adj, bool add, int ns, int nd, float *s, float *d);
void planesignoi_close(void);

#endif
