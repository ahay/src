#ifndef _imospray_h
#define _imospray_h

#include <rsf.h>

void imospray_init (float slow, float y0, float dy, float z0, float dz, 
		    int nz, int ny);
void imospray_lop(bool adj, bool add, int n1, int n2, 
		   float *stack, float *gather);
void imospray_close(void);

#endif
