#ifndef _int22_h
#define _int22_h

#include "bool.h"

typedef void (*interpolator)(float,int,float*);

void  int22_init (float** coord, float o[], float d[], int n[], 
		  interpolator interp, int nf_in, int nd_in);
void  int22_lop (bool adj, bool add, int nm, int nd, float* x, float* ord);
void int22_close (void);

#endif
