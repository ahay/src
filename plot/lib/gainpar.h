#ifndef _vp_gainpar_h
#define _vp_gainpar_h

#include <rsf.h>

void gainpar (sf_file in, float **data, int n1, int n2,int step,
	      float tpow, float o1, float pclip,float phalf,
	      float *clip, float *gpow, float bias, float d1, 
	      int n3, int panel);

#endif
