#ifndef _vp_gainpar_h
#define _vp_gainpar_h

#include <rsf.h>

void gainpar (sf_file in, float **data, int n1, int n2,int step,
	      float pclip,float phalf,
	      float *clip, float *gpow, float bias, int n3, int panel);

#endif
