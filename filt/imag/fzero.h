#ifndef _fzero_h
#define _fzero_h

#include <rsf.h>

float fzero (float (*func)(float), 
	     float a, float b, float fa, float fb,
	     float toler, bool verb);

#endif
