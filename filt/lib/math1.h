#ifndef _sf_math_h
#define _sf_math_h

#include "file.h"

#define SF_PI (3.141592653589793)

void sf_math_evaluate (int len, int nbuf, float** fbuf, float** fst);
int sf_math_parse (char* output, sf_file out);

#endif

/* 	$Id$	 */
