#ifndef _fzero_h
#define _fzero_h

#include <rsf.h>

float fzero (float (*func)(float), 
	     float a, float b, float fa, float fb,
	     float toler, bool verb);

#endif

/* 	$Id: fzero.h,v 1.2 2003/09/30 14:30:52 fomels Exp $	 */
