#ifndef _window2_h
#define _window2_h

#include <rsf.h>

void window2_init (int* w_in, int* nw_in, int* n_in, float* h_in);
void window2_apply (const int* i, float** dat, 
		    bool left, bool right, bool top, bool bottom,
		    int* i0, float** win);
#endif

/* 	$Id: window2.h,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */
