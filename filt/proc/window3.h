#ifndef _window3_h
#define _window3_h

#include <rsf.h>

void window3_init (int* w_in, int* nw_in, int* n_in, float* h_in);
void window3_apply (const int* i, float*** dat, 
		    bool near, bool far, 
		    bool left, bool right, 
		    bool top, bool bottom,
		    int* i0, float*** win);
#endif

/* 	$Id: window3.h,v 1.1 2004/05/06 00:03:13 fomels Exp $	 */
