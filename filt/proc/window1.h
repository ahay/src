#ifndef _window1_h
#define _window1_h

void window1_init (float h_in, float dw_in);
int window1_apply (int iw, int w, const float* dat, 
		   bool left, bool right, float *win);

#endif

/* 	$Id$	 */
