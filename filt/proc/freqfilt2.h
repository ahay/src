#ifndef _freqfilt2_h
#define _freqfilt2_h

#include <rsf.h>

void freqfilt2_init(int n1, int n2, int nw);
void freqfilt2_set(float** filt);
void freqfilt2_close(void);
void freqfilt2_spec (const float* x, float** y);
void freqfilt2_lop (bool adj, bool add, int nx, int ny, float* x, float* y);

#endif

/* 	$Id$	 */
