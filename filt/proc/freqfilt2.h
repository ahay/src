#ifndef _freqfilt2_h
#define _freqfilt2_h

#include <rsf.h>

void freqfilt2_init(int n1, int n2, int nfft, int nw, int nk);
void freqfilt2_set(float** filt);
void freqfilt2_close(void);
void freqfilt2_spec (const float* x, float** y);
void freqfilt2_lop (bool adj, bool add, int nx, int ny, float* x, float* y);
void twofreqfilt2_lop (bool adj, bool add, int nx, int ny, float* x, float* y);

#endif

/* 	$Id: freqfilt2.h,v 1.2 2004/02/27 20:59:57 fomels Exp $	 */
