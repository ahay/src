#ifndef _freqfilt_h
#define _freqfilt_h

#include <rsf.h>

void freqfilt_init(int nfft, int nw);
void freqfilt_set(float* filt);
void freqfilt_close(void);
void freqfilt_lop (bool adj, bool add, int nx, int ny, float* x, float* y);

#endif

/* 	$Id$	 */
