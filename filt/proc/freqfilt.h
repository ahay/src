#ifndef _freqfilt_h
#define _freqfilt_h

#include <rsf.h>

void freqfilt_init(int n1, int nfft, int nw);
void freqfilt_set(float* filt);
void freqfilt_close(void);
void freqfilt_lop (bool adj, bool add, int nx, int ny, float* x, float* y);

#endif

/* 	$Id: freqfilt.h,v 1.1 2004/04/02 02:30:49 fomels Exp $	 */
