#ifndef _kolmog_h
#define _kolmog_h

#include <rsf.h>

void kolmog(int nfft, float *trace);
void kolmog2(int nfft, int nw, float *trace, float complex *fft);

#endif
