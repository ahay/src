#ifndef _kolmog_h
#define _kolmog_h

#include <rsf.h>

void kolmog_init(int n1);
void kolmog_close(void);
void kolmog(float *trace);
void kolmog2(float *trace, float complex *fft);

#endif
