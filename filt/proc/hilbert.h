#ifndef _hilbert_h
#define _hilbert_h

void hilbert_init(int nt1, int n1, float c1);
void hilbert_free(void);
void hilbert (const float* t1, float *t2);
void hilbert4 (const float* t1, float *t2);
void deriv (const float* trace, float* trace2);

#endif
