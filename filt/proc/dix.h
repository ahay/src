#ifndef _dix_h
#define _dix_h

int dix_init(int ndim, int *rect, int *ndat);
void dix_close(void);
void dix(int niter, float* weight, float* vrms, float* vint);

#endif
