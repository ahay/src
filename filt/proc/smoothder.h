#ifndef _smoothder_h
#define _smoothder_h

int smoothder_init(int ndim, int *rect, int *ndat);
void smoothder_close(void);
void smoothder(int niter, float* weight, float* data, float* der);

#endif
