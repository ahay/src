#ifndef _smoothder_h
#define _smoothder_h

int smoothder_init(int ndim, int *rect, int *ndat, bool diff, bool dip);
void smoothder_close(void);
void smoothder(int niter, float* weight, float* data, float* der);
void smoothdiff(int niter, int ncycle, float* weight, float* data, float* der);
void smoothdip(int niter, float** dip, float* weight, float* data, float* der);

#endif
