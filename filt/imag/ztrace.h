#ifndef _ztrace_h
#define _ztrace_h

#define NS 3

void ztrace_init (int order, int iorder,
		  int nx, int nz, int np, int nt,
		  float dx, float dz, float dp,
		  float x0, float z0, float p0,
		  float** vel, float** slice_in);

void ztrace_close (void);

void ztrace_step (int kz);

#endif
