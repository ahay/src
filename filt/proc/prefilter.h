#ifndef _prefilter_h
#define _prefilter_h

/* 
   File: prefilter.h
   -----------------
   Convert data to B-spline coefficients
   by fast B-spline transform
*/

/*
  Function: prefilter_init
  ------------------------
  Initialize internal storage
  nw  - spline order (interpolator length)
  nt  - length for temporary storage
  pad - length for padding
*/
void prefilter_init (int nw, int nt, int pad);

/*
  Function: prefilter_apply
  -------------------------
  Convert 1-D data to spline coefficients
  nd      - data length
  dat[nd] - data/coefficients (input/output)
*/
void prefilter_apply (int nd, float* dat);

/*
  Function: prefilter
  -------------------
  Convert N-D data to spline coefficients
  dim    - dimensionality
  n[dim] - data dimensions
  dat    - data/coefficients as 1-D array (input/output)
*/
void prefilter (int dim, int* n, float* dat);

/*
  Function: prefilter_close
  -------------------------
  Free temporary storage
*/
void prefilter_close( void);

#endif
