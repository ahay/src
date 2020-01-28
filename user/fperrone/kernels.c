#include <rsf.h>

#ifndef _KERNELS_H

typedef struct wfl_struct wfl_struct_t;
/*^*/

typedef struct acq_struct acq_struct_t;
/*^*/

typedef struct mod_struct mod_struct_t;
/*^*/

struct wfl_struct{
  float *pc;
  float *pp;
  float *pa;
  int n1;
  int n2;
  float d1;
  float d2;
  float o1;
  float o2;
};
/*^*/

struct acq_struct{
  int ns;
  float *nr;
  float *xs;
  float *ys;
  float *zs;
  float *xr;
  float *yr;
  float *zr;
};
/*^*/

struct mod_struct{
  float *vmod;
  float *dmod;
};
/*^*/

#endif

void fwdextrap2d(wfl_struct_t * wfl, acq_struct_t const * acq, mod_struct_t const * mod)
/*< extrapolation kernel 2d - forward operator >*/
{

}

void adjextrap2d(wfl_struct_t * wfl, acq_struct_t const * acq, mod_struct_t const * mod)
/*< extrapolation kernel 2d - adjoint operator  >*/
{

}


