/* Muting the first arrival */

#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "muting.h"

int mutingf(int nt, int nx, float dt, float dx, float dz, int isx, int isz, int gpz, float vel, int wd, float **dat)
/*< muting function >*/
{
  float z2; /* source location */
  float tt;
  int *t;
  int ix,it,cut;

  t = sf_intalloc(nx);

  z2= powf((isz-gpz)*dz,2);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,tt,cut)
#endif
  for (ix = 0; ix < nx; ix++) {
    tt = sqrtf( powf((isx-ix)*dx,2) + z2 ) / vel;
    cut = (int) (tt/dt) + wd; /* why cannot use int to convert */
    if (cut > nt) cut = nt;
    t[ix] = cut;
  }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,it)
#endif
  for (ix = 0; ix < nx; ix++) {
    for (it = 0; it < t[ix]; it ++) {
      dat[ix][it] = 0.0f;
    }
  }

  return 0;
}

int mutingc(int nt, int nx, float dt, float dx, float dz, int isx, int isz, int gpz, float vel, int wd, sf_complex **dat)
/*< muting function >*/
{
  float z2; /* source location */
  float tt;
  int *t;
  int ix,it,cut;

  t = sf_intalloc(nx);

  z2= powf((isz-gpz)*dz,2);

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,tt,cut)
#endif
  for (ix = 0; ix < nx; ix++) {
    tt = sqrtf( powf((isx-ix)*dx,2) + z2 ) / vel;
    cut = (int) (tt/dt) + wd; /* why cannot use int to convert */
    if (cut > nt) cut = nt;
    t[ix] = cut;
  }

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(ix,it)
#endif
  for (ix = 0; ix < nx; ix++) {
    for (it = 0; it < t[ix]; it ++) {
      dat[ix][it] = sf_cmplx(0.0f,0.0f);
    }
  }

  return 0;
}
