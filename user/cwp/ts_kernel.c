/* Finite-difference time-stepping kernel */
#ifndef _TS_KERNEL_H
#define _TS_KERNEL_H
#include <stddef.h>

#ifndef SF_HAS_SSE
#undef __SSE__
#endif
#ifndef SF_HAS_AVX
#undef __AVX__
#endif

#ifdef sun
#define restrict
#endif

#if defined __SSE__ || defined __AVX__
#include <immintrin.h>
#ifdef __SSE__
#define vType __m128
#define prefix _mm
#elif defined __AVX__
#define vType __m256
#define prefix _mm256
#endif
#if defined(__STDC__) || defined(__cplusplus) || defined(c_plusplus)
#define _STR_JOIN(a, b) a##b
#else
#define _STR_JOIN(a, b) a /**/ b
#endif
#define STR_JOIN(a, b) _STR_JOIN(a, b)
#define _prefix_loadu_ps STR_JOIN(prefix, _loadu_ps)
#define _prefix_storeu_ps STR_JOIN(prefix, _storeu_ps)
#define _prefix_add_ps STR_JOIN(prefix, _add_ps)
#define _prefix_sub_ps STR_JOIN(prefix, _sub_ps)
#define _prefix_mul_ps STR_JOIN(prefix, _mul_ps)
#define _prefix_div_ps STR_JOIN(prefix, _div_ps)
#define _prefix_set1_ps STR_JOIN(prefix, _set1_ps)
#endif

#ifdef __SSE__
static const int vSize = 4;
#elif defined __AVX__
static const int vSize = 8;
#else
static const int vSize = 1;
#endif

void
step_forward_2d(float** restrict u0, float** restrict u1, float** restrict vel,
                float** restrict rho, float* restrict fdcoef_d2,
                float* restrict fdcoef_d1, int nop, int nzpad, int nxpad)
/*< 2D FD time-stepping with SIMD vectorization support >*/
{
  int ix, iz, iop;
  float c0 = fdcoef_d2[0];
  float* cz = &fdcoef_d2[0];
  float* cx = &fdcoef_d2[nop];
  float *bz, *bx;
  float du_z = 0.f, du_x = 0.f;
  float drho_z = 0.f, drho_x = 0.f;
  float lap;
  float drho_dot_du;

  if (rho != NULL) { /* variable density */
    bz = &fdcoef_d1[0];
    bx = &fdcoef_d1[nop];
  }

#ifdef _OPENMP
#pragma omp parallel for \
  schedule(static, 1) \
  shared(nxpad, nzpad, u0, u1, vel, c0, cx, cz) \
  private(ix, iz, iop, lap, du_z, du_x, drho_z, drho_x, drho_dot_du)
#endif
  for (ix = nop; ix < nxpad - nop; ix++) {
    for (iz = nop; iz < nzpad - nop; iz += vSize) {
#if defined(__SSE__) || defined(__AVX__)
      /* load u0 u1 vel from memory to register */
      vType vec_u0 = _prefix_loadu_ps(&u0[ix][iz]);
      vType vec_u1 = _prefix_loadu_ps(&u1[ix][iz]);
      vec_u0 =
        _prefix_sub_ps(_prefix_mul_ps(_prefix_set1_ps(2.f), vec_u1), vec_u0);

      vType vec_lap = _prefix_mul_ps(vec_u1, _prefix_set1_ps(c0));
      for (iop = 1; iop <= nop; iop++) {
        /* z axis derivative <1st axis> */
        vec_lap = _prefix_add_ps(
          vec_lap,
          _prefix_mul_ps(_prefix_set1_ps(cz[iop]),
                         _prefix_add_ps(_prefix_loadu_ps(&u1[ix][iz + iop]),
                                        _prefix_loadu_ps(&u1[ix][iz - iop]))));
        /* x axis derivative <2nd axis> */
        vec_lap = _prefix_add_ps(
          vec_lap,
          _prefix_mul_ps(_prefix_set1_ps(cx[iop]),
                         _prefix_add_ps(_prefix_loadu_ps(&u1[ix + iop][iz]),
                                        _prefix_loadu_ps(&u1[ix - iop][iz]))));
      }

      if (rho != NULL) {
        vType vec_du_z = _prefix_set1_ps(0.f);
        vType vec_du_x = _prefix_set1_ps(0.f);
        vType vec_drho_z = _prefix_set1_ps(0.f);
        vType vec_drho_x = _prefix_set1_ps(0.f);
        for (iop = 1; iop <= nop; iop++) {
          vec_du_z = _prefix_add_ps(
            vec_du_z, _prefix_mul_ps(
                        _prefix_set1_ps(bz[iop]),
                        _prefix_sub_ps(_prefix_loadu_ps(&u1[ix][iz + iop]),
                                       _prefix_loadu_ps(&u1[ix][iz - iop]))));
          vec_drho_z = _prefix_add_ps(
            vec_drho_z,
            _prefix_mul_ps(
              _prefix_set1_ps(bz[iop]),
              _prefix_sub_ps(_prefix_loadu_ps(&rho[ix][iz + iop]),
                             _prefix_loadu_ps(&rho[ix][iz - iop]))));
          vec_du_x = _prefix_add_ps(
            vec_du_x, _prefix_mul_ps(
                        _prefix_set1_ps(bx[iop]),
                        _prefix_sub_ps(_prefix_loadu_ps(&u1[ix + iop][iz]),
                                       _prefix_loadu_ps(&u1[ix - iop][iz]))));
          vec_drho_x = _prefix_add_ps(
            vec_drho_x,
            _prefix_mul_ps(
              _prefix_set1_ps(bx[iop]),
              _prefix_sub_ps(_prefix_loadu_ps(&rho[ix + iop][iz]),
                             _prefix_loadu_ps(&rho[ix - iop][iz]))));
        }
        vec_lap = _prefix_sub_ps(
          vec_lap,
          _prefix_div_ps(_prefix_add_ps(_prefix_mul_ps(vec_du_z, vec_drho_z),
                                        _prefix_mul_ps(vec_du_x, vec_drho_x)),
                         _prefix_loadu_ps(&rho[ix][iz])));
      }
      vec_u0 = _prefix_add_ps(
        vec_u0, _prefix_mul_ps(_prefix_loadu_ps(&vel[ix][iz]), vec_lap));
      _prefix_storeu_ps(&u0[ix][iz], vec_u0);
#else
      lap = u1[ix][iz] * c0;
      for (iop = 1; iop <= nop; iop++) {
        lap += (u1[ix][iz - iop] + u1[ix][iz + iop]) * cz[iop] +
               (u1[ix - iop][iz] + u1[ix + iop][iz]) * cx[iop];
      }
      if (rho != NULL) { /* variable density term */
        du_z = du_x = drho_z = drho_x = 0.f;
        for (iop = 1; iop <= nop; iop++) {
          du_z += (u1[ix][iz + iop] - u1[ix][iz - iop]) * bz[iop];
          du_x += (u1[ix + iop][iz] - u1[ix - iop][iz]) * bx[iop];
          drho_z += (rho[ix][iz + iop] - rho[ix][iz - iop]) * bz[iop];
          drho_x += (rho[ix + iop][iz] - rho[ix - iop][iz]) * bx[iop];
        }
        drho_dot_du = (du_z * drho_z + du_x * drho_x) / rho[ix][iz];
        lap -= drho_dot_du;
      }
      u0[ix][iz] = 2. * u1[ix][iz] - u0[ix][iz] + vel[ix][iz] * lap;
#endif
    }
  }
  return;
}

void
step_forward_3d(float*** restrict u0, float*** restrict u1,
                float*** restrict vel, float*** restrict rho,
                float* restrict fdcoef_d2, float* restrict fdcoef_d1, int nop,
                int nzpad, int nxpad, int nypad)
/*< 3D FD time-stepping with SIMD vectorization support >*/
{
  int ix, iy, iz, iop;
  float c0 = fdcoef_d2[0];
  float* cz = &fdcoef_d2[0];
  float* cx = &fdcoef_d2[nop];
  float* cy = &fdcoef_d2[nop + nop];
  float *bz, *bx, *by;
  float drho_dot_du;
  float du_z = 0.f, du_x = 0.f, du_y = 0.f;
  float drho_z = 0.f, drho_x = 0.f, drho_y = 0.f;
  float lap;

  if (rho != NULL) { /* variable density */
    bz = &fdcoef_d1[0];
    bx = &fdcoef_d1[nop];
    by = &fdcoef_d1[nop + nop];
  }

#ifdef _OPENMP
#pragma omp parallel for \
  schedule(static, 1) \
  shared(nxpad, nypad, nzpad, u0, u1, vel, c0, cx, cy, cz) \
  private(ix, iy, iz, iop, drho_dot_du, du_z, du_x, du_y, drho_z, drho_x, drho_y, lap)
#endif
  for (iy = nop; iy < nypad - nop; iy++) {
    for (ix = nop; ix < nxpad - nop; ix++) {
      for (iz = nop; iz < nzpad - nop; iz += vSize) {
#if defined(__SSE__) || defined(__AVX__)
        /* load u0 u1 vel from memory to register */
        vType vec_u0 = _prefix_loadu_ps(&u0[iy][ix][iz]);
        vType vec_u1 = _prefix_loadu_ps(&u1[iy][ix][iz]);
        vec_u0 =
          _prefix_sub_ps(_prefix_mul_ps(_prefix_set1_ps(2.f), vec_u1), vec_u0);

        vType vec_lap = _prefix_mul_ps(vec_u1, _prefix_set1_ps(c0));
        for (iop = 1; iop <= nop; iop++) {
          /* z axis derivative <1st axis> */
          vec_lap = _prefix_add_ps(
            vec_lap,
            _prefix_mul_ps(
              _prefix_set1_ps(cz[iop]),
              _prefix_add_ps(_prefix_loadu_ps(&u1[iy][ix][iz + iop]),
                             _prefix_loadu_ps(&u1[iy][ix][iz - iop]))));
          /* x axis derivative <2nd axis> */
          vec_lap = _prefix_add_ps(
            vec_lap,
            _prefix_mul_ps(
              _prefix_set1_ps(cx[iop]),
              _prefix_add_ps(_prefix_loadu_ps(&u1[iy][ix + iop][iz]),
                             _prefix_loadu_ps(&u1[iy][ix - iop][iz]))));
          /* y axis derivative <3rd axis> */
          vec_lap = _prefix_add_ps(
            vec_lap,
            _prefix_mul_ps(
              _prefix_set1_ps(cy[iop]),
              _prefix_add_ps(_prefix_loadu_ps(&u1[iy + iop][ix][iz]),
                             _prefix_loadu_ps(&u1[iy - iop][ix][iz]))));
        }

        if (rho != NULL) {
          vType vec_du_z = _prefix_set1_ps(0.f);
          vType vec_du_x = _prefix_set1_ps(0.f);
          vType vec_du_y = _prefix_set1_ps(0.f);
          vType vec_drho_z = _prefix_set1_ps(0.f);
          vType vec_drho_x = _prefix_set1_ps(0.f);
          vType vec_drho_y = _prefix_set1_ps(0.f);
          for (iop = 1; iop <= nop; iop++) {
            vec_du_z = _prefix_add_ps(
              vec_du_z,
              _prefix_mul_ps(
                _prefix_set1_ps(bz[iop]),
                _prefix_sub_ps(_prefix_loadu_ps(&u1[iy][ix][iz + iop]),
                               _prefix_loadu_ps(&u1[iy][ix][iz - iop]))));
            vec_drho_z = _prefix_add_ps(
              vec_drho_z,
              _prefix_mul_ps(
                _prefix_set1_ps(bz[iop]),
                _prefix_sub_ps(_prefix_loadu_ps(&rho[iy][ix][iz + iop]),
                               _prefix_loadu_ps(&rho[iy][ix][iz - iop]))));
            vec_du_x = _prefix_add_ps(
              vec_du_x,
              _prefix_mul_ps(
                _prefix_set1_ps(bx[iop]),
                _prefix_sub_ps(_prefix_loadu_ps(&u1[iy][ix + iop][iz]),
                               _prefix_loadu_ps(&u1[iy][ix - iop][iz]))));
            vec_drho_x = _prefix_add_ps(
              vec_drho_x,
              _prefix_mul_ps(
                _prefix_set1_ps(bx[iop]),
                _prefix_sub_ps(_prefix_loadu_ps(&rho[iy][ix + iop][iz]),
                               _prefix_loadu_ps(&rho[iy][ix - iop][iz]))));
            vec_du_y = _prefix_add_ps(
              vec_du_y,
              _prefix_mul_ps(
                _prefix_set1_ps(by[iop]),
                _prefix_sub_ps(_prefix_loadu_ps(&u1[iy + iop][ix][iz]),
                               _prefix_loadu_ps(&u1[iy - iop][ix][iz]))));
            vec_drho_y = _prefix_add_ps(
              vec_drho_y,
              _prefix_mul_ps(
                _prefix_set1_ps(by[iop]),
                _prefix_sub_ps(_prefix_loadu_ps(&rho[iy + iop][ix][iz]),
                               _prefix_loadu_ps(&rho[iy - iop][ix][iz]))));
          }
          vec_lap = _prefix_sub_ps(
            vec_lap, _prefix_div_ps(
                       _prefix_add_ps(
                         _prefix_mul_ps(vec_du_z, vec_drho_z),
                         _prefix_add_ps(_prefix_mul_ps(vec_du_x, vec_drho_x),
                                        _prefix_mul_ps(vec_du_y, vec_drho_y))),
                       _prefix_loadu_ps(&rho[iy][ix][iz])));
        }
        vec_u0 = _prefix_add_ps(
          vec_u0, _prefix_mul_ps(_prefix_loadu_ps(&vel[iy][ix][iz]), vec_lap));
        _prefix_storeu_ps(&u0[iy][ix][iz], vec_u0);

#else
        lap = u1[iy][ix][iz] * c0;
        for (iop = 1; iop <= nop; iop++) {
          lap += (u1[iy][ix][iz - iop] + u1[iy][ix][iz + iop]) * cz[iop] +
                 (u1[iy][ix - iop][iz] + u1[iy][ix + iop][iz]) * cx[iop] +
                 (u1[iy - iop][ix][iz] + u1[iy + iop][ix][iz]) * cy[iop];
        }
        if (rho != NULL) { /* variable density term */
          du_z = du_x = du_y = drho_z = drho_x = drho_y = 0.f;
          for (iop = 1; iop <= nop; iop++) {
            du_z += (u1[iy][ix][iz + iop] - u1[iy][ix][iz - iop]) * bz[iop];
            du_x += (u1[iy][ix + iop][iz] - u1[iy][ix - iop][iz]) * bx[iop];
            du_y += (u1[iy + iop][ix][iz] - u1[iy - iop][ix][iz]) * by[iop];
            drho_z += (rho[iy][ix][iz + iop] - rho[iy][ix][iz - iop]) * bz[iop];
            drho_x += (rho[iy][ix + iop][iz] - rho[iy][ix - iop][iz]) * bx[iop];
            drho_y += (rho[iy + iop][ix][iz] - rho[iy - iop][ix][iz]) * by[iop];
          }
          drho_dot_du =
            (du_z * drho_z + du_x * drho_x + du_y * drho_y) / rho[iy][ix][iz];
          lap -= drho_dot_du;
        }
        u0[iy][ix][iz] =
          2. * u1[iy][ix][iz] - u0[iy][ix][iz] + vel[iy][ix][iz] * lap;
#endif
      }
    }
  }

  return;
}

#endif
