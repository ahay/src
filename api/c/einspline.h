/* Library for creating and evaluating B-splines */
/*
  Copyright (C) 2009 University of Texas at Austin
  Copyright (C) 2007 Kenneth P. Esler, Jr.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef EINSPLINE_H
#define EINSPLINE_H

/*
 * bspline_base.h
 */

#ifndef NO_COMPLEX
#ifndef SF_HAS_COMPLEX_H
#ifndef KISS_FFT_H
#ifdef __cplusplus
#include <complex>
typedef std::complex<float>  sf_complex;
typedef std::complex<double> sf_double_complex;
#else
#include <complex.h>
typedef complex float  sf_complex;
typedef complex double sf_double_complex;
#endif /* __cplusplus */
#endif /* KISS_FFT_H */
#endif /* SF_HAS_COMPLEX_H */
#endif /* NO_COMPLEX */

typedef enum { PERIODIC, DERIV1, DERIV2, FLAT, NATURAL, ANTIPERIODIC } bc_code;
typedef enum { U1D, U2D, U3D } spline_code;
typedef enum { SINGLE_REAL, DOUBLE_REAL, SINGLE_COMPLEX, DOUBLE_COMPLEX } type_code;

typedef struct 
{
  bc_code lCode, rCode;
  float lVal, rVal;
} BCtype_s;

typedef struct 
{
  bc_code lCode, rCode;
  double lVal, rVal;
} BCtype_d;

typedef struct 
{
  bc_code lCode, rCode;
  float lVal_r, lVal_i, rVal_r, rVal_i;
} BCtype_c;

typedef struct 
{
  bc_code lCode, rCode;
  double lVal_r, lVal_i, rVal_r, rVal_i;
} BCtype_z;

typedef struct
{
  double start, end;
  int num;

  /* private */
  double delta, delta_inv;
} Ugrid;

typedef struct
{
  spline_code sp_code;
  type_code   t_code;
  void *coefs;
} Bspline;

void destroy_Bspline (void *spline);

/*
 * bspline_structs.h
 */

/*************************
 * Single precision real *
 *************************/

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  float *coefs;
  Ugrid x_grid;
  BCtype_s xBC;
} UBspline_1d_s;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  float *coefs;
  int x_stride;
  Ugrid x_grid, y_grid;
  BCtype_s xBC, yBC;
} UBspline_2d_s;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  float *coefs;
  int x_stride, y_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_s xBC, yBC, zBC;
} UBspline_3d_s;

/*************************
 * Double precision real *
 *************************/

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  double *coefs;
  Ugrid x_grid;
  BCtype_d xBC;
} UBspline_1d_d;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  double *coefs;
  int x_stride;
  Ugrid x_grid, y_grid;
  BCtype_d xBC, yBC;
} UBspline_2d_d;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  double *coefs;
  int x_stride, y_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_d xBC, yBC, zBC;
} UBspline_3d_d;

#ifndef NO_COMPLEX
/****************************
 * Single precision complex *
 ****************************/

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  sf_complex *coefs;
  Ugrid x_grid;
  BCtype_c xBC;
} UBspline_1d_c;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  sf_complex *coefs;
  int x_stride;
  Ugrid x_grid, y_grid;
  BCtype_c xBC, yBC;
} UBspline_2d_c;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  sf_complex *coefs;
  int x_stride, y_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_c xBC, yBC, zBC;

} UBspline_3d_c;

/****************************
 * Double precision complex *
 ****************************/

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  sf_double_complex *coefs;
  Ugrid x_grid;
  BCtype_z xBC;
} UBspline_1d_z;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  sf_double_complex *coefs;
  int x_stride;
  Ugrid x_grid, y_grid;
  BCtype_z xBC, yBC;
} UBspline_2d_z;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  sf_double_complex *coefs;
  int x_stride, y_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_z xBC, yBC, zBC;
} UBspline_3d_z;

/*
 * bspline_eval_std_s|d|c|z.h
 */
#endif /* NO_COMPLEX */

/************************************************************/
/* 1D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_1d_s (UBspline_1d_s * spline, 
                         double x, float* val);

/* Value and first derivative */
void eval_UBspline_1d_s_vg (UBspline_1d_s * spline, double x, 
                            float* val, float* grad);

/* Value, first derivative, and second derivative */
void eval_UBspline_1d_s_vgl (UBspline_1d_s * spline, double x, 
                             float* val, float* grad,
                             float* lapl);
void eval_UBspline_1d_s_vgh (UBspline_1d_s * spline, double x, 
                             float* val, float* grad,
                             float* hess);

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_2d_s (UBspline_2d_s * spline, 
                         double x, double y, float* val);

/* Value and gradient */
void eval_UBspline_2d_s_vg (UBspline_2d_s * spline, 
                            double x, double y, 
                            float* val, float* grad);

/* Value, gradient, and laplacian */
void eval_UBspline_2d_s_vgl (UBspline_2d_s * spline, 
                             double x, double y, float* val, 
                             float* grad, float* lapl);

/* Value, gradient, and Hessian */
void eval_UBspline_2d_s_vgh (UBspline_2d_s * spline, 
                             double x, double y, float* val, 
                             float* grad, float* hess);

/************************************************************/
/* 3D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_3d_s (UBspline_3d_s * spline, 
                         double x, double y, double z,
                         float* val);

/* Value and gradient */
void eval_UBspline_3d_s_vg (UBspline_3d_s * spline, 
                            double x, double y, double z,
                            float* val, float* grad);

/* Value, gradient, and laplacian */
void eval_UBspline_3d_s_vgl (UBspline_3d_s * spline, 
                             double x, double y, double z,
                             float* val, float* grad, float* lapl);

/* Value, gradient, and Hessian */
void eval_UBspline_3d_s_vgh (UBspline_3d_s * spline, 
                             double x, double y, double z,
                             float* val, float* grad, float* hess);

/************************************************************/
/* 1D double-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_1d_d (UBspline_1d_d * spline, 
                         double x, double* val);

/* Value and first derivative */
void eval_UBspline_1d_d_vg (UBspline_1d_d * spline, double x, 
                            double* val, double* grad);

/* Value, first derivative, and second derivative */
void eval_UBspline_1d_d_vgl (UBspline_1d_d * spline, double x, 
                             double* val, double* grad,
                             double* lapl);

void eval_UBspline_1d_d_vgh (UBspline_1d_d * spline, double x, 
                             double* val, double* grad,
                             double* hess);

/************************************************************/
/* 2D double-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_2d_d (UBspline_2d_d * spline, 
                         double x, double y, double* val);

/* Value and gradient */
void eval_UBspline_2d_d_vg (UBspline_2d_d * spline, 
                            double x, double y, 
                            double* val, double* grad);

/* Value, gradient, and laplacian */
void eval_UBspline_2d_d_vgl (UBspline_2d_d * spline, 
                             double x, double y, double* val, 
                             double* grad, double* lapl);

/* Value, gradient, and Hessian */
void eval_UBspline_2d_d_vgh (UBspline_2d_d * spline, 
                             double x, double y, double* val, 
                             double* grad, double* hess);

/************************************************************/
/* 3D double-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_3d_d (UBspline_3d_d * spline, 
                         double x, double y, double z,
                         double* val);

/* Value and gradient */
void eval_UBspline_3d_d_vg (UBspline_3d_d * spline, 
                            double x, double y, double z,
                            double* val, double* grad);

/* Value, gradient, and laplacian */
void eval_UBspline_3d_d_vgl (UBspline_3d_d * spline, 
                             double x, double y, double z,
                             double* val, double* grad, double* lapl);

/* Value, gradient, and Hessian */
void eval_UBspline_3d_d_vgh (UBspline_3d_d * spline, 
                             double x, double y, double z,
                             double* val, double* grad, double* hess);

#ifndef NO_COMPLEX
/************************************************************/
/* 1D single-precision, complex evaulation functions        */
/************************************************************/

/* Value only */
void eval_UBspline_1d_c (UBspline_1d_c * spline, 
                         double x, sf_complex* val);

/* Value and gradient */
void eval_UBspline_1d_c_vg (UBspline_1d_c * spline, double x, 
                            sf_complex* val, sf_complex* grad);

/* Value, gradient, and laplacian */
void eval_UBspline_1d_c_vgl (UBspline_1d_c * spline, double x, 
                             sf_complex* val, sf_complex* grad,
                             sf_complex* lapl);

/* Value, gradient, and Hessian */
void eval_UBspline_1d_c_vgh (UBspline_1d_c * spline, double x, 
                             sf_complex* val, 
                             sf_complex* grad,
                             sf_complex* hess);

/************************************************************/
/* 2D single-precision, complex evaulation functions        */
/************************************************************/

/* Value only */
void eval_UBspline_2d_c (UBspline_2d_c * spline, 
                         double x, double y, sf_complex* val);

/* Value and gradient */
void eval_UBspline_2d_c_vg (UBspline_2d_c * spline, 
                            double x, double y, 
                            sf_complex* val, sf_complex* grad);

/* Value, gradient, and laplacian */
void eval_UBspline_2d_c_vgl (UBspline_2d_c * spline, 
                             double x, double y, sf_complex* val, 
                             sf_complex* grad, sf_complex* lapl);

/* Value, gradient, and Hessian */
void eval_UBspline_2d_c_vgh (UBspline_2d_c * spline, 
                             double x, double y, sf_complex* val, 
                             sf_complex* grad, sf_complex* hess);

/************************************************************/
/* 3D single-precision, complex evaulation functions        */
/************************************************************/

/* Value only */
void eval_UBspline_3d_c (UBspline_3d_c * spline, 
                         double x, double y, double z,
                         sf_complex* val);

/* Value and gradient */
void eval_UBspline_3d_c_vg (UBspline_3d_c * spline, 
                            double x, double y, double z,
                            sf_complex* val, sf_complex* grad);

/* Value, gradient, and laplacian */
void eval_UBspline_3d_c_vgl (UBspline_3d_c * spline, 
                             double x, double y, double z,
                             sf_complex* val, sf_complex* grad, 
                             sf_complex* lapl);

/* Value, gradient, and Hessian */
void eval_UBspline_3d_c_vgh (UBspline_3d_c * spline, 
                             double x, double y, double z,
                             sf_complex* val, sf_complex* grad, 
                             sf_complex* hess);

/************************************************************/
/* 1D double-precision, complex evaulation functions        */
/************************************************************/

/* Value only */
void eval_UBspline_1d_z (UBspline_1d_z * spline, 
                         double x, sf_double_complex* val);

/* Value and gradient */
void eval_UBspline_1d_z_vg (UBspline_1d_z * spline, double x, 
                            sf_double_complex* val, 
                            sf_double_complex* grad);

/* Value, gradient, and laplacian */
void eval_UBspline_1d_z_vgl (UBspline_1d_z * spline, double x, 
                             sf_double_complex* val, sf_double_complex* grad,
                             sf_double_complex* lapl);

/* Value, gradient, and Hessian */
void eval_UBspline_1d_z_vgh (UBspline_1d_z * spline, double x, 
                             sf_double_complex* val, 
                             sf_double_complex* grad,
                             sf_double_complex* hess);

/************************************************************/
/* 2D double-precision, complex evaulation functions        */
/************************************************************/

/* Value only */
void eval_UBspline_2d_z (UBspline_2d_z * spline, 
                         double x, double y, sf_double_complex* val);

/* Value and gradient */
void eval_UBspline_2d_z_vg (UBspline_2d_z * spline, 
                            double x, double y, 
                            sf_double_complex* val, 
                            sf_double_complex* grad);

/* Value, gradient, and laplacian */
void eval_UBspline_2d_z_vgl (UBspline_2d_z * spline, 
                             double x, double y, sf_double_complex* val, 
                             sf_double_complex* grad, 
                             sf_double_complex* lapl);

/* Value, gradient, and Hessian */
void eval_UBspline_2d_z_vgh (UBspline_2d_z * spline, 
                             double x, double y, sf_double_complex* val, 
                             sf_double_complex* grad, 
                             sf_double_complex* hess);

/************************************************************/
/* 3D double-precision, complex evaulation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_3d_z (UBspline_3d_z * spline, 
                         double x, double y, double z,
                         sf_double_complex* val);

/* Value and gradient */
void eval_UBspline_3d_z_vg (UBspline_3d_z * spline, 
                            double x, double y, double z,
                            sf_double_complex* val, 
                            sf_double_complex* grad);

/* Value, gradient, and laplacian */
void eval_UBspline_3d_z_vgl (UBspline_3d_z * spline, 
                             double x, double y, double z,
                             sf_double_complex* val, 
                             sf_double_complex* grad, 
                             sf_double_complex* lapl);

/* Value, gradient, and Hessian */
void eval_UBspline_3d_z_vgh (UBspline_3d_z * spline, 
                             double x, double y, double z,
                             sf_double_complex* val, 
                             sf_double_complex* grad, 
                             sf_double_complex* hess);
#endif /* NO_COMPLEX */

/*
 * bspline_create.h
 */

/***********************************************************
 ***********************************************************
 ****             Spline creation functions             ****
 ***********************************************************
 ***********************************************************/

/*************************************
 *************************************
 ** Uniform, single precision, real **
 *************************************
 *************************************/

/* Create 1D uniform single-precision, real Bspline */
UBspline_1d_s* create_UBspline_1d_s (Ugrid x_grid, BCtype_s xBC, float *data);

/* Create 2D uniform single-precision, real Bspline */
UBspline_2d_s* create_UBspline_2d_s (Ugrid x_grid,   Ugrid y_grid,
                                     BCtype_s   xBC, BCtype_s   yBC,
                                     float *data);

/* Create 3D uniform single-precision, real Bspline */
UBspline_3d_s* create_UBspline_3d_s (Ugrid x_grid,   Ugrid y_grid,   Ugrid z_grid,
                                     BCtype_s  xBC,  BCtype_s   yBC, BCtype_s   zBC,
                                     float *data);

void recompute_UBspline_1d_s (UBspline_1d_s* spline, float *data);

void recompute_UBspline_2d_s (UBspline_2d_s* spline, float *data);

void recompute_UBspline_3d_s (UBspline_3d_s* spline, float *data);

/*************************************
 *************************************
 ** Uniform, double precision, real **
 *************************************
 *************************************/

/* Create 1D uniform single-precision, real Bspline */
UBspline_1d_d* create_UBspline_1d_d (Ugrid x_grid, BCtype_d xBC, double *data);

/* Create 2D uniform single-precision, real Bspline */
UBspline_2d_d* create_UBspline_2d_d (Ugrid x_grid,   Ugrid y_grid,
                                     BCtype_d   xBC, BCtype_d   yBC,
                                     double *data);

/* Create 3D uniform single-precision, real Bspline */
UBspline_3d_d* create_UBspline_3d_d (Ugrid x_grid,   Ugrid   y_grid,   Ugrid z_grid,
                                     BCtype_d  xBC,  BCtype_d   yBC, BCtype_d   zBC,
                                     double *data);

void recompute_UBspline_1d_d (UBspline_1d_d* spline, double *data);

void recompute_UBspline_2d_d (UBspline_2d_d* spline, double *data);

void recompute_UBspline_3d_d (UBspline_3d_d* spline, double *data);

#ifndef NO_COMPLEX
/****************************************
 ****************************************
 ** Uniform, single precision, complex **
 ****************************************
 ****************************************/

/* Create 1D uniform single-precision, real Bspline */
UBspline_1d_c* create_UBspline_1d_c (Ugrid x_grid, BCtype_c xBC, sf_complex *data);

/* Create 2D uniform single-precision, real Bspline */
UBspline_2d_c* create_UBspline_2d_c (Ugrid   x_grid, Ugrid   y_grid,
                                     BCtype_c   xBC, BCtype_c   yBC,
                                     sf_complex *data);

/* Create 3D uniform single-precision, real Bspline */
UBspline_3d_c* create_UBspline_3d_c (Ugrid  x_grid, Ugrid y_grid, Ugrid z_grid,
                                     BCtype_c  xBC, BCtype_c yBC, BCtype_c zBC,
                                     sf_complex *data);

void recompute_UBspline_1d_c (UBspline_1d_c* spline, sf_complex *data);

void recompute_UBspline_2d_c (UBspline_2d_c* spline, sf_complex *data);

void recompute_UBspline_3d_c (UBspline_3d_c* spline, sf_complex *data);

/****************************************
 ****************************************
 ** Uniform, double precision, complex **
 ****************************************
 ****************************************/

/* Create 1D uniform double-precision, complex Bspline */
UBspline_1d_z* create_UBspline_1d_z (Ugrid x_grid, BCtype_z xBC, sf_double_complex *data);

/* Create 2D uniform double-precision, complex Bspline */
UBspline_2d_z* create_UBspline_2d_z (Ugrid x_grid, Ugrid y_grid,
                                     BCtype_z   xBC, BCtype_z   yBC,
                                     sf_double_complex *data);

/* Create 3D uniform double-precision, complex Bspline */
UBspline_3d_z* create_UBspline_3d_z (Ugrid  x_grid, Ugrid   y_grid, Ugrid z_grid,
                                     BCtype_z  xBC, BCtype_z   yBC, BCtype_z zBC,
                                     sf_double_complex *data);

void recompute_UBspline_1d_z (UBspline_1d_z* spline, sf_double_complex *data);

void recompute_UBspline_2d_z (UBspline_2d_z* spline, sf_double_complex *data);

void recompute_UBspline_3d_z (UBspline_3d_z* spline, sf_double_complex *data);
#endif /* NO_COMPLEX */

#endif /* EINSPLINE_H */
