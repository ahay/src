/* Library for creating and evaluating B-splines */
/*

  Derived from einspline-0.9.2 - http://einspline.sf.net/
  Original code by Kenneth P. Esler, Jr.
  Licensed under GPL

*/

#ifndef EINSPLINE_H
#define EINSPLINE_H

#include <inttypes.h>

#if defined(__sun) && !defined(__GNUC__)
#define restrict _Restrict
#endif

/*
 * bspline_base.h
 */

typedef enum { PERIODIC, DERIV1, DERIV2, FLAT, NATURAL, ANTIPERIODIC } bc_code;
typedef enum { U1D       , U2D       , U3D      , 
               MULTI_U1D , MULTI_U2D , MULTI_U3D } spline_code;
typedef enum { SINGLE_REAL }
  type_code;

typedef struct 
{
  bc_code lCode, rCode;
  float lVal, rVal;
} BCtype_s;

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
  void *restrict coefs;
} Bspline;


void destroy_Bspline (void *spline);

/*************************************
 *************************************
 **        Uniform splines          **
 *************************************
 *************************************/

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
  float* restrict coefs;
  Ugrid x_grid;
  BCtype_s xBC;
  intptr_t nc;
} UBspline_1d_s;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  float* restrict coefs;
  int x_stride;
  Ugrid x_grid, y_grid;
  BCtype_s xBC, yBC;
  intptr_t nc;
} UBspline_2d_s;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  float* restrict coefs;
  int x_stride, y_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_s xBC, yBC, zBC;
  intptr_t nc;
} UBspline_3d_s;

/************************************************************/
/* 1D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_1d_s (UBspline_1d_s * restrict spline,
                         double x, float* restrict val);

/* Value and first derivative */
void eval_UBspline_1d_s_vg (UBspline_1d_s * restrict spline, double x,
                            float* restrict val, float* restrict grad);

/* Value, first derivative, and second derivative */
void eval_UBspline_1d_s_vgl (UBspline_1d_s * restrict spline, double x, 
                             float* restrict val, float* restrict grad,
                             float* restrict lapl);
void eval_UBspline_1d_s_vgh (UBspline_1d_s * restrict spline, double x, 
                             float* restrict val, float* restrict grad,
                             float* restrict hess);

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_2d_s (UBspline_2d_s * restrict spline, 
                         double x, double y, float* restrict val);

/* Value and gradient */
void eval_UBspline_2d_s_vg (UBspline_2d_s * restrict spline, 
                            double x, double y, 
                            float* restrict val, float* restrict grad);

/* Value, gradient, and laplacian */
void eval_UBspline_2d_s_vgl (UBspline_2d_s * restrict spline, 
                             double x, double y, float* restrict val, 
                             float* restrict grad, float* restrict lapl);

/* Value, gradient, and Hessian */
void eval_UBspline_2d_s_vgh (UBspline_2d_s * restrict spline, 
                             double x, double y, float* restrict val, 
                             float* restrict grad, float* restrict hess);

/************************************************************/
/* 3D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_3d_s (UBspline_3d_s * restrict spline, 
                         double x, double y, double z,
                         float* restrict val);

/* Value and gradient */
void eval_UBspline_3d_s_vg (UBspline_3d_s * restrict spline, 
                            double x, double y, double z,
                            float* restrict val, float* restrict grad);

/* Value, gradient, and laplacian */
void eval_UBspline_3d_s_vgl (UBspline_3d_s * restrict spline, 
                             double x, double y, double z,
                             float* restrict val, float* restrict grad, float* restrict lapl);

/* Value, gradient, and Hessian */
void eval_UBspline_3d_s_vgh (UBspline_3d_s * restrict spline, 
                             double x, double y, double z,
                             float* restrict val, float* restrict grad, float* restrict hess);

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
 **    Uniform, multiple splines    **
 *************************************
 *************************************/


/*
 * multiple_bspline_structs.h
 */

/*************************
 * Single precision real *
 *************************/

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  float* restrict coefs;
  intptr_t x_stride;
  Ugrid x_grid;
  BCtype_s xBC;
  int num_splines;
  intptr_t nc;
} multi_UBspline_1d_s;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  float* restrict coefs;
  intptr_t x_stride, y_stride;
  Ugrid x_grid, y_grid;
  BCtype_s xBC, yBC;
  int num_splines;
  intptr_t nc;
} multi_UBspline_2d_s;

typedef struct
{
  spline_code spcode;
  type_code    tcode;
  float* restrict coefs;
  intptr_t x_stride, y_stride, z_stride;
  Ugrid x_grid, y_grid, z_grid;
  BCtype_s xBC, yBC, zBC;
  int num_splines;
  intptr_t nc;
} multi_UBspline_3d_s;

/************************************************************/
/* 1D single-precision, real evaulation functions           */
/************************************************************/
void eval_multi_UBspline_1d_s (multi_UBspline_1d_s *spline,
                               double x,
                               float* restrict vals);

void eval_multi_UBspline_1d_s_vg (multi_UBspline_1d_s *spline,
                                  double x,
                                  float* restrict vals,
                                  float* restrict grads);

void eval_multi_UBspline_1d_s_vgl (multi_UBspline_1d_s *spline,
                                   double x,
                                   float* restrict vals,
                                   float* restrict grads,
                                   float* restrict lapl);

void eval_multi_UBspline_1d_s_vgh (multi_UBspline_1d_s *spline,
                                   double x,
                                   float* restrict vals,
                                   float* restrict grads,
                                   float* restrict hess);

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/
void eval_multi_UBspline_2d_s (multi_UBspline_2d_s *spline,
                               double x, double y,
                               float* restrict vals);

void eval_multi_UBspline_2d_s_vg (multi_UBspline_2d_s *spline,
                                  double x, double y,
                                  float* restrict vals,
                                  float* restrict grads);

void eval_multi_UBspline_2d_s_vgl (multi_UBspline_2d_s *spline,
                                   double x, double y,
                                   float* restrict vals,
                                   float* restrict grads,
                                   float* restrict lapl);

void eval_multi_UBspline_2d_s_vgh (multi_UBspline_2d_s *spline,
                                   double x, double y,
                                   float* restrict vals,
                                   float* restrict grads,
                                   float* restrict hess);

/************************************************************/
/* 3D single-precision, real evaulation functions           */
/************************************************************/
void eval_multi_UBspline_3d_s (multi_UBspline_3d_s *spline,
                               double x, double y, double z,
                               float* restrict vals);

void eval_multi_UBspline_3d_s_vg (multi_UBspline_3d_s *spline,
                                  double x, double y, double z,
                                  float* restrict vals,
                                  float* restrict grads);

void eval_multi_UBspline_3d_s_vgl (multi_UBspline_3d_s *spline,
                                   double x, double y, double z,
                                   float* restrict vals,
                                   float* restrict grads,
                                   float* restrict lapl);

void eval_multi_UBspline_3d_s_vgh (multi_UBspline_3d_s *spline,
                                   double x, double y, double z,
                                   float* restrict vals,
                                   float* restrict grads,
                                   float* restrict hess);

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
multi_UBspline_1d_s *
create_multi_UBspline_1d_s (Ugrid x_grid, BCtype_s xBC, int num_splines);

/* Create 2D uniform single-precision, real Bspline */
multi_UBspline_2d_s *
create_multi_UBspline_2d_s (Ugrid x_grid,   Ugrid y_grid,
                            BCtype_s   xBC, BCtype_s   yBC,
                            int num_splines);

/* Create 3D uniform single-precision, real Bspline */
multi_UBspline_3d_s *
create_multi_UBspline_3d_s (Ugrid x_grid,   Ugrid y_grid,   Ugrid z_grid,
                            BCtype_s  xBC,  BCtype_s   yBC, BCtype_s   zBC,
                            int num_splines);
  
/* Set the data for the splines, and compute spline coefficients */
void set_multi_UBspline_1d_s (multi_UBspline_1d_s *spline, 
                              int spline_num, float *data);

void set_multi_UBspline_2d_s (multi_UBspline_2d_s *spline, 
                              int spline_num, float *data);

void set_multi_UBspline_3d_s (multi_UBspline_3d_s *spline, 
                              int spline_num, float *data);

void init_einspline (void);

#endif /* EINSPLINE_H */

