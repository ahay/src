/* Library for creating and evaluating B-splines */
/*

  Derived from einspline-0.9.2 - http://einspline.sf.net/
  Original code by Kenneth P. Esler, Jr.
  Licensed under GPL

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "einspline.h"

/*
 * bspline_data.c
 */

/********************
 * Single precision *
 ********************/
const float A44f[16] = 
  { -1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
     3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
    -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
     1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0 };
const float* restrict Af = A44f;

const float dA44f[16] =
  {  0.0, -0.5,  1.0, -0.5,
     0.0,  1.5, -2.0,  0.0,
     0.0, -1.5,  1.0,  0.5,
     0.0,  0.5,  0.0,  0.0 };
const float* restrict dAf = dA44f;

const float d2A44f[16] = 
  {  0.0, 0.0, -1.0,  1.0,
     0.0, 0.0,  3.0, -2.0,
     0.0, 0.0, -3.0,  1.0,
     0.0, 0.0,  1.0,  0.0 };
const float* restrict d2Af = d2A44f;

/*
 * bspline_create.c
 */

static void solve_deriv_interp_1d_s (float bands[], float coefs[],
                                      int M, int cstride) {
  int row;
  /* Solve interpolating equations
     First and last rows are different */
  bands[4*(0)+1] /= bands[4*(0)+0];
  bands[4*(0)+2] /= bands[4*(0)+0];
  bands[4*(0)+3] /= bands[4*(0)+0];
  bands[4*(0)+0] = 1.0;
  bands[4*(1)+1] -= bands[4*(1)+0]*bands[4*(0)+1];
  bands[4*(1)+2] -= bands[4*(1)+0]*bands[4*(0)+2];
  bands[4*(1)+3] -= bands[4*(1)+0]*bands[4*(0)+3];
  bands[4*(0)+0] = 0.0;
  bands[4*(1)+2] /= bands[4*(1)+1];
  bands[4*(1)+3] /= bands[4*(1)+1];
  bands[4*(1)+1] = 1.0;

  /* Now do rows 2 through M+1 */
  for (row=2; row < (M+1); row++) {
    bands[4*(row)+1] -= bands[4*(row)+0]*bands[4*(row-1)+2];
    bands[4*(row)+3] -= bands[4*(row)+0]*bands[4*(row-1)+3];
    bands[4*(row)+2] /= bands[4*(row)+1];
    bands[4*(row)+3] /= bands[4*(row)+1];
    bands[4*(row)+0] = 0.0;
    bands[4*(row)+1] = 1.0;
  }

  /* Do last row */
  bands[4*(M+1)+1] -= bands[4*(M+1)+0]*bands[4*(M-1)+2];
  bands[4*(M+1)+3] -= bands[4*(M+1)+0]*bands[4*(M-1)+3];
  bands[4*(M+1)+2] -= bands[4*(M+1)+1]*bands[4*(M)+2];
  bands[4*(M+1)+3] -= bands[4*(M+1)+1]*bands[4*(M)+3];
  bands[4*(M+1)+3] /= bands[4*(M+1)+2];
  bands[4*(M+1)+2] = 1.0;

  coefs[(M+1)*cstride] = bands[4*(M+1)+3];
  /* Now back substitute up */
  for (row=M; row>0; row--)
    coefs[row*cstride] = bands[4*(row)+3] - bands[4*(row)+2]*coefs[cstride*(row+1)];

  /* Finish with first row */
  coefs[0] = bands[4*(0)+3] - bands[4*(0)+1]*coefs[1*cstride] - bands[4*(0)+2]*coefs[2*cstride];
}

/* On input, bands should be filled with:
   row 0   :  abcdInitial from boundary conditions
   rows 1:M:  basis functions in first 3 cols, data in last
   row M+1 :  abcdFinal   from boundary conditions
   cstride gives the stride between values in coefs.
   On exit, coefs with contain interpolating B-spline coefs */
static void solve_periodic_interp_1d_s (float bands[], float coefs[],
                                        int M, int cstride) {
  float *lastCol;
  int row;

  lastCol = (float*) malloc(M*sizeof(float));

  /* Now solve:
     First and last rows are different */
  bands[4*(0)+2] /= bands[4*(0)+1];
  bands[4*(0)+0] /= bands[4*(0)+1];
  bands[4*(0)+3] /= bands[4*(0)+1];
  bands[4*(0)+1]  = 1.0;
  bands[4*(M-1)+1] -= bands[4*(M-1)+2]*bands[4*(0)+0];
  bands[4*(M-1)+3] -= bands[4*(M-1)+2]*bands[4*(0)+3];
  bands[4*(M-1)+2]  = -bands[4*(M-1)+2]*bands[4*(0)+2];
  lastCol[0] = bands[4*(0)+0];

  for (row=1; row < (M-1); row++) {
    bands[4*(row)+1] -= bands[4*(row)+0] * bands[4*(row-1)+2];
    bands[4*(row)+3] -= bands[4*(row)+0] * bands[4*(row-1)+3];
    lastCol[row]   = -bands[4*(row)+0] * lastCol[row-1];
    bands[4*(row)+0] = 0.0;
    bands[4*(row)+2] /= bands[4*(row)+1];
    bands[4*(row)+3] /= bands[4*(row)+1];
    lastCol[row]  /= bands[4*(row)+1];
    bands[4*(row)+1]  = 1.0;
    if (row < (M-2)) {
      bands[4*(M-1)+3] -= bands[4*(M-1)+2]*bands[4*(row)+3];
      bands[4*(M-1)+1] -= bands[4*(M-1)+2]*lastCol[row];
      bands[4*(M-1)+2] = -bands[4*(M-1)+2]*bands[4*(row)+2];
    }
  }

  /* Now do last row
     The [2] element and [0] element are now on top of each other */
  bands[4*(M-1)+0] += bands[4*(M-1)+2];
  bands[4*(M-1)+1] -= bands[4*(M-1)+0] * (bands[4*(M-2)+2]+lastCol[M-2]);
  bands[4*(M-1)+3] -= bands[4*(M-1)+0] *  bands[4*(M-2)+3];
  bands[4*(M-1)+3] /= bands[4*(M-1)+1];
  coefs[M*cstride] = bands[4*(M-1)+3];
  for (row=M-2; row>=0; row--) 
    coefs[(row+1)*cstride] = 
      bands[4*(row)+3] - bands[4*(row)+2]*coefs[(row+2)*cstride] - lastCol[row]*coefs[M*cstride];

  coefs[0*cstride] = coefs[M*cstride];
  coefs[(M+1)*cstride] = coefs[1*cstride];
  coefs[(M+2)*cstride] = coefs[2*cstride];

  free(lastCol);
}

/* On input, bands should be filled with:
   row 0   :  abcdInitial from boundary conditions
   rows 1:M:  basis functions in first 3 cols, data in last
   row M+1 :  abcdFinal   from boundary conditions
   cstride gives the stride between values in coefs.
   On exit, coefs with contain interpolating B-spline coefs */
static void solve_antiperiodic_interp_1d_s (float bands[], float coefs[],
                                            int M, int cstride)
{
  bands[4*0+0]     *= -1.0;
  bands[4*(M-1)+2] *= -1.0;
  int row;

  float lastCol[M];
  /* Now solve:
     First and last rows are different */
  bands[4*(0)+2] /= bands[4*(0)+1];
  bands[4*(0)+0] /= bands[4*(0)+1];
  bands[4*(0)+3] /= bands[4*(0)+1];
  bands[4*(0)+1]  = 1.0;
  bands[4*(M-1)+1] -= bands[4*(M-1)+2]*bands[4*(0)+0];
  bands[4*(M-1)+3] -= bands[4*(M-1)+2]*bands[4*(0)+3];
  bands[4*(M-1)+2]  = -bands[4*(M-1)+2]*bands[4*(0)+2];
  lastCol[0] = bands[4*(0)+0];

  for (row=1; row < (M-1); row++) {
    bands[4*(row)+1] -= bands[4*(row)+0] * bands[4*(row-1)+2];
    bands[4*(row)+3] -= bands[4*(row)+0] * bands[4*(row-1)+3];
    lastCol[row]   = -bands[4*(row)+0] * lastCol[row-1];
    bands[4*(row)+0] = 0.0;
    bands[4*(row)+2] /= bands[4*(row)+1];
    bands[4*(row)+3] /= bands[4*(row)+1];
    lastCol[row]  /= bands[4*(row)+1];
    bands[4*(row)+1]  = 1.0;
    if (row < (M-2)) {
      bands[4*(M-1)+3] -= bands[4*(M-1)+2]*bands[4*(row)+3];
      bands[4*(M-1)+1] -= bands[4*(M-1)+2]*lastCol[row];
      bands[4*(M-1)+2] = -bands[4*(M-1)+2]*bands[4*(row)+2];
    }
  }

  /* Now do last row
     The [2] element and [0] element are now on top of each other */
  bands[4*(M-1)+0] += bands[4*(M-1)+2];
  bands[4*(M-1)+1] -= bands[4*(M-1)+0] * (bands[4*(M-2)+2]+lastCol[M-2]);
  bands[4*(M-1)+3] -= bands[4*(M-1)+0] *  bands[4*(M-2)+3];
  bands[4*(M-1)+3] /= bands[4*(M-1)+1];
  coefs[M*cstride] = bands[4*(M-1)+3];
  for (row=M-2; row>=0; row--) 
    coefs[(row+1)*cstride] = 
      bands[4*(row)+3] - bands[4*(row)+2]*coefs[(row+2)*cstride] - lastCol[row]*coefs[M*cstride];
  
  coefs[0*cstride]     = -coefs[M*cstride];
  coefs[(M+1)*cstride] = -coefs[1*cstride];
  coefs[(M+2)*cstride] = -coefs[2*cstride];


}

static void find_coefs_1d_s (Ugrid grid, BCtype_s bc, 
                             float *data,  intptr_t dstride,
                             float *coefs, intptr_t cstride) {
  int M = grid.num, i, j;
  float basis[4] = {1.0/6.0, 2.0/3.0, 1.0/6.0, 0.0};
  float *bands;

  if (bc.lCode == PERIODIC || bc.lCode == ANTIPERIODIC) {
      bands = (float*) malloc(4*M*sizeof(float));
      for (i=0; i<M; i++) {
          bands[4*i+0] = basis[0];
          bands[4*i+1] = basis[1];
          bands[4*i+2] = basis[2];
          bands[4*i+3] = data[i*dstride];
      }
      if (bc.lCode == PERIODIC)
        solve_periodic_interp_1d_s (bands, coefs, M, cstride);
      else
        solve_antiperiodic_interp_1d_s (bands, coefs, M, cstride);
      free (bands);
  } else {
    /* Setup boundary conditions */
    float abcd_left[4], abcd_right[4];
    /* Left boundary */
    if (bc.lCode == FLAT || bc.lCode == NATURAL)
      bc.lVal = 0.0;
    if (bc.lCode == FLAT || bc.lCode == DERIV1) {
      abcd_left[0] = -0.5 * grid.delta_inv;
      abcd_left[1] =  0.0 * grid.delta_inv; 
      abcd_left[2] =  0.5 * grid.delta_inv;
      abcd_left[3] =  bc.lVal;
    }
    if (bc.lCode == NATURAL || bc.lCode == DERIV2) {
      abcd_left[0] = 1.0 * grid.delta_inv * grid.delta_inv;
      abcd_left[1] =-2.0 * grid.delta_inv * grid.delta_inv;
      abcd_left[2] = 1.0 * grid.delta_inv * grid.delta_inv;
      abcd_left[3] = bc.lVal;
    }

    /* Right boundary */
    if (bc.rCode == FLAT || bc.rCode == NATURAL)
      bc.rVal = 0.0;
    if (bc.rCode == FLAT || bc.rCode == DERIV1) {
      abcd_right[0] = -0.5 * grid.delta_inv;
      abcd_right[1] =  0.0 * grid.delta_inv; 
      abcd_right[2] =  0.5 * grid.delta_inv;
      abcd_right[3] =  bc.rVal;
    }
    if (bc.rCode == NATURAL || bc.rCode == DERIV2) {
      abcd_right[0] = 1.0 *grid.delta_inv * grid.delta_inv;
      abcd_right[1] =-2.0 *grid.delta_inv * grid.delta_inv;
      abcd_right[2] = 1.0 *grid.delta_inv * grid.delta_inv;
      abcd_right[3] = bc.rVal;
    }
    bands = (float*) malloc ((M+2)*4*sizeof(float));
    for (i=0; i<4; i++) {
      bands[4*( 0 )+i]   = abcd_left[i];
      bands[4*(M+1)+i] = abcd_right[i];
    }
    for (i=0; i<M; i++) {
      for (j=0; j<3; j++)
        bands[4*(i+1)+j] = basis[j];
      bands[4*(i+1)+3] = data[i*dstride];
    }   
    /* Now, solve for coefficients */
    solve_deriv_interp_1d_s (bands, coefs, M, cstride);
    free (bands);
  }
}

/***********************************************************

        Single-Precision, Real Creation Routines

 ***********************************************************

  On input, bands should be filled with:
  row 0   :  abcdInitial from boundary conditions
  rows 1:M:  basis functions in first 3 cols, data in last
  row M+1 :  abcdFinal   from boundary conditions
  cstride gives the stride between values in coefs.
  On exit, coefs with contain interpolating B-spline coefs */

UBspline_1d_s* create_UBspline_1d_s (Ugrid x_grid, BCtype_s xBC, float *data) {
  /* Create new spline */
  int M, N;
  UBspline_1d_s* restrict spline = malloc (sizeof(UBspline_1d_s));
  spline->spcode = U1D;
  spline->tcode  = SINGLE_REAL;
  spline->xBC = xBC; spline->x_grid = x_grid;

  /* Setup internal variables */
  M = x_grid.num;

  if (xBC.lCode == PERIODIC || xBC.lCode == ANTIPERIODIC) {
    x_grid.delta     = (x_grid.end-x_grid.start)/(double)(x_grid.num);
    N = M+3;
  }
  else {
    x_grid.delta     = (x_grid.end-x_grid.start)/(double)(x_grid.num-1);
    N = M+2;
  }

  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;
  spline->coefs = malloc ((size_t)sizeof(float)*(size_t)N);
  spline->nc = (size_t)sizeof(float)*(size_t)N;
  find_coefs_1d_s (spline->x_grid, xBC, data, 1, spline->coefs, 1);

  return spline;
}

void recompute_UBspline_1d_s (UBspline_1d_s* spline, float *data) {
  find_coefs_1d_s (spline->x_grid, spline->xBC, data, 1, spline->coefs, 1);
}

UBspline_2d_s* create_UBspline_2d_s (Ugrid x_grid, Ugrid y_grid,
                                     BCtype_s xBC, BCtype_s yBC, float *data) {
  /* Create new spline */
  int Mx, My;
  int Nx, Ny, iy, ix;

  UBspline_2d_s* restrict spline = malloc (sizeof(UBspline_2d_s));
  spline->spcode = U2D;
  spline->tcode  = SINGLE_REAL;
  spline->xBC = xBC; 
  spline->yBC = yBC; 
  /* Setup internal variables */
  Mx = x_grid.num;
  My = y_grid.num;


  if (xBC.lCode == PERIODIC || xBC.lCode == ANTIPERIODIC)
    Nx = Mx+3;
  else
    Nx = Mx+2;
  x_grid.delta = (x_grid.end - x_grid.start)/(double)(Nx-3);
  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC || yBC.lCode == ANTIPERIODIC)
    Ny = My+3;
  else
    Ny = My+2;
  y_grid.delta = (y_grid.end - y_grid.start)/(double)(Ny-3);
  y_grid.delta_inv = 1.0/y_grid.delta;
  spline->y_grid   = y_grid;
  spline->x_stride = Ny;
  spline->coefs = malloc ((size_t)sizeof(float)*(size_t)Nx*(size_t)Ny);
  spline->nc = (size_t)sizeof(float)*(size_t)Nx*(size_t)Ny;

  /* First, solve in the X-direction */
  for (iy=0; iy<My; iy++) {
    intptr_t doffset = iy;
    intptr_t coffset = iy;
    find_coefs_1d_s (spline->x_grid, spline->xBC, data+doffset, My,
                     spline->coefs+coffset, Ny);
  }

  /* Now, solve in the Y-direction */
  for (ix=0; ix<Nx; ix++) {
    intptr_t doffset = ix*Ny;
    intptr_t coffset = ix*Ny;
    find_coefs_1d_s (spline->y_grid, spline->yBC, spline->coefs+doffset, 1, 
                     spline->coefs+coffset, 1);
  }
  return spline;
}

void recompute_UBspline_2d_s (UBspline_2d_s* spline, float *data) {
  int Mx = spline->x_grid.num;
  int My = spline->y_grid.num;
  int Nx, Ny, iy, ix;

  if (spline->xBC.lCode == PERIODIC || spline->xBC.lCode == ANTIPERIODIC)
    Nx = Mx+3;
  else
    Nx = Mx+2;
  if (spline->yBC.lCode == PERIODIC || spline->yBC.lCode == ANTIPERIODIC)
    Ny = My+3;
  else
    Ny = My+2;

  /* First, solve in the X-direction */
  for (iy=0; iy<My; iy++) {
    intptr_t doffset = iy;
    intptr_t coffset = iy;
    find_coefs_1d_s (spline->x_grid, spline->xBC, data+doffset, My,
                     spline->coefs+coffset, Ny);
  }

  /* Now, solve in the Y-direction */
  for (ix=0; ix<Nx; ix++) {
    intptr_t doffset = ix*Ny;
    intptr_t coffset = ix*Ny;
    find_coefs_1d_s (spline->y_grid, spline->yBC, spline->coefs+doffset, 1, 
                     spline->coefs+coffset, 1);
  }
}

UBspline_3d_s* create_UBspline_3d_s (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
                                     BCtype_s xBC, BCtype_s yBC, BCtype_s zBC,
                                     float *data) {
  /* Create new spline */
  int Mx, My, Mz;
  int Nx, Ny, Nz, iy, ix, iz;
    
  UBspline_3d_s* restrict spline = malloc (sizeof(UBspline_3d_s));
  spline->spcode = U3D;
  spline->tcode  = SINGLE_REAL;
  spline->xBC = xBC; 
  spline->yBC = yBC; 
  spline->zBC = zBC; 
  /* Setup internal variables */
  Mx = x_grid.num;
  My = y_grid.num;
  Mz = z_grid.num;

  if (xBC.lCode == PERIODIC || xBC.lCode == ANTIPERIODIC)
    Nx = Mx+3;
  else
    Nx = Mx+2;
  x_grid.delta = (x_grid.end - x_grid.start)/(double)(Nx-3);
  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC || yBC.lCode == ANTIPERIODIC)
    Ny = My+3;
  else
    Ny = My+2;
  y_grid.delta = (y_grid.end - y_grid.start)/(double)(Ny-3);
  y_grid.delta_inv = 1.0/y_grid.delta;
  spline->y_grid   = y_grid;

  if (zBC.lCode == PERIODIC || zBC.lCode == ANTIPERIODIC)
    Nz = Mz+3;
  else
    Nz = Mz+2;
  z_grid.delta = (z_grid.end - z_grid.start)/(double)(Nz-3);
  z_grid.delta_inv = 1.0/z_grid.delta;
  spline->z_grid   = z_grid;

  spline->x_stride = Ny*Nz;
  spline->y_stride = Nz;

  spline->coefs = malloc ((size_t)sizeof(float)*(size_t)Nx*(size_t)Ny*(size_t)Nz);
  spline->nc = (size_t)sizeof(float)*(size_t)Nx*(size_t)Ny*(size_t)Nz;

  /* First, solve in the X-direction */
  for (iy=0; iy<My; iy++) 
    for (iz=0; iz<Mz; iz++) {
      intptr_t doffset = iy*Mz+iz;
      intptr_t coffset = iy*Nz+iz;
      find_coefs_1d_s (spline->x_grid, xBC, data+doffset, My*Mz,
                       spline->coefs+coffset, Ny*Nz);
    }

  /* Now, solve in the Y-direction */
  for (ix=0; ix<Nx; ix++) 
    for (iz=0; iz<Nz; iz++) {
      intptr_t doffset = ix*Ny*Nz + iz;
      intptr_t coffset = ix*Ny*Nz + iz;
      find_coefs_1d_s (spline->y_grid, yBC, spline->coefs+doffset, Nz, 
                       spline->coefs+coffset, Nz);
    }

  /* Now, solve in the Z-direction */
  for (ix=0; ix<Nx; ix++) 
    for (iy=0; iy<Ny; iy++) {
      intptr_t doffset = (ix*Ny+iy)*Nz;
      intptr_t coffset = (ix*Ny+iy)*Nz;
      find_coefs_1d_s (spline->z_grid, zBC, spline->coefs+doffset, 1, 
                       spline->coefs+coffset, 1);
    }
  return spline;
}

void recompute_UBspline_3d_s (UBspline_3d_s* spline, float *data) {
  int Mx = spline->x_grid.num;
  int My = spline->y_grid.num;
  int Mz = spline->z_grid.num;
  int Nx, Ny, Nz, ix, iy, iz;

  if (spline->xBC.lCode == PERIODIC || spline->xBC.lCode == ANTIPERIODIC)
    Nx = Mx+3;
  else
    Nx = Mx+2;
  if (spline->yBC.lCode == PERIODIC || spline->yBC.lCode == ANTIPERIODIC)
    Ny = My+3;
  else
    Ny = My+2;
  if (spline->zBC.lCode == PERIODIC || spline->zBC.lCode == ANTIPERIODIC)
    Nz = Mz+3;
  else
    Nz = Mz+2;

  /* First, solve in the X-direction */
  for (iy=0; iy<My; iy++) 
    for (iz=0; iz<Mz; iz++) {
      intptr_t doffset = iy*Mz+iz;
      intptr_t coffset = iy*Nz+iz;
      find_coefs_1d_s (spline->x_grid, spline->xBC, data+doffset, My*Mz,
                       spline->coefs+coffset, Ny*Nz);
    }

  /* Now, solve in the Y-direction */
  for (ix=0; ix<Nx; ix++) 
    for (iz=0; iz<Nz; iz++) {
      intptr_t doffset = ix*Ny*Nz + iz;
      intptr_t coffset = ix*Ny*Nz + iz;
      find_coefs_1d_s (spline->y_grid, spline->yBC, spline->coefs+doffset, Nz, 
                       spline->coefs+coffset, Nz);
    }

  /* Now, solve in the Z-direction */
  for (ix=0; ix<Nx; ix++) 
    for (iy=0; iy<Ny; iy++) {
      intptr_t doffset = (ix*Ny+iy)*Nz;
      intptr_t coffset = (ix*Ny+iy)*Nz;
      find_coefs_1d_s (spline->z_grid, spline->zBC, spline->coefs+doffset, 1, 
                       spline->coefs+coffset, 1);
    }
}

static void destroy_UBspline (Bspline *spline) {
  free (spline->coefs);
  free (spline);
}

/*
 * bspline_eval_std_s.h
 */

/************************************************************/
/* 1D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_1d_s (UBspline_1d_s * restrict spline, 
                         double x, float* restrict val)
{
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;
  
  float tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  float* restrict coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
}

/* Value and first derivative */
void eval_UBspline_1d_s_vg (UBspline_1d_s * restrict spline, double x, 
                            float* restrict val, float* restrict grad)
{
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;
  
  float tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  float* restrict coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
  *grad = spline->x_grid.delta_inv * 
    (coefs[i+0]*(dAf[ 1]*tp[1] + dAf[ 2]*tp[2] + dAf[ 3]*tp[3])+
     coefs[i+1]*(dAf[ 5]*tp[1] + dAf[ 6]*tp[2] + dAf[ 7]*tp[3])+
     coefs[i+2]*(dAf[ 9]*tp[1] + dAf[10]*tp[2] + dAf[11]*tp[3])+
     coefs[i+3]*(dAf[13]*tp[1] + dAf[14]*tp[2] + dAf[15]*tp[3]));
}

/* Value, first derivative, and second derivative */
void eval_UBspline_1d_s_vgl (UBspline_1d_s * restrict spline, double x, 
                             float* restrict val, float* restrict grad,
                             float* restrict lapl)
{
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;
  
  float tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  float* restrict coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
  *grad = spline->x_grid.delta_inv * 
    (coefs[i+0]*(dAf[ 1]*tp[1] + dAf[ 2]*tp[2] + dAf[ 3]*tp[3])+
     coefs[i+1]*(dAf[ 5]*tp[1] + dAf[ 6]*tp[2] + dAf[ 7]*tp[3])+
     coefs[i+2]*(dAf[ 9]*tp[1] + dAf[10]*tp[2] + dAf[11]*tp[3])+
     coefs[i+3]*(dAf[13]*tp[1] + dAf[14]*tp[2] + dAf[15]*tp[3]));
  *lapl = spline->x_grid.delta_inv * spline->x_grid.delta_inv * 
    (coefs[i+0]*(d2Af[ 2]*tp[2] + d2Af[ 3]*tp[3])+
     coefs[i+1]*(d2Af[ 6]*tp[2] + d2Af[ 7]*tp[3])+
     coefs[i+2]*(d2Af[10]*tp[2] + d2Af[11]*tp[3])+
     coefs[i+3]*(d2Af[14]*tp[2] + d2Af[15]*tp[3]));
}

void eval_UBspline_1d_s_vgh (UBspline_1d_s * restrict spline, double x, 
                             float* restrict val, float* restrict grad,
                             float* restrict hess)
{
  eval_UBspline_1d_s_vgl (spline, x, val, grad, hess);
}

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_2d_s (UBspline_2d_s * restrict spline, 
                         double x, double y, float* restrict val)
{
  int xs;
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);
  ty = modff (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;
  
  float tpx[4], tpy[4], a[4], b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0] = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1] = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2] = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3] = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);

  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  
  xs = spline->x_stride;
#define C(i,j) coefs[(ix+(i))*xs+iy+(j)]
  *val = (a[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
          a[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
          a[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
          a[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));
#undef C

}

/* Value and gradient */
void eval_UBspline_2d_s_vg (UBspline_2d_s * restrict spline, 
                            double x, double y, 
                            float* restrict val, float* restrict grad)
{
  int xs;
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);
  ty = modff (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;
  
  float tpx[4], tpy[4], a[4], b[4], da[4], db[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0] = (dAf[ 1]*tpx[1] + dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1] = (dAf[ 5]*tpx[1] + dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2] = (dAf[ 9]*tpx[1] + dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3] = (dAf[13]*tpx[1] + dAf[14]*tpx[2] + dAf[15]*tpx[3]);

  b[0]  = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1]  = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2]  = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3]  = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  
  xs = spline->x_stride;
#define C(i,j) coefs[(ix+(i))*xs+iy+(j)]
  *val =    
    (a[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
     a[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
     a[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
     a[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));
  grad[0] = spline->x_grid.delta_inv *
    (da[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
     da[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
     da[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
     da[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));
  grad[1] = spline->y_grid.delta_inv * 
    (a[0]*(C(0,0)*db[0]+C(0,1)*db[1]+C(0,2)*db[2]+C(0,3)*db[3])+
     a[1]*(C(1,0)*db[0]+C(1,1)*db[1]+C(1,2)*db[2]+C(1,3)*db[3])+
     a[2]*(C(2,0)*db[0]+C(2,1)*db[1]+C(2,2)*db[2]+C(2,3)*db[3])+
     a[3]*(C(3,0)*db[0]+C(3,1)*db[1]+C(3,2)*db[2]+C(3,3)*db[3]));
#undef C

}

/* Value, gradient, and laplacian */
void eval_UBspline_2d_s_vgl (UBspline_2d_s * restrict spline, 
                             double x, double y, float* restrict val, 
                             float* restrict grad, float* restrict lapl)
{
  int xs;
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);
  ty = modff (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;
  
  float tpx[4], tpy[4], a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0]   = (  Af[ 0]*tpx[0] +   Af[ 1]*tpx[1] +  Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]   = (  Af[ 4]*tpx[0] +   Af[ 5]*tpx[1] +  Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]   = (  Af[ 8]*tpx[0] +   Af[ 9]*tpx[1] +  Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]   = (  Af[12]*tpx[0] +   Af[13]*tpx[1] +  Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0]  = ( dAf[ 1]*tpx[1] +  dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1]  = ( dAf[ 5]*tpx[1] +  dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2]  = ( dAf[ 9]*tpx[1] +  dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3]  = ( dAf[13]*tpx[1] +  dAf[14]*tpx[2] + dAf[15]*tpx[3]);
  d2a[0] = (d2Af[ 2]*tpx[2] + d2Af[ 3]*tpx[3]);
  d2a[1] = (d2Af[ 6]*tpx[2] + d2Af[ 7]*tpx[3]);
  d2a[2] = (d2Af[10]*tpx[2] + d2Af[11]*tpx[3]);
  d2a[3] = (d2Af[14]*tpx[2] + d2Af[15]*tpx[3]);

  b[0]  = ( Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1]  = ( Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2]  = ( Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3]  = ( Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  d2b[0] = (d2Af[ 2]*tpy[2] + d2Af[ 3]*tpy[3]);
  d2b[1] = (d2Af[ 6]*tpy[2] + d2Af[ 7]*tpy[3]);
  d2b[2] = (d2Af[10]*tpy[2] + d2Af[11]*tpy[3]);
  d2b[3] = (d2Af[14]*tpy[2] + d2Af[15]*tpy[3]);
  
  xs = spline->x_stride;
#define C(i,j) coefs[(ix+(i))*xs+iy+(j)]
  *val =    
    (a[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
     a[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
     a[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
     a[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));
  grad[0] = spline->x_grid.delta_inv *
    (da[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
     da[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
     da[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
     da[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));
  grad[1] = spline->y_grid.delta_inv *
    (a[0]*(C(0,0)*db[0]+C(0,1)*db[1]+C(0,2)*db[2]+C(0,3)*db[3])+
     a[1]*(C(1,0)*db[0]+C(1,1)*db[1]+C(1,2)*db[2]+C(1,3)*db[3])+
     a[2]*(C(2,0)*db[0]+C(2,1)*db[1]+C(2,2)*db[2]+C(2,3)*db[3])+
     a[3]*(C(3,0)*db[0]+C(3,1)*db[1]+C(3,2)*db[2]+C(3,3)*db[3]));
  *lapl   = 
    spline->y_grid.delta_inv * spline->y_grid.delta_inv *
    (a[0]*(C(0,0)*d2b[0]+C(0,1)*d2b[1]+C(0,2)*d2b[2]+C(0,3)*d2b[3])+
     a[1]*(C(1,0)*d2b[0]+C(1,1)*d2b[1]+C(1,2)*d2b[2]+C(1,3)*d2b[3])+
     a[2]*(C(2,0)*d2b[0]+C(2,1)*d2b[1]+C(2,2)*d2b[2]+C(2,3)*d2b[3])+
     a[3]*(C(3,0)*d2b[0]+C(3,1)*d2b[1]+C(3,2)*d2b[2]+C(3,3)*d2b[3])) + 
    spline->x_grid.delta_inv * spline->x_grid.delta_inv *
     (d2a[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
      d2a[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
      d2a[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
      d2a[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));

#undef C

}

/* Value, gradient, and Hessian */
void eval_UBspline_2d_s_vgh (UBspline_2d_s * restrict spline, 
                             double x, double y, float* restrict val, 
                             float* restrict grad, float* restrict hess)
{
  int xs;
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);
  ty = modff (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;
  
  float tpx[4], tpy[4], a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0]   = (  Af[ 0]*tpx[0] +   Af[ 1]*tpx[1] +  Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]   = (  Af[ 4]*tpx[0] +   Af[ 5]*tpx[1] +  Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]   = (  Af[ 8]*tpx[0] +   Af[ 9]*tpx[1] +  Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]   = (  Af[12]*tpx[0] +   Af[13]*tpx[1] +  Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0]  = ( dAf[ 1]*tpx[1] +  dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1]  = ( dAf[ 5]*tpx[1] +  dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2]  = ( dAf[ 9]*tpx[1] +  dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3]  = ( dAf[13]*tpx[1] +  dAf[14]*tpx[2] + dAf[15]*tpx[3]);
  d2a[0] = (d2Af[ 2]*tpx[2] + d2Af[ 3]*tpx[3]);
  d2a[1] = (d2Af[ 6]*tpx[2] + d2Af[ 7]*tpx[3]);
  d2a[2] = (d2Af[10]*tpx[2] + d2Af[11]*tpx[3]);
  d2a[3] = (d2Af[14]*tpx[2] + d2Af[15]*tpx[3]);

  b[0]   = (  Af[ 0]*tpy[0] +   Af[ 1]*tpy[1] +  Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1]   = (  Af[ 4]*tpy[0] +   Af[ 5]*tpy[1] +  Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2]   = (  Af[ 8]*tpy[0] +   Af[ 9]*tpy[1] +  Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3]   = (  Af[12]*tpy[0] +   Af[13]*tpy[1] +  Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0]  = ( dAf[ 1]*tpy[1] +  dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1]  = ( dAf[ 5]*tpy[1] +  dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2]  = ( dAf[ 9]*tpy[1] +  dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3]  = ( dAf[13]*tpy[1] +  dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  d2b[0] = (d2Af[ 2]*tpy[2] + d2Af[ 3]*tpy[3]);
  d2b[1] = (d2Af[ 6]*tpy[2] + d2Af[ 7]*tpy[3]);
  d2b[2] = (d2Af[10]*tpy[2] + d2Af[11]*tpy[3]);
  d2b[3] = (d2Af[14]*tpy[2] + d2Af[15]*tpy[3]);
  
  xs = spline->x_stride;
#define C(i,j) coefs[(ix+(i))*xs+iy+(j)]
  *val =    
    (  a[0]*(C(0,0)*  b[0]+C(0,1)*  b[1]+C(0,2)*  b[2]+C(0,3)*  b[3])+
       a[1]*(C(1,0)*  b[0]+C(1,1)*  b[1]+C(1,2)*  b[2]+C(1,3)*  b[3])+
       a[2]*(C(2,0)*  b[0]+C(2,1)*  b[1]+C(2,2)*  b[2]+C(2,3)*  b[3])+
       a[3]*(C(3,0)*  b[0]+C(3,1)*  b[1]+C(3,2)*  b[2]+C(3,3)*  b[3]));
  grad[0] = spline->x_grid.delta_inv *
    ( da[0]*(C(0,0)*  b[0]+C(0,1)*  b[1]+C(0,2)*  b[2]+C(0,3)*  b[3])+
      da[1]*(C(1,0)*  b[0]+C(1,1)*  b[1]+C(1,2)*  b[2]+C(1,3)*  b[3])+
      da[2]*(C(2,0)*  b[0]+C(2,1)*  b[1]+C(2,2)*  b[2]+C(2,3)*  b[3])+
      da[3]*(C(3,0)*  b[0]+C(3,1)*  b[1]+C(3,2)*  b[2]+C(3,3)*  b[3]));
  grad[1] = spline->y_grid.delta_inv *
    (  a[0]*(C(0,0)* db[0]+C(0,1)* db[1]+C(0,2)* db[2]+C(0,3)* db[3])+
       a[1]*(C(1,0)* db[0]+C(1,1)* db[1]+C(1,2)* db[2]+C(1,3)* db[3])+
       a[2]*(C(2,0)* db[0]+C(2,1)* db[1]+C(2,2)* db[2]+C(2,3)* db[3])+
       a[3]*(C(3,0)* db[0]+C(3,1)* db[1]+C(3,2)* db[2]+C(3,3)* db[3]));
  hess[0] = spline->x_grid.delta_inv * spline->x_grid.delta_inv *
    (d2a[0]*(C(0,0)*  b[0]+C(0,1)*  b[1]+C(0,2)*  b[2]+C(0,3)*  b[3])+
     d2a[1]*(C(1,0)*  b[0]+C(1,1)*  b[1]+C(1,2)*  b[2]+C(1,3)*  b[3])+
     d2a[2]*(C(2,0)*  b[0]+C(2,1)*  b[1]+C(2,2)*  b[2]+C(2,3)*  b[3])+
     d2a[3]*(C(3,0)*  b[0]+C(3,1)*  b[1]+C(3,2)*  b[2]+C(3,3)*  b[3]));
  hess[1] = spline->x_grid.delta_inv * spline->y_grid.delta_inv *
    ( da[0]*(C(0,0)* db[0]+C(0,1)* db[1]+C(0,2)* db[2]+C(0,3)* db[3])+
      da[1]*(C(1,0)* db[0]+C(1,1)* db[1]+C(1,2)* db[2]+C(1,3)* db[3])+
      da[2]*(C(2,0)* db[0]+C(2,1)* db[1]+C(2,2)* db[2]+C(2,3)* db[3])+
      da[3]*(C(3,0)* db[0]+C(3,1)* db[1]+C(3,2)* db[2]+C(3,3)* db[3]));
  hess[3] = spline->y_grid.delta_inv * spline->y_grid.delta_inv *
    (  a[0]*(C(0,0)*d2b[0]+C(0,1)*d2b[1]+C(0,2)*d2b[2]+C(0,3)*d2b[3])+
       a[1]*(C(1,0)*d2b[0]+C(1,1)*d2b[1]+C(1,2)*d2b[2]+C(1,3)*d2b[3])+
       a[2]*(C(2,0)*d2b[0]+C(2,1)*d2b[1]+C(2,2)*d2b[2]+C(2,3)*d2b[3])+
       a[3]*(C(3,0)*d2b[0]+C(3,1)*d2b[1]+C(3,2)*d2b[2]+C(3,3)*d2b[3]));
  hess[2] = hess[1];

#undef C

}

/************************************************************/
/* 3D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_3d_s (UBspline_3d_s * restrict spline, 
                         double x, double y, double z,
                         float* restrict val)
{
  int xs, ys;
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  ty = modff (uy, &iparty);  int iy = (int) iparty;
  tz = modff (uz, &ipartz);  int iz = (int) ipartz;

  float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0] = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1] = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2] = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3] = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);

  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);

  c[0] = (Af[ 0]*tpz[0] + Af[ 1]*tpz[1] + Af[ 2]*tpz[2] + Af[ 3]*tpz[3]);
  c[1] = (Af[ 4]*tpz[0] + Af[ 5]*tpz[1] + Af[ 6]*tpz[2] + Af[ 7]*tpz[3]);
  c[2] = (Af[ 8]*tpz[0] + Af[ 9]*tpz[1] + Af[10]*tpz[2] + Af[11]*tpz[3]);
  c[3] = (Af[12]*tpz[0] + Af[13]*tpz[1] + Af[14]*tpz[2] + Af[15]*tpz[3]);
  
  xs = spline->x_stride;
  ys = spline->y_stride;
#define P(i,j,k) coefs[(ix+(i))*xs+(iy+(j))*ys+(iz+(k))]
  *val = (a[0]*(b[0]*(P(0,0,0)*c[0]+P(0,0,1)*c[1]+P(0,0,2)*c[2]+P(0,0,3)*c[3])+
                b[1]*(P(0,1,0)*c[0]+P(0,1,1)*c[1]+P(0,1,2)*c[2]+P(0,1,3)*c[3])+
                b[2]*(P(0,2,0)*c[0]+P(0,2,1)*c[1]+P(0,2,2)*c[2]+P(0,2,3)*c[3])+
                b[3]*(P(0,3,0)*c[0]+P(0,3,1)*c[1]+P(0,3,2)*c[2]+P(0,3,3)*c[3]))+
          a[1]*(b[0]*(P(1,0,0)*c[0]+P(1,0,1)*c[1]+P(1,0,2)*c[2]+P(1,0,3)*c[3])+
                b[1]*(P(1,1,0)*c[0]+P(1,1,1)*c[1]+P(1,1,2)*c[2]+P(1,1,3)*c[3])+
                b[2]*(P(1,2,0)*c[0]+P(1,2,1)*c[1]+P(1,2,2)*c[2]+P(1,2,3)*c[3])+
                b[3]*(P(1,3,0)*c[0]+P(1,3,1)*c[1]+P(1,3,2)*c[2]+P(1,3,3)*c[3]))+
          a[2]*(b[0]*(P(2,0,0)*c[0]+P(2,0,1)*c[1]+P(2,0,2)*c[2]+P(2,0,3)*c[3])+
                b[1]*(P(2,1,0)*c[0]+P(2,1,1)*c[1]+P(2,1,2)*c[2]+P(2,1,3)*c[3])+
                b[2]*(P(2,2,0)*c[0]+P(2,2,1)*c[1]+P(2,2,2)*c[2]+P(2,2,3)*c[3])+
                b[3]*(P(2,3,0)*c[0]+P(2,3,1)*c[1]+P(2,3,2)*c[2]+P(2,3,3)*c[3]))+
          a[3]*(b[0]*(P(3,0,0)*c[0]+P(3,0,1)*c[1]+P(3,0,2)*c[2]+P(3,0,3)*c[3])+
                b[1]*(P(3,1,0)*c[0]+P(3,1,1)*c[1]+P(3,1,2)*c[2]+P(3,1,3)*c[3])+
                b[2]*(P(3,2,0)*c[0]+P(3,2,1)*c[1]+P(3,2,2)*c[2]+P(3,2,3)*c[3])+
                b[3]*(P(3,3,0)*c[0]+P(3,3,1)*c[1]+P(3,3,2)*c[2]+P(3,3,3)*c[3])));
#undef P

}

/* Value and gradient */
void eval_UBspline_3d_s_vg (UBspline_3d_s * restrict spline, 
                            double x, double y, double z,
                            float* restrict val, float* restrict grad)
{
  int xs, ys;
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;  
  ty = modff (uy, &iparty);  int iy = (int) iparty; 
  tz = modff (uz, &ipartz);  int iz = (int) ipartz; 
  
  float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], da[4], db[4], dc[4], 
    cP[16], bcP[4], dbcP[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0]   = (  Af[ 0]*tpx[0] +   Af[ 1]*tpx[1] +  Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]   = (  Af[ 4]*tpx[0] +   Af[ 5]*tpx[1] +  Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]   = (  Af[ 8]*tpx[0] +   Af[ 9]*tpx[1] +  Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]   = (  Af[12]*tpx[0] +   Af[13]*tpx[1] +  Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0]  = ( dAf[ 1]*tpx[1] +  dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1]  = ( dAf[ 5]*tpx[1] +  dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2]  = ( dAf[ 9]*tpx[1] +  dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3]  = ( dAf[13]*tpx[1] +  dAf[14]*tpx[2] + dAf[15]*tpx[3]);

  b[0]  = ( Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1]  = ( Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2]  = ( Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3]  = ( Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);

  c[0]  = ( Af[ 0]*tpz[0] + Af[ 1]*tpz[1] + Af[ 2]*tpz[2] + Af[ 3]*tpz[3]);
  c[1]  = ( Af[ 4]*tpz[0] + Af[ 5]*tpz[1] + Af[ 6]*tpz[2] + Af[ 7]*tpz[3]);
  c[2]  = ( Af[ 8]*tpz[0] + Af[ 9]*tpz[1] + Af[10]*tpz[2] + Af[11]*tpz[3]);
  c[3]  = ( Af[12]*tpz[0] + Af[13]*tpz[1] + Af[14]*tpz[2] + Af[15]*tpz[3]);
  dc[0] = (dAf[ 1]*tpz[1] + dAf[ 2]*tpz[2] + dAf[ 3]*tpz[3]);
  dc[1] = (dAf[ 5]*tpz[1] + dAf[ 6]*tpz[2] + dAf[ 7]*tpz[3]);
  dc[2] = (dAf[ 9]*tpz[1] + dAf[10]*tpz[2] + dAf[11]*tpz[3]);
  dc[3] = (dAf[13]*tpz[1] + dAf[14]*tpz[2] + dAf[15]*tpz[3]);
  
  xs = spline->x_stride;
  ys = spline->y_stride;
#define P(i,j,k) coefs[(ix+(i))*xs+(iy+(j))*ys+(iz+(k))]
  cP[ 0] = (P(0,0,0)*c[0]+P(0,0,1)*c[1]+P(0,0,2)*c[2]+P(0,0,3)*c[3]);
  cP[ 1] = (P(0,1,0)*c[0]+P(0,1,1)*c[1]+P(0,1,2)*c[2]+P(0,1,3)*c[3]);
  cP[ 2] = (P(0,2,0)*c[0]+P(0,2,1)*c[1]+P(0,2,2)*c[2]+P(0,2,3)*c[3]);
  cP[ 3] = (P(0,3,0)*c[0]+P(0,3,1)*c[1]+P(0,3,2)*c[2]+P(0,3,3)*c[3]);
  cP[ 4] = (P(1,0,0)*c[0]+P(1,0,1)*c[1]+P(1,0,2)*c[2]+P(1,0,3)*c[3]);
  cP[ 5] = (P(1,1,0)*c[0]+P(1,1,1)*c[1]+P(1,1,2)*c[2]+P(1,1,3)*c[3]);
  cP[ 6] = (P(1,2,0)*c[0]+P(1,2,1)*c[1]+P(1,2,2)*c[2]+P(1,2,3)*c[3]);
  cP[ 7] = (P(1,3,0)*c[0]+P(1,3,1)*c[1]+P(1,3,2)*c[2]+P(1,3,3)*c[3]);
  cP[ 8] = (P(2,0,0)*c[0]+P(2,0,1)*c[1]+P(2,0,2)*c[2]+P(2,0,3)*c[3]);
  cP[ 9] = (P(2,1,0)*c[0]+P(2,1,1)*c[1]+P(2,1,2)*c[2]+P(2,1,3)*c[3]);
  cP[10] = (P(2,2,0)*c[0]+P(2,2,1)*c[1]+P(2,2,2)*c[2]+P(2,2,3)*c[3]);
  cP[11] = (P(2,3,0)*c[0]+P(2,3,1)*c[1]+P(2,3,2)*c[2]+P(2,3,3)*c[3]);
  cP[12] = (P(3,0,0)*c[0]+P(3,0,1)*c[1]+P(3,0,2)*c[2]+P(3,0,3)*c[3]);
  cP[13] = (P(3,1,0)*c[0]+P(3,1,1)*c[1]+P(3,1,2)*c[2]+P(3,1,3)*c[3]);
  cP[14] = (P(3,2,0)*c[0]+P(3,2,1)*c[1]+P(3,2,2)*c[2]+P(3,2,3)*c[3]);
  cP[15] = (P(3,3,0)*c[0]+P(3,3,1)*c[1]+P(3,3,2)*c[2]+P(3,3,3)*c[3]);

  bcP[0] = ( b[0]*cP[ 0] + b[1]*cP[ 1] + b[2]*cP[ 2] + b[3]*cP[ 3]);
  bcP[1] = ( b[0]*cP[ 4] + b[1]*cP[ 5] + b[2]*cP[ 6] + b[3]*cP[ 7]);
  bcP[2] = ( b[0]*cP[ 8] + b[1]*cP[ 9] + b[2]*cP[10] + b[3]*cP[11]);
  bcP[3] = ( b[0]*cP[12] + b[1]*cP[13] + b[2]*cP[14] + b[3]*cP[15]);

  dbcP[0] = ( db[0]*cP[ 0] + db[1]*cP[ 1] + db[2]*cP[ 2] + db[3]*cP[ 3]);
  dbcP[1] = ( db[0]*cP[ 4] + db[1]*cP[ 5] + db[2]*cP[ 6] + db[3]*cP[ 7]);
  dbcP[2] = ( db[0]*cP[ 8] + db[1]*cP[ 9] + db[2]*cP[10] + db[3]*cP[11]);
  dbcP[3] = ( db[0]*cP[12] + db[1]*cP[13] + db[2]*cP[14] + db[3]*cP[15]);

  *val    = ( a[0]*bcP[0] +  a[1]*bcP[1] +  a[2]*bcP[2] +  a[3]*bcP[3]);
  grad[0] = spline->x_grid.delta_inv * 
    (da[0]*bcP[0] + da[1]*bcP[1] + da[2]*bcP[2] + da[3]*bcP[3]);
  grad[1] = spline->y_grid.delta_inv * 
    (a[0]*dbcP[0] + a[1]*dbcP[1] + a[2]*dbcP[2] + a[3]*dbcP[3]);
  grad[2] = spline->z_grid.delta_inv * 
    (a[0]*(b[0]*(P(0,0,0)*dc[0]+P(0,0,1)*dc[1]+P(0,0,2)*dc[2]+P(0,0,3)*dc[3])+
           b[1]*(P(0,1,0)*dc[0]+P(0,1,1)*dc[1]+P(0,1,2)*dc[2]+P(0,1,3)*dc[3])+
           b[2]*(P(0,2,0)*dc[0]+P(0,2,1)*dc[1]+P(0,2,2)*dc[2]+P(0,2,3)*dc[3])+
           b[3]*(P(0,3,0)*dc[0]+P(0,3,1)*dc[1]+P(0,3,2)*dc[2]+P(0,3,3)*dc[3]))+
     a[1]*(b[0]*(P(1,0,0)*dc[0]+P(1,0,1)*dc[1]+P(1,0,2)*dc[2]+P(1,0,3)*dc[3])+
           b[1]*(P(1,1,0)*dc[0]+P(1,1,1)*dc[1]+P(1,1,2)*dc[2]+P(1,1,3)*dc[3])+
           b[2]*(P(1,2,0)*dc[0]+P(1,2,1)*dc[1]+P(1,2,2)*dc[2]+P(1,2,3)*dc[3])+
           b[3]*(P(1,3,0)*dc[0]+P(1,3,1)*dc[1]+P(1,3,2)*dc[2]+P(1,3,3)*dc[3]))+
     a[2]*(b[0]*(P(2,0,0)*dc[0]+P(2,0,1)*dc[1]+P(2,0,2)*dc[2]+P(2,0,3)*dc[3])+
           b[1]*(P(2,1,0)*dc[0]+P(2,1,1)*dc[1]+P(2,1,2)*dc[2]+P(2,1,3)*dc[3])+
           b[2]*(P(2,2,0)*dc[0]+P(2,2,1)*dc[1]+P(2,2,2)*dc[2]+P(2,2,3)*dc[3])+
           b[3]*(P(2,3,0)*dc[0]+P(2,3,1)*dc[1]+P(2,3,2)*dc[2]+P(2,3,3)*dc[3]))+
     a[3]*(b[0]*(P(3,0,0)*dc[0]+P(3,0,1)*dc[1]+P(3,0,2)*dc[2]+P(3,0,3)*dc[3])+
           b[1]*(P(3,1,0)*dc[0]+P(3,1,1)*dc[1]+P(3,1,2)*dc[2]+P(3,1,3)*dc[3])+
           b[2]*(P(3,2,0)*dc[0]+P(3,2,1)*dc[1]+P(3,2,2)*dc[2]+P(3,2,3)*dc[3])+
           b[3]*(P(3,3,0)*dc[0]+P(3,3,1)*dc[1]+P(3,3,2)*dc[2]+P(3,3,3)*dc[3])));
#undef P

}

/* Value, gradient, and laplacian */
void eval_UBspline_3d_s_vgl (UBspline_3d_s * restrict spline, 
                             double x, double y, double z,
                             float* restrict val, float* restrict grad, float* restrict lapl)
{
  int xs, ys;
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;  
  ty = modff (uy, &iparty);  int iy = (int) iparty; 
  tz = modff (uz, &ipartz);  int iz = (int) ipartz; 
  
  float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], da[4], db[4], dc[4], 
    d2a[4], d2b[4], d2c[4], cP[16], dcP[16], bcP[4], dbcP[4], d2bcP[4], bdcP[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0]   = (  Af[ 0]*tpx[0] +   Af[ 1]*tpx[1] +  Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]   = (  Af[ 4]*tpx[0] +   Af[ 5]*tpx[1] +  Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]   = (  Af[ 8]*tpx[0] +   Af[ 9]*tpx[1] +  Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]   = (  Af[12]*tpx[0] +   Af[13]*tpx[1] +  Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0]  = ( dAf[ 1]*tpx[1] +  dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1]  = ( dAf[ 5]*tpx[1] +  dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2]  = ( dAf[ 9]*tpx[1] +  dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3]  = ( dAf[13]*tpx[1] +  dAf[14]*tpx[2] + dAf[15]*tpx[3]);
  d2a[0] = (d2Af[ 2]*tpx[2] + d2Af[ 3]*tpx[3]);
  d2a[1] = (d2Af[ 6]*tpx[2] + d2Af[ 7]*tpx[3]);
  d2a[2] = (d2Af[10]*tpx[2] + d2Af[11]*tpx[3]);
  d2a[3] = (d2Af[14]*tpx[2] + d2Af[15]*tpx[3]);

  b[0]  = ( Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1]  = ( Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2]  = ( Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3]  = ( Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  d2b[0] = (d2Af[ 2]*tpy[2] + d2Af[ 3]*tpy[3]);
  d2b[1] = (d2Af[ 6]*tpy[2] + d2Af[ 7]*tpy[3]);
  d2b[2] = (d2Af[10]*tpy[2] + d2Af[11]*tpy[3]);
  d2b[3] = (d2Af[14]*tpy[2] + d2Af[15]*tpy[3]);

  c[0]  = ( Af[ 0]*tpz[0] + Af[ 1]*tpz[1] + Af[ 2]*tpz[2] + Af[ 3]*tpz[3]);
  c[1]  = ( Af[ 4]*tpz[0] + Af[ 5]*tpz[1] + Af[ 6]*tpz[2] + Af[ 7]*tpz[3]);
  c[2]  = ( Af[ 8]*tpz[0] + Af[ 9]*tpz[1] + Af[10]*tpz[2] + Af[11]*tpz[3]);
  c[3]  = ( Af[12]*tpz[0] + Af[13]*tpz[1] + Af[14]*tpz[2] + Af[15]*tpz[3]);
  dc[0] = (dAf[ 1]*tpz[1] + dAf[ 2]*tpz[2] + dAf[ 3]*tpz[3]);
  dc[1] = (dAf[ 5]*tpz[1] + dAf[ 6]*tpz[2] + dAf[ 7]*tpz[3]);
  dc[2] = (dAf[ 9]*tpz[1] + dAf[10]*tpz[2] + dAf[11]*tpz[3]);
  dc[3] = (dAf[13]*tpz[1] + dAf[14]*tpz[2] + dAf[15]*tpz[3]);
  d2c[0] = (d2Af[ 2]*tpz[2] + d2Af[ 3]*tpz[3]);
  d2c[1] = (d2Af[ 6]*tpz[2] + d2Af[ 7]*tpz[3]);
  d2c[2] = (d2Af[10]*tpz[2] + d2Af[11]*tpz[3]);
  d2c[3] = (d2Af[14]*tpz[2] + d2Af[15]*tpz[3]);
  
  xs = spline->x_stride;
  ys = spline->y_stride;
#define P(i,j,k) coefs[(ix+(i))*xs+(iy+(j))*ys+(iz+(k))]
  cP[ 0] = (P(0,0,0)*c[0]+P(0,0,1)*c[1]+P(0,0,2)*c[2]+P(0,0,3)*c[3]);
  cP[ 1] = (P(0,1,0)*c[0]+P(0,1,1)*c[1]+P(0,1,2)*c[2]+P(0,1,3)*c[3]);
  cP[ 2] = (P(0,2,0)*c[0]+P(0,2,1)*c[1]+P(0,2,2)*c[2]+P(0,2,3)*c[3]);
  cP[ 3] = (P(0,3,0)*c[0]+P(0,3,1)*c[1]+P(0,3,2)*c[2]+P(0,3,3)*c[3]);
  cP[ 4] = (P(1,0,0)*c[0]+P(1,0,1)*c[1]+P(1,0,2)*c[2]+P(1,0,3)*c[3]);
  cP[ 5] = (P(1,1,0)*c[0]+P(1,1,1)*c[1]+P(1,1,2)*c[2]+P(1,1,3)*c[3]);
  cP[ 6] = (P(1,2,0)*c[0]+P(1,2,1)*c[1]+P(1,2,2)*c[2]+P(1,2,3)*c[3]);
  cP[ 7] = (P(1,3,0)*c[0]+P(1,3,1)*c[1]+P(1,3,2)*c[2]+P(1,3,3)*c[3]);
  cP[ 8] = (P(2,0,0)*c[0]+P(2,0,1)*c[1]+P(2,0,2)*c[2]+P(2,0,3)*c[3]);
  cP[ 9] = (P(2,1,0)*c[0]+P(2,1,1)*c[1]+P(2,1,2)*c[2]+P(2,1,3)*c[3]);
  cP[10] = (P(2,2,0)*c[0]+P(2,2,1)*c[1]+P(2,2,2)*c[2]+P(2,2,3)*c[3]);
  cP[11] = (P(2,3,0)*c[0]+P(2,3,1)*c[1]+P(2,3,2)*c[2]+P(2,3,3)*c[3]);
  cP[12] = (P(3,0,0)*c[0]+P(3,0,1)*c[1]+P(3,0,2)*c[2]+P(3,0,3)*c[3]);
  cP[13] = (P(3,1,0)*c[0]+P(3,1,1)*c[1]+P(3,1,2)*c[2]+P(3,1,3)*c[3]);
  cP[14] = (P(3,2,0)*c[0]+P(3,2,1)*c[1]+P(3,2,2)*c[2]+P(3,2,3)*c[3]);
  cP[15] = (P(3,3,0)*c[0]+P(3,3,1)*c[1]+P(3,3,2)*c[2]+P(3,3,3)*c[3]);

  dcP[ 0] = (P(0,0,0)*dc[0]+P(0,0,1)*dc[1]+P(0,0,2)*dc[2]+P(0,0,3)*dc[3]);
  dcP[ 1] = (P(0,1,0)*dc[0]+P(0,1,1)*dc[1]+P(0,1,2)*dc[2]+P(0,1,3)*dc[3]);
  dcP[ 2] = (P(0,2,0)*dc[0]+P(0,2,1)*dc[1]+P(0,2,2)*dc[2]+P(0,2,3)*dc[3]);
  dcP[ 3] = (P(0,3,0)*dc[0]+P(0,3,1)*dc[1]+P(0,3,2)*dc[2]+P(0,3,3)*dc[3]);
  dcP[ 4] = (P(1,0,0)*dc[0]+P(1,0,1)*dc[1]+P(1,0,2)*dc[2]+P(1,0,3)*dc[3]);
  dcP[ 5] = (P(1,1,0)*dc[0]+P(1,1,1)*dc[1]+P(1,1,2)*dc[2]+P(1,1,3)*dc[3]);
  dcP[ 6] = (P(1,2,0)*dc[0]+P(1,2,1)*dc[1]+P(1,2,2)*dc[2]+P(1,2,3)*dc[3]);
  dcP[ 7] = (P(1,3,0)*dc[0]+P(1,3,1)*dc[1]+P(1,3,2)*dc[2]+P(1,3,3)*dc[3]);
  dcP[ 8] = (P(2,0,0)*dc[0]+P(2,0,1)*dc[1]+P(2,0,2)*dc[2]+P(2,0,3)*dc[3]);
  dcP[ 9] = (P(2,1,0)*dc[0]+P(2,1,1)*dc[1]+P(2,1,2)*dc[2]+P(2,1,3)*dc[3]);
  dcP[10] = (P(2,2,0)*dc[0]+P(2,2,1)*dc[1]+P(2,2,2)*dc[2]+P(2,2,3)*dc[3]);
  dcP[11] = (P(2,3,0)*dc[0]+P(2,3,1)*dc[1]+P(2,3,2)*dc[2]+P(2,3,3)*dc[3]);
  dcP[12] = (P(3,0,0)*dc[0]+P(3,0,1)*dc[1]+P(3,0,2)*dc[2]+P(3,0,3)*dc[3]);
  dcP[13] = (P(3,1,0)*dc[0]+P(3,1,1)*dc[1]+P(3,1,2)*dc[2]+P(3,1,3)*dc[3]);
  dcP[14] = (P(3,2,0)*dc[0]+P(3,2,1)*dc[1]+P(3,2,2)*dc[2]+P(3,2,3)*dc[3]);
  dcP[15] = (P(3,3,0)*dc[0]+P(3,3,1)*dc[1]+P(3,3,2)*dc[2]+P(3,3,3)*dc[3]);

  bcP[0] = ( b[0]*cP[ 0] + b[1]*cP[ 1] + b[2]*cP[ 2] + b[3]*cP[ 3]);
  bcP[1] = ( b[0]*cP[ 4] + b[1]*cP[ 5] + b[2]*cP[ 6] + b[3]*cP[ 7]);
  bcP[2] = ( b[0]*cP[ 8] + b[1]*cP[ 9] + b[2]*cP[10] + b[3]*cP[11]);
  bcP[3] = ( b[0]*cP[12] + b[1]*cP[13] + b[2]*cP[14] + b[3]*cP[15]);

  dbcP[0] = ( db[0]*cP[ 0] + db[1]*cP[ 1] + db[2]*cP[ 2] + db[3]*cP[ 3]);
  dbcP[1] = ( db[0]*cP[ 4] + db[1]*cP[ 5] + db[2]*cP[ 6] + db[3]*cP[ 7]);
  dbcP[2] = ( db[0]*cP[ 8] + db[1]*cP[ 9] + db[2]*cP[10] + db[3]*cP[11]);
  dbcP[3] = ( db[0]*cP[12] + db[1]*cP[13] + db[2]*cP[14] + db[3]*cP[15]);

  bdcP[0] = ( b[0]*dcP[ 0] + b[1]*dcP[ 1] + b[2]*dcP[ 2] + b[3]*dcP[ 3]);
  bdcP[1] = ( b[0]*dcP[ 4] + b[1]*dcP[ 5] + b[2]*dcP[ 6] + b[3]*dcP[ 7]);
  bdcP[2] = ( b[0]*dcP[ 8] + b[1]*dcP[ 9] + b[2]*dcP[10] + b[3]*dcP[11]);
  bdcP[3] = ( b[0]*dcP[12] + b[1]*dcP[13] + b[2]*dcP[14] + b[3]*dcP[15]);

  d2bcP[0] = ( d2b[0]*cP[ 0] + d2b[1]*cP[ 1] + d2b[2]*cP[ 2] + d2b[3]*cP[ 3]);
  d2bcP[1] = ( d2b[0]*cP[ 4] + d2b[1]*cP[ 5] + d2b[2]*cP[ 6] + d2b[3]*cP[ 7]);
  d2bcP[2] = ( d2b[0]*cP[ 8] + d2b[1]*cP[ 9] + d2b[2]*cP[10] + d2b[3]*cP[11]);
  d2bcP[3] = ( d2b[0]*cP[12] + d2b[1]*cP[13] + d2b[2]*cP[14] + d2b[3]*cP[15]);


  *val    = 
    ( a[0]*bcP[0] +  a[1]*bcP[1] +  a[2]*bcP[2] +  a[3]*bcP[3]);

  grad[0] = spline->x_grid.delta_inv *
    (da[0]*bcP[0] + da[1]*bcP[1] + da[2]*bcP[2] + da[3]*bcP[3]);
  grad[1] = spline->y_grid.delta_inv * 
    (a[0]*dbcP[0] + a[1]*dbcP[1] + a[2]*dbcP[2] + a[3]*dbcP[3]);
  grad[2] = spline->z_grid.delta_inv * 
    (a[0]*bdcP[0] + a[1]*bdcP[1] + a[2]*bdcP[2] + a[3]*bdcP[3]);

  *lapl = 
    spline->x_grid.delta_inv * spline->x_grid.delta_inv * 
    (d2a[0]*bcP[0] + d2a[1]*bcP[1] + d2a[2]*bcP[2] + d2a[3]*bcP[3])
    
    + spline->y_grid.delta_inv * spline->y_grid.delta_inv * 
    (a[0]*d2bcP[0] + a[1]*d2bcP[1] + a[2]*d2bcP[2] + a[3]*d2bcP[3]) +
    
    + spline->z_grid.delta_inv * spline->z_grid.delta_inv * 
    (a[0]*(b[0]*(P(0,0,0)*d2c[0]+P(0,0,1)*d2c[1]+P(0,0,2)*d2c[2]+P(0,0,3)*d2c[3])+    
           b[1]*(P(0,1,0)*d2c[0]+P(0,1,1)*d2c[1]+P(0,1,2)*d2c[2]+P(0,1,3)*d2c[3])+
           b[2]*(P(0,2,0)*d2c[0]+P(0,2,1)*d2c[1]+P(0,2,2)*d2c[2]+P(0,2,3)*d2c[3])+
           b[3]*(P(0,3,0)*d2c[0]+P(0,3,1)*d2c[1]+P(0,3,2)*d2c[2]+P(0,3,3)*d2c[3]))+
     a[1]*(b[0]*(P(1,0,0)*d2c[0]+P(1,0,1)*d2c[1]+P(1,0,2)*d2c[2]+P(1,0,3)*d2c[3])+
           b[1]*(P(1,1,0)*d2c[0]+P(1,1,1)*d2c[1]+P(1,1,2)*d2c[2]+P(1,1,3)*d2c[3])+
           b[2]*(P(1,2,0)*d2c[0]+P(1,2,1)*d2c[1]+P(1,2,2)*d2c[2]+P(1,2,3)*d2c[3])+
           b[3]*(P(1,3,0)*d2c[0]+P(1,3,1)*d2c[1]+P(1,3,2)*d2c[2]+P(1,3,3)*d2c[3]))+
     a[2]*(b[0]*(P(2,0,0)*d2c[0]+P(2,0,1)*d2c[1]+P(2,0,2)*d2c[2]+P(2,0,3)*d2c[3])+
           b[1]*(P(2,1,0)*d2c[0]+P(2,1,1)*d2c[1]+P(2,1,2)*d2c[2]+P(2,1,3)*d2c[3])+
           b[2]*(P(2,2,0)*d2c[0]+P(2,2,1)*d2c[1]+P(2,2,2)*d2c[2]+P(2,2,3)*d2c[3])+
           b[3]*(P(2,3,0)*d2c[0]+P(2,3,1)*d2c[1]+P(2,3,2)*d2c[2]+P(2,3,3)*d2c[3]))+
     a[3]*(b[0]*(P(3,0,0)*d2c[0]+P(3,0,1)*d2c[1]+P(3,0,2)*d2c[2]+P(3,0,3)*d2c[3])+
           b[1]*(P(3,1,0)*d2c[0]+P(3,1,1)*d2c[1]+P(3,1,2)*d2c[2]+P(3,1,3)*d2c[3])+
           b[2]*(P(3,2,0)*d2c[0]+P(3,2,1)*d2c[1]+P(3,2,2)*d2c[2]+P(3,2,3)*d2c[3])+
           b[3]*(P(3,3,0)*d2c[0]+P(3,3,1)*d2c[1]+P(3,3,2)*d2c[2]+P(3,3,3)*d2c[3])));
#undef P

}

/* Value, gradient, and Hessian */
void eval_UBspline_3d_s_vgh (UBspline_3d_s * restrict spline, 
                             double x, double y, double z,
                             float* restrict val, float* restrict grad, float* restrict hess)
{
  int xs, ys, offmax;
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;  
  z -= spline->z_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
  uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
  uz = fmin (uz, (double)(spline->z_grid.num)-1.0e-5);
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  ty = modff (uy, &iparty);  int iy = (int) iparty;
  tz = modff (uz, &ipartz);  int iz = (int) ipartz;

/*   if ((ix >= spline->x_grid.num))    x = spline->x_grid.num;
     if ((ix < 0))                      x = 0;                 
     if ((iy >= spline->y_grid.num))    y = spline->y_grid.num;
     if ((iy < 0))                      y = 0;                 
     if ((iz >= spline->z_grid.num))    z = spline->z_grid.num;
     if ((iz < 0))                      z = 0;                 */

  float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], da[4], db[4], dc[4], 
    d2a[4], d2b[4], d2c[4], cP[16], dcP[16], d2cP[16], bcP[4], dbcP[4],
    d2bcP[4], dbdcP[4], bd2cP[4], bdcP[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0]   = (  Af[ 0]*tpx[0] +   Af[ 1]*tpx[1] +  Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]   = (  Af[ 4]*tpx[0] +   Af[ 5]*tpx[1] +  Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]   = (  Af[ 8]*tpx[0] +   Af[ 9]*tpx[1] +  Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]   = (  Af[12]*tpx[0] +   Af[13]*tpx[1] +  Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0]  = ( dAf[ 1]*tpx[1] +  dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1]  = ( dAf[ 5]*tpx[1] +  dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2]  = ( dAf[ 9]*tpx[1] +  dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3]  = ( dAf[13]*tpx[1] +  dAf[14]*tpx[2] + dAf[15]*tpx[3]);
  d2a[0] = (d2Af[ 2]*tpx[2] + d2Af[ 3]*tpx[3]);
  d2a[1] = (d2Af[ 6]*tpx[2] + d2Af[ 7]*tpx[3]);
  d2a[2] = (d2Af[10]*tpx[2] + d2Af[11]*tpx[3]);
  d2a[3] = (d2Af[14]*tpx[2] + d2Af[15]*tpx[3]);

  b[0]  = ( Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1]  = ( Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2]  = ( Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3]  = ( Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  d2b[0] = (d2Af[ 2]*tpy[2] + d2Af[ 3]*tpy[3]);
  d2b[1] = (d2Af[ 6]*tpy[2] + d2Af[ 7]*tpy[3]);
  d2b[2] = (d2Af[10]*tpy[2] + d2Af[11]*tpy[3]);
  d2b[3] = (d2Af[14]*tpy[2] + d2Af[15]*tpy[3]);

  c[0]  = ( Af[ 0]*tpz[0] + Af[ 1]*tpz[1] + Af[ 2]*tpz[2] + Af[ 3]*tpz[3]);
  c[1]  = ( Af[ 4]*tpz[0] + Af[ 5]*tpz[1] + Af[ 6]*tpz[2] + Af[ 7]*tpz[3]);
  c[2]  = ( Af[ 8]*tpz[0] + Af[ 9]*tpz[1] + Af[10]*tpz[2] + Af[11]*tpz[3]);
  c[3]  = ( Af[12]*tpz[0] + Af[13]*tpz[1] + Af[14]*tpz[2] + Af[15]*tpz[3]);
  dc[0] = (dAf[ 1]*tpz[1] + dAf[ 2]*tpz[2] + dAf[ 3]*tpz[3]);
  dc[1] = (dAf[ 5]*tpz[1] + dAf[ 6]*tpz[2] + dAf[ 7]*tpz[3]);
  dc[2] = (dAf[ 9]*tpz[1] + dAf[10]*tpz[2] + dAf[11]*tpz[3]);
  dc[3] = (dAf[13]*tpz[1] + dAf[14]*tpz[2] + dAf[15]*tpz[3]);
  d2c[0] = (d2Af[ 2]*tpz[2] + d2Af[ 3]*tpz[3]);
  d2c[1] = (d2Af[ 6]*tpz[2] + d2Af[ 7]*tpz[3]);
  d2c[2] = (d2Af[10]*tpz[2] + d2Af[11]*tpz[3]);
  d2c[3] = (d2Af[14]*tpz[2] + d2Af[15]*tpz[3]);
  
  xs = spline->x_stride;
  ys = spline->y_stride;
  offmax = (ix+3)*xs + (iy+3)*ys + iz+3;
/*   if (offmax > spline->coef_size) {
        fprintf (stderr, "Outside bounds in spline evalutation.\n"
                 "offmax = %d  csize = %d\n", offmax, spline->csize);
        fprintf (stderr, "ix=%d   iy=%d   iz=%d\n", ix,iy,iz);
     }*/
#define P(i,j,k) coefs[(ix+(i))*xs+(iy+(j))*ys+(iz+(k))]
  cP[ 0] = (P(0,0,0)*c[0]+P(0,0,1)*c[1]+P(0,0,2)*c[2]+P(0,0,3)*c[3]);
  cP[ 1] = (P(0,1,0)*c[0]+P(0,1,1)*c[1]+P(0,1,2)*c[2]+P(0,1,3)*c[3]);
  cP[ 2] = (P(0,2,0)*c[0]+P(0,2,1)*c[1]+P(0,2,2)*c[2]+P(0,2,3)*c[3]);
  cP[ 3] = (P(0,3,0)*c[0]+P(0,3,1)*c[1]+P(0,3,2)*c[2]+P(0,3,3)*c[3]);
  cP[ 4] = (P(1,0,0)*c[0]+P(1,0,1)*c[1]+P(1,0,2)*c[2]+P(1,0,3)*c[3]);
  cP[ 5] = (P(1,1,0)*c[0]+P(1,1,1)*c[1]+P(1,1,2)*c[2]+P(1,1,3)*c[3]);
  cP[ 6] = (P(1,2,0)*c[0]+P(1,2,1)*c[1]+P(1,2,2)*c[2]+P(1,2,3)*c[3]);
  cP[ 7] = (P(1,3,0)*c[0]+P(1,3,1)*c[1]+P(1,3,2)*c[2]+P(1,3,3)*c[3]);
  cP[ 8] = (P(2,0,0)*c[0]+P(2,0,1)*c[1]+P(2,0,2)*c[2]+P(2,0,3)*c[3]);
  cP[ 9] = (P(2,1,0)*c[0]+P(2,1,1)*c[1]+P(2,1,2)*c[2]+P(2,1,3)*c[3]);
  cP[10] = (P(2,2,0)*c[0]+P(2,2,1)*c[1]+P(2,2,2)*c[2]+P(2,2,3)*c[3]);
  cP[11] = (P(2,3,0)*c[0]+P(2,3,1)*c[1]+P(2,3,2)*c[2]+P(2,3,3)*c[3]);
  cP[12] = (P(3,0,0)*c[0]+P(3,0,1)*c[1]+P(3,0,2)*c[2]+P(3,0,3)*c[3]);
  cP[13] = (P(3,1,0)*c[0]+P(3,1,1)*c[1]+P(3,1,2)*c[2]+P(3,1,3)*c[3]);
  cP[14] = (P(3,2,0)*c[0]+P(3,2,1)*c[1]+P(3,2,2)*c[2]+P(3,2,3)*c[3]);
  cP[15] = (P(3,3,0)*c[0]+P(3,3,1)*c[1]+P(3,3,2)*c[2]+P(3,3,3)*c[3]);

  dcP[ 0] = (P(0,0,0)*dc[0]+P(0,0,1)*dc[1]+P(0,0,2)*dc[2]+P(0,0,3)*dc[3]);
  dcP[ 1] = (P(0,1,0)*dc[0]+P(0,1,1)*dc[1]+P(0,1,2)*dc[2]+P(0,1,3)*dc[3]);
  dcP[ 2] = (P(0,2,0)*dc[0]+P(0,2,1)*dc[1]+P(0,2,2)*dc[2]+P(0,2,3)*dc[3]);
  dcP[ 3] = (P(0,3,0)*dc[0]+P(0,3,1)*dc[1]+P(0,3,2)*dc[2]+P(0,3,3)*dc[3]);
  dcP[ 4] = (P(1,0,0)*dc[0]+P(1,0,1)*dc[1]+P(1,0,2)*dc[2]+P(1,0,3)*dc[3]);
  dcP[ 5] = (P(1,1,0)*dc[0]+P(1,1,1)*dc[1]+P(1,1,2)*dc[2]+P(1,1,3)*dc[3]);
  dcP[ 6] = (P(1,2,0)*dc[0]+P(1,2,1)*dc[1]+P(1,2,2)*dc[2]+P(1,2,3)*dc[3]);
  dcP[ 7] = (P(1,3,0)*dc[0]+P(1,3,1)*dc[1]+P(1,3,2)*dc[2]+P(1,3,3)*dc[3]);
  dcP[ 8] = (P(2,0,0)*dc[0]+P(2,0,1)*dc[1]+P(2,0,2)*dc[2]+P(2,0,3)*dc[3]);
  dcP[ 9] = (P(2,1,0)*dc[0]+P(2,1,1)*dc[1]+P(2,1,2)*dc[2]+P(2,1,3)*dc[3]);
  dcP[10] = (P(2,2,0)*dc[0]+P(2,2,1)*dc[1]+P(2,2,2)*dc[2]+P(2,2,3)*dc[3]);
  dcP[11] = (P(2,3,0)*dc[0]+P(2,3,1)*dc[1]+P(2,3,2)*dc[2]+P(2,3,3)*dc[3]);
  dcP[12] = (P(3,0,0)*dc[0]+P(3,0,1)*dc[1]+P(3,0,2)*dc[2]+P(3,0,3)*dc[3]);
  dcP[13] = (P(3,1,0)*dc[0]+P(3,1,1)*dc[1]+P(3,1,2)*dc[2]+P(3,1,3)*dc[3]);
  dcP[14] = (P(3,2,0)*dc[0]+P(3,2,1)*dc[1]+P(3,2,2)*dc[2]+P(3,2,3)*dc[3]);
  dcP[15] = (P(3,3,0)*dc[0]+P(3,3,1)*dc[1]+P(3,3,2)*dc[2]+P(3,3,3)*dc[3]);

  d2cP[ 0] = (P(0,0,0)*d2c[0]+P(0,0,1)*d2c[1]+P(0,0,2)*d2c[2]+P(0,0,3)*d2c[3]);
  d2cP[ 1] = (P(0,1,0)*d2c[0]+P(0,1,1)*d2c[1]+P(0,1,2)*d2c[2]+P(0,1,3)*d2c[3]);
  d2cP[ 2] = (P(0,2,0)*d2c[0]+P(0,2,1)*d2c[1]+P(0,2,2)*d2c[2]+P(0,2,3)*d2c[3]);
  d2cP[ 3] = (P(0,3,0)*d2c[0]+P(0,3,1)*d2c[1]+P(0,3,2)*d2c[2]+P(0,3,3)*d2c[3]);
  d2cP[ 4] = (P(1,0,0)*d2c[0]+P(1,0,1)*d2c[1]+P(1,0,2)*d2c[2]+P(1,0,3)*d2c[3]);
  d2cP[ 5] = (P(1,1,0)*d2c[0]+P(1,1,1)*d2c[1]+P(1,1,2)*d2c[2]+P(1,1,3)*d2c[3]);
  d2cP[ 6] = (P(1,2,0)*d2c[0]+P(1,2,1)*d2c[1]+P(1,2,2)*d2c[2]+P(1,2,3)*d2c[3]);
  d2cP[ 7] = (P(1,3,0)*d2c[0]+P(1,3,1)*d2c[1]+P(1,3,2)*d2c[2]+P(1,3,3)*d2c[3]);
  d2cP[ 8] = (P(2,0,0)*d2c[0]+P(2,0,1)*d2c[1]+P(2,0,2)*d2c[2]+P(2,0,3)*d2c[3]);
  d2cP[ 9] = (P(2,1,0)*d2c[0]+P(2,1,1)*d2c[1]+P(2,1,2)*d2c[2]+P(2,1,3)*d2c[3]);
  d2cP[10] = (P(2,2,0)*d2c[0]+P(2,2,1)*d2c[1]+P(2,2,2)*d2c[2]+P(2,2,3)*d2c[3]);
  d2cP[11] = (P(2,3,0)*d2c[0]+P(2,3,1)*d2c[1]+P(2,3,2)*d2c[2]+P(2,3,3)*d2c[3]);
  d2cP[12] = (P(3,0,0)*d2c[0]+P(3,0,1)*d2c[1]+P(3,0,2)*d2c[2]+P(3,0,3)*d2c[3]);
  d2cP[13] = (P(3,1,0)*d2c[0]+P(3,1,1)*d2c[1]+P(3,1,2)*d2c[2]+P(3,1,3)*d2c[3]);
  d2cP[14] = (P(3,2,0)*d2c[0]+P(3,2,1)*d2c[1]+P(3,2,2)*d2c[2]+P(3,2,3)*d2c[3]);
  d2cP[15] = (P(3,3,0)*d2c[0]+P(3,3,1)*d2c[1]+P(3,3,2)*d2c[2]+P(3,3,3)*d2c[3]);

  bcP[0] = ( b[0]*cP[ 0] + b[1]*cP[ 1] + b[2]*cP[ 2] + b[3]*cP[ 3]);
  bcP[1] = ( b[0]*cP[ 4] + b[1]*cP[ 5] + b[2]*cP[ 6] + b[3]*cP[ 7]);
  bcP[2] = ( b[0]*cP[ 8] + b[1]*cP[ 9] + b[2]*cP[10] + b[3]*cP[11]);
  bcP[3] = ( b[0]*cP[12] + b[1]*cP[13] + b[2]*cP[14] + b[3]*cP[15]);

  dbcP[0] = ( db[0]*cP[ 0] + db[1]*cP[ 1] + db[2]*cP[ 2] + db[3]*cP[ 3]);
  dbcP[1] = ( db[0]*cP[ 4] + db[1]*cP[ 5] + db[2]*cP[ 6] + db[3]*cP[ 7]);
  dbcP[2] = ( db[0]*cP[ 8] + db[1]*cP[ 9] + db[2]*cP[10] + db[3]*cP[11]);
  dbcP[3] = ( db[0]*cP[12] + db[1]*cP[13] + db[2]*cP[14] + db[3]*cP[15]);

  bdcP[0] = ( b[0]*dcP[ 0] + b[1]*dcP[ 1] + b[2]*dcP[ 2] + b[3]*dcP[ 3]);
  bdcP[1] = ( b[0]*dcP[ 4] + b[1]*dcP[ 5] + b[2]*dcP[ 6] + b[3]*dcP[ 7]);
  bdcP[2] = ( b[0]*dcP[ 8] + b[1]*dcP[ 9] + b[2]*dcP[10] + b[3]*dcP[11]);
  bdcP[3] = ( b[0]*dcP[12] + b[1]*dcP[13] + b[2]*dcP[14] + b[3]*dcP[15]);

  bd2cP[0] = ( b[0]*d2cP[ 0] + b[1]*d2cP[ 1] + b[2]*d2cP[ 2] + b[3]*d2cP[ 3]);
  bd2cP[1] = ( b[0]*d2cP[ 4] + b[1]*d2cP[ 5] + b[2]*d2cP[ 6] + b[3]*d2cP[ 7]);
  bd2cP[2] = ( b[0]*d2cP[ 8] + b[1]*d2cP[ 9] + b[2]*d2cP[10] + b[3]*d2cP[11]);
  bd2cP[3] = ( b[0]*d2cP[12] + b[1]*d2cP[13] + b[2]*d2cP[14] + b[3]*d2cP[15]);

  d2bcP[0] = ( d2b[0]*cP[ 0] + d2b[1]*cP[ 1] + d2b[2]*cP[ 2] + d2b[3]*cP[ 3]);
  d2bcP[1] = ( d2b[0]*cP[ 4] + d2b[1]*cP[ 5] + d2b[2]*cP[ 6] + d2b[3]*cP[ 7]);
  d2bcP[2] = ( d2b[0]*cP[ 8] + d2b[1]*cP[ 9] + d2b[2]*cP[10] + d2b[3]*cP[11]);
  d2bcP[3] = ( d2b[0]*cP[12] + d2b[1]*cP[13] + d2b[2]*cP[14] + d2b[3]*cP[15]);

  dbdcP[0] = ( db[0]*dcP[ 0] + db[1]*dcP[ 1] + db[2]*dcP[ 2] + db[3]*dcP[ 3]);
  dbdcP[1] = ( db[0]*dcP[ 4] + db[1]*dcP[ 5] + db[2]*dcP[ 6] + db[3]*dcP[ 7]);
  dbdcP[2] = ( db[0]*dcP[ 8] + db[1]*dcP[ 9] + db[2]*dcP[10] + db[3]*dcP[11]);
  dbdcP[3] = ( db[0]*dcP[12] + db[1]*dcP[13] + db[2]*dcP[14] + db[3]*dcP[15]);

  *val = a[0]*bcP[0] + a[1]*bcP[1] + a[2]*bcP[2] + a[3]*bcP[3];
  grad[0] = spline->x_grid.delta_inv *
    (da[0] *bcP[0] + da[1]*bcP[1] + da[2]*bcP[2] + da[3]*bcP[3]);
  grad[1] = spline->y_grid.delta_inv *
    (a[0]*dbcP[0] + a[1]*dbcP[1] + a[2]*dbcP[2] + a[3]*dbcP[3]);
  grad[2] = spline->z_grid.delta_inv *
    (a[0]*bdcP[0] + a[1]*bdcP[1] + a[2]*bdcP[2] + a[3]*bdcP[3]);
  /* d2x */
  hess[0] = spline->x_grid.delta_inv * spline->x_grid.delta_inv *
    (d2a[0]*bcP[0] + d2a[1]*bcP[1] + d2a[2]*bcP[2] + d2a[3]*bcP[3]);
  /* dx dy */
  hess[1] = spline->x_grid.delta_inv * spline->y_grid.delta_inv *
    (da[0]*dbcP[0] + da[1]*dbcP[1] + da[2]*dbcP[2] + da[3]*dbcP[3]);
  hess[3] = hess[1];
  /* dx dz */
  hess[2] = spline->x_grid.delta_inv * spline->z_grid.delta_inv *
    (da[0]*bdcP[0] + da[1]*bdcP[1] + da[2]*bdcP[2] + da[3]*bdcP[3]);
  hess[6] = hess[2];
  /* d2y */
  hess[4] = spline->y_grid.delta_inv * spline->y_grid.delta_inv *
    (a[0]*d2bcP[0] + a[1]*d2bcP[1] + a[2]*d2bcP[2] + a[3]*d2bcP[3]);
  /* dy dz */
  hess[5] = spline->y_grid.delta_inv * spline->z_grid.delta_inv *
    (a[0]*dbdcP[0] + a[1]*dbdcP[1] + a[2]*dbdcP[2] + a[3]*dbdcP[3]);
  hess[7] = hess[5];
  /* d2z */
  hess[8] = spline->z_grid.delta_inv * spline->z_grid.delta_inv *
    (a[0]*bd2cP[0] + a[1]*bd2cP[1] + a[2]*bd2cP[2] + a[3]*bd2cP[3]);
#undef P

}

/*************************************
 *************************************
 **    Uniform, multiple splines    **
 *************************************
 *************************************/


/***********************************************************

        Single-Precision, Real Creation Routines

 ***********************************************************

   On input, bands should be filled with:
   row 0   :  abcdInitial from boundary conditions
   rows 1:M:  basis functions in first 3 cols, data in last
   row M+1 :  abcdFinal   from boundary conditions
   cstride gives the stride between values in coefs.
   On exit, coefs with contain interpolating B-spline coefs */
multi_UBspline_1d_s*
create_multi_UBspline_1d_s (Ugrid x_grid, BCtype_s xBC, int num_splines)
{
  int Mx, Nx, N;
  /* Create new spline */
  multi_UBspline_1d_s* restrict spline = malloc (sizeof(multi_UBspline_1d_s));
  spline->spcode = MULTI_U1D;
  spline->tcode  = SINGLE_REAL;
  spline->xBC = xBC; spline->x_grid = x_grid;
  spline->num_splines = num_splines;

  /* Setup internal variables */
  Mx = x_grid.num;

  if (xBC.lCode == PERIODIC || xBC.lCode == ANTIPERIODIC) {
    x_grid.delta     = (x_grid.end-x_grid.start)/(double)(x_grid.num);
    Nx = Mx+3;
  }
  else {
    x_grid.delta     = (x_grid.end-x_grid.start)/(double)(x_grid.num-1);
    Nx = Mx+2;
  }

  N = num_splines;
  spline->x_stride = N;
  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;
  spline->coefs = malloc ((size_t)sizeof(float)*(size_t)Nx*(size_t)N);
  spline->nc = (size_t)sizeof(float)*(size_t)Nx*(size_t)N;

  return spline;
}

void set_multi_UBspline_1d_s (multi_UBspline_1d_s *spline, int num,
                              float *data)
{
  float *coefs = spline->coefs + num;
  int xs = spline->x_stride;
  find_coefs_1d_s (spline->x_grid, spline->xBC, data, 1, 
                   coefs, xs);
}

multi_UBspline_2d_s*
create_multi_UBspline_2d_s (Ugrid x_grid, Ugrid y_grid,
                            BCtype_s xBC, BCtype_s yBC, int num_splines)
{
  int Mx, My, Nx, Ny, N;
  /* Create new spline */
  multi_UBspline_2d_s* restrict spline = malloc (sizeof(multi_UBspline_2d_s));
  spline->spcode = MULTI_U2D;
  spline->tcode  = SINGLE_REAL;
  spline->xBC = xBC; 
  spline->yBC = yBC; 
  spline->num_splines = num_splines;
  /* Setup internal variables */
  Mx = x_grid.num;
  My = y_grid.num;

  if (xBC.lCode == PERIODIC || xBC.lCode == ANTIPERIODIC)     
    Nx = Mx+3;
  else                           
    Nx = Mx+2;
  x_grid.delta = (x_grid.end - x_grid.start)/(double)(Nx-3);

  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC || yBC.lCode == ANTIPERIODIC)     
    Ny = My+3;
  else                           
    Ny = My+2;
  y_grid.delta = (y_grid.end - y_grid.start)/(double)(Ny-3);
  y_grid.delta_inv = 1.0/y_grid.delta;
  spline->y_grid   = y_grid;

  N = num_splines;
  spline->x_stride = Ny*N;
  spline->y_stride = N;
  spline->coefs = malloc ((size_t)sizeof(float)*(size_t)Nx*(size_t)Ny*(size_t)N);
  spline->nc = (size_t)sizeof(float)*(size_t)Nx*(size_t)Ny*(size_t)N;

  return spline;
}

void set_multi_UBspline_2d_s (multi_UBspline_2d_s* spline, int num, float *data)
{
  int Mx = spline->x_grid.num;
  int My = spline->y_grid.num;
  int Nx, Ny, ys, iy, ix;

  if (spline->xBC.lCode == PERIODIC || spline->xBC.lCode == ANTIPERIODIC)     
    Nx = Mx+3;
  else                                   
    Nx = Mx+2;
  if (spline->yBC.lCode == PERIODIC || spline->yBC.lCode == ANTIPERIODIC)     
    Ny = My+3;
  else                                   
    Ny = My+2;

  float *coefs = spline->coefs + num;
  ys = spline->y_stride;
  /* First, solve in the X-direction */
  for (iy=0; iy<My; iy++) {
    intptr_t doffset = iy;
    intptr_t coffset = iy*ys;
    find_coefs_1d_s (spline->x_grid, spline->xBC, data+doffset, (intptr_t)My,
                     coefs+coffset, (intptr_t)Ny*ys);
  }
  
  /* Now, solve in the Y-direction */
  for (ix=0; ix<Nx; ix++) {
    intptr_t doffset = ix*Ny*ys;
    intptr_t coffset = ix*Ny*ys;
    find_coefs_1d_s (spline->y_grid, spline->yBC, coefs+doffset, (intptr_t)ys, 
                     coefs+coffset, (intptr_t)ys);
  }
}

multi_UBspline_3d_s*
create_multi_UBspline_3d_s (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
                            BCtype_s xBC, BCtype_s yBC, BCtype_s zBC,
                            int num_splines)
{
  int Mx, My, Mz, Nx, Ny, Nz, N;
  /* Create new spline */
  multi_UBspline_3d_s* restrict spline = malloc (sizeof(multi_UBspline_3d_s));
  spline->spcode = MULTI_U3D;
  spline->tcode  = SINGLE_REAL;
  spline->xBC = xBC; 
  spline->yBC = yBC; 
  spline->zBC = zBC; 
  spline->num_splines = num_splines;
  /* Setup internal variables */
  Mx = x_grid.num; My = y_grid.num; Mz = z_grid.num;

  if (xBC.lCode == PERIODIC || xBC.lCode == ANTIPERIODIC)     
    Nx = Mx+3;
  else                           
    Nx = Mx+2;
  x_grid.delta = (x_grid.end - x_grid.start)/(double)(Nx-3);
  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC || yBC.lCode == ANTIPERIODIC)     
    Ny = My+3;
  else                           
    Ny = My+2;
  y_grid.delta = (y_grid.end - y_grid.start)/(double)(Ny-3);
  y_grid.delta_inv = 1.0/y_grid.delta;
  spline->y_grid   = y_grid;

  if (zBC.lCode == PERIODIC || zBC.lCode == ANTIPERIODIC)     
    Nz = Mz+3;
  else                           
    Nz = Mz+2;
  z_grid.delta = (z_grid.end - z_grid.start)/(double)(Nz-3);
  z_grid.delta_inv = 1.0/z_grid.delta;
  spline->z_grid   = z_grid;

  N = num_splines;
  spline->x_stride      = Ny*Nz*N;
  spline->y_stride      = Nz*N;
  spline->z_stride      = N;
  spline->coefs = malloc ((size_t)sizeof(float)*(size_t)Nx*(size_t)Ny*(size_t)Nz*(size_t)N);
  spline->nc = (size_t)sizeof(float)*(size_t)Nx*(size_t)Ny*(size_t)Nz*(size_t)N;

  return spline;
}

void set_multi_UBspline_3d_s (multi_UBspline_3d_s* spline, int num, float *data)
{
  int Mx = spline->x_grid.num;
  int My = spline->y_grid.num;
  int Mz = spline->z_grid.num;
  int Nx, Ny, Nz, zs, ix, iy, iz;

  if (spline->xBC.lCode == PERIODIC || spline->xBC.lCode == ANTIPERIODIC)     
    Nx = Mx+3;
  else                                   
    Nx = Mx+2;
  if (spline->yBC.lCode == PERIODIC || spline->yBC.lCode == ANTIPERIODIC)     
    Ny = My+3;
  else                                   
    Ny = My+2;
  if (spline->zBC.lCode == PERIODIC || spline->zBC.lCode == ANTIPERIODIC)     
    Nz = Mz+3;
  else                                   
    Nz = Mz+2;

  float *coefs = spline->coefs + num;

  zs = spline->z_stride;
  /* First, solve in the X-direction */
  for (iy=0; iy<My; iy++) 
    for (iz=0; iz<Mz; iz++) {
      intptr_t doffset = iy*Mz+iz;
      intptr_t coffset = (iy*Nz+iz)*zs;
      find_coefs_1d_s (spline->x_grid, spline->xBC, data+doffset, (intptr_t)My*Mz,
                       coefs+coffset, (intptr_t)Ny*Nz*zs);
    }
  
  /* Now, solve in the Y-direction */
  for (ix=0; ix<Nx; ix++) 
    for (iz=0; iz<Nz; iz++) {
      intptr_t doffset = (ix*Ny*Nz + iz)*zs;
      intptr_t coffset = (ix*Ny*Nz + iz)*zs;
      find_coefs_1d_s (spline->y_grid, spline->yBC, coefs+doffset, (intptr_t)Nz*zs, 
                       coefs+coffset, (intptr_t)Nz*zs);
    }

  /* Now, solve in the Z-direction */
  for (ix=0; ix<Nx; ix++) 
    for (iy=0; iy<Ny; iy++) {
      intptr_t doffset = ((ix*Ny+iy)*Nz)*zs;
      intptr_t coffset = ((ix*Ny+iy)*Nz)*zs;
      find_coefs_1d_s (spline->z_grid, spline->zBC, coefs+doffset, 
                       (intptr_t)zs, coefs+coffset, (intptr_t)zs);
    }
}

void destroy_multi_UBspline (Bspline *spline)
{
  free (spline->coefs);
  free (spline);
}

/*
 * multi_bspline_eval_std_s.h
 */

/************************************************************/
/* 1D single precision, real evaulation functions           */
/************************************************************/
void eval_multi_UBspline_1d_s (multi_UBspline_1d_s *spline,
                               double x,
                               float* restrict vals)
{
  int n, i;
  x -= spline->x_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float ipartx, tx;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  
  float tpx[4], a[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  float* restrict coefs = spline->coefs;
  
  a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);

  intptr_t xs = spline->x_stride;

  for (n=0; n<spline->num_splines; n++) 
    vals[n]  = 0.0;

  for (i=0; i<4; i++) {
    coefs = spline->coefs + ((ix+i)*xs);
    for (n=0; n<spline->num_splines; n++) 
      vals[n]  +=   a[i] * coefs[n];
  }
}

void eval_multi_UBspline_1d_s_vg (multi_UBspline_1d_s *spline,
                                  double x,
                                  float* restrict vals,
                                  float* restrict grads)
{
  int n, i;
  x -= spline->x_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float ipartx, tx;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  
  float tpx[4], a[4], da[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0] = (dAf[ 0]*tpx[0] + dAf[ 1]*tpx[1] + dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1] = (dAf[ 4]*tpx[0] + dAf[ 5]*tpx[1] + dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2] = (dAf[ 8]*tpx[0] + dAf[ 9]*tpx[1] + dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3] = (dAf[12]*tpx[0] + dAf[13]*tpx[1] + dAf[14]*tpx[2] + dAf[15]*tpx[3]);

  intptr_t xs = spline->x_stride;

  for (n=0; n<spline->num_splines; n++) {
    vals[n]  = 0.0;
    grads[n] = 0.0;
  }

  for (i=0; i<4; i++) { 
    coefs = spline->coefs + ((ix+i)*xs);
    for (n=0; n<spline->num_splines; n++) {
      vals[n]  +=   a[i] * coefs[n];
      grads[n] +=  da[i] * coefs[n];
    }
  }
  
  float dxInv = spline->x_grid.delta_inv;
  for (n=0; n<spline->num_splines; n++) 
    grads[n] *= dxInv;
}

void eval_multi_UBspline_1d_s_vgl (multi_UBspline_1d_s *spline,
                                   double x,
                                   float* restrict vals,
                                   float* restrict grads,
                                   float* restrict lapl)          
{
  int n, i;
  x -= spline->x_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float ipartx, tx;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  
  float tpx[4], a[4], da[4], d2a[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0] = (dAf[ 0]*tpx[0] + dAf[ 1]*tpx[1] + dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1] = (dAf[ 4]*tpx[0] + dAf[ 5]*tpx[1] + dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2] = (dAf[ 8]*tpx[0] + dAf[ 9]*tpx[1] + dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3] = (dAf[12]*tpx[0] + dAf[13]*tpx[1] + dAf[14]*tpx[2] + dAf[15]*tpx[3]);
  d2a[0] = (d2Af[ 0]*tpx[0] + d2Af[ 1]*tpx[1] + d2Af[ 2]*tpx[2] + d2Af[ 3]*tpx[3]);
  d2a[1] = (d2Af[ 4]*tpx[0] + d2Af[ 5]*tpx[1] + d2Af[ 6]*tpx[2] + d2Af[ 7]*tpx[3]);
  d2a[2] = (d2Af[ 8]*tpx[0] + d2Af[ 9]*tpx[1] + d2Af[10]*tpx[2] + d2Af[11]*tpx[3]);
  d2a[3] = (d2Af[12]*tpx[0] + d2Af[13]*tpx[1] + d2Af[14]*tpx[2] + d2Af[15]*tpx[3]);

  intptr_t xs = spline->x_stride;

  for (n=0; n<spline->num_splines; n++) {
    vals[n]  = 0.0;
    grads[n] = 0.0;
    lapl[n]  = 0.0;
  }

  for (i=0; i<4; i++) {      
    coefs = spline->coefs + ((ix+i)*xs);
    for (n=0; n<spline->num_splines; n++) {
      vals[n]  +=   a[i] * coefs[n];
      grads[n] +=  da[i] * coefs[n];
      lapl[n]  += d2a[i] * coefs[n];
    }
  }

  float dxInv = spline->x_grid.delta_inv;
  for (n=0; n<spline->num_splines; n++) {
    grads[n] *= dxInv;
    lapl [n] *= dxInv*dxInv;
  }
}

void eval_multi_UBspline_1d_s_vgh (multi_UBspline_1d_s *spline,
                                   double x,
                                   float* restrict vals,
                                   float* restrict grads,
                                   float* restrict hess)
{
  eval_multi_UBspline_1d_s_vgl (spline, x, vals, grads, hess);
}

/************************************************************/
/* 2D single precision, real evaulation functions           */
/************************************************************/
void eval_multi_UBspline_2d_s (multi_UBspline_2d_s *spline,
                               double x, double y,
                               float* restrict vals)
{
  int n, i, j;
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  ty = modff (uy, &iparty);  int iy = (int) iparty;
  
  float tpx[4], tpy[4], a[4], b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0] = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1] = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2] = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3] = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);

  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;

  for (n=0; n<spline->num_splines; n++)
    vals[n] = 0.0;

  for (i=0; i<4; i++)
    for (j=0; j<4; j++) {
      float prefactor = a[i]*b[j];
      coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
      for (n=0; n<spline->num_splines; n++) 
        vals[n] += prefactor*coefs[n];
    }
}

void eval_multi_UBspline_2d_s_vg (multi_UBspline_2d_s *spline,
                                  double x, double y,
                                  float* restrict vals,
                                  float* restrict grads)
{
  int n, i, j;
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  ty = modff (uy, &iparty);  int iy = (int) iparty;
  
  float tpx[4], tpy[4], a[4], b[4], da[4], db[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0] = (dAf[ 0]*tpx[0] + dAf[ 1]*tpx[1] + dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1] = (dAf[ 4]*tpx[0] + dAf[ 5]*tpx[1] + dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2] = (dAf[ 8]*tpx[0] + dAf[ 9]*tpx[1] + dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3] = (dAf[12]*tpx[0] + dAf[13]*tpx[1] + dAf[14]*tpx[2] + dAf[15]*tpx[3]);

  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 0]*tpy[0] + dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 4]*tpy[0] + dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 8]*tpy[0] + dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[12]*tpy[0] + dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;

  for (n=0; n<spline->num_splines; n++) {
    vals[n] = 0.0;
    grads[2*n+0] = grads[2*n+1] = grads[2*n+2] = 0.0;
  }

  for (i=0; i<4; i++)
    for (j=0; j<4; j++) {
      float ab = a[i]*b[j];
      float dab[2];
      dab[0] = da[i]* b[j];
      dab[1] =  a[i]*db[j];
      
      coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
      for (n=0; n<spline->num_splines; n++) {
        vals [n]     +=   ab   *coefs[n];
        grads[2*n+0] +=  dab[0]*coefs[n];
        grads[2*n+1] +=  dab[1]*coefs[n];
      }
    }

  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  for (n=0; n<spline->num_splines; n++) {
    grads[2*n+0] *= dxInv;
    grads[2*n+1] *= dyInv;
  }
}

void eval_multi_UBspline_2d_s_vgl (multi_UBspline_2d_s *spline,
                                   double x, double y,
                                   float* restrict vals,
                                   float* restrict grads,
                                   float* restrict lapl)          
{
  int n, i, j;
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  ty = modff (uy, &iparty);  int iy = (int) iparty;
  
  float tpx[4], tpy[4], a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0] = (dAf[ 0]*tpx[0] + dAf[ 1]*tpx[1] + dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1] = (dAf[ 4]*tpx[0] + dAf[ 5]*tpx[1] + dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2] = (dAf[ 8]*tpx[0] + dAf[ 9]*tpx[1] + dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3] = (dAf[12]*tpx[0] + dAf[13]*tpx[1] + dAf[14]*tpx[2] + dAf[15]*tpx[3]);
  d2a[0] = (d2Af[ 0]*tpx[0] + d2Af[ 1]*tpx[1] + d2Af[ 2]*tpx[2] + d2Af[ 3]*tpx[3]);
  d2a[1] = (d2Af[ 4]*tpx[0] + d2Af[ 5]*tpx[1] + d2Af[ 6]*tpx[2] + d2Af[ 7]*tpx[3]);
  d2a[2] = (d2Af[ 8]*tpx[0] + d2Af[ 9]*tpx[1] + d2Af[10]*tpx[2] + d2Af[11]*tpx[3]);
  d2a[3] = (d2Af[12]*tpx[0] + d2Af[13]*tpx[1] + d2Af[14]*tpx[2] + d2Af[15]*tpx[3]);

  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 0]*tpy[0] + dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 4]*tpy[0] + dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 8]*tpy[0] + dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[12]*tpy[0] + dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  d2b[0] = (d2Af[ 0]*tpy[0] + d2Af[ 1]*tpy[1] + d2Af[ 2]*tpy[2] + d2Af[ 3]*tpy[3]);
  d2b[1] = (d2Af[ 4]*tpy[0] + d2Af[ 5]*tpy[1] + d2Af[ 6]*tpy[2] + d2Af[ 7]*tpy[3]);
  d2b[2] = (d2Af[ 8]*tpy[0] + d2Af[ 9]*tpy[1] + d2Af[10]*tpy[2] + d2Af[11]*tpy[3]);
  d2b[3] = (d2Af[12]*tpy[0] + d2Af[13]*tpy[1] + d2Af[14]*tpy[2] + d2Af[15]*tpy[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;

  float lapl2[2*spline->num_splines];
  for (n=0; n<spline->num_splines; n++) {
    vals[n] = 0.0;
    grads[2*n+0] = grads[2*n+1] = 0.0;
    lapl2[2*n+0] = lapl2[2*n+1] = 0.0;
  }

  for (i=0; i<4; i++)
    for (j=0; j<4; j++) {
      float ab = a[i]*b[j];
      float dab[2], d2ab[2];
      dab[0] = da[i]* b[j];
      dab[1] =  a[i]*db[j];
      d2ab[0] = d2a[i]*  b[j];
      d2ab[1] =   a[i]*d2b[j];
      
      coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
      for (n=0; n<spline->num_splines; n++) {
        vals[n]      +=   ab   *coefs[n];
        grads[2*n+0] +=  dab[0]*coefs[n];
        grads[2*n+1] +=  dab[1]*coefs[n];
        lapl2[2*n+0] += d2ab[0]*coefs[n];
        lapl2[2*n+1] += d2ab[1]*coefs[n];
      }
    }

  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  for (n=0; n<spline->num_splines; n++) {
    grads[2*n+0] *= dxInv;
    grads[2*n+1] *= dyInv;
    lapl2[2*n+0] *= dxInv*dxInv;
    lapl2[2*n+1] *= dyInv*dyInv;
    lapl[n] = lapl2[2*n+0] + lapl2[2*n+1];
  }
}

void eval_multi_UBspline_2d_s_vgh (multi_UBspline_2d_s *spline,
                                   double x, double y,
                                   float* restrict vals,
                                   float* restrict grads,
                                   float* restrict hess)          
{
  int n, i, j;
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float ipartx, iparty, tx, ty;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  ty = modff (uy, &iparty);  int iy = (int) iparty;
  
  float tpx[4], tpy[4], a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0] = (dAf[ 0]*tpx[0] + dAf[ 1]*tpx[1] + dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1] = (dAf[ 4]*tpx[0] + dAf[ 5]*tpx[1] + dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2] = (dAf[ 8]*tpx[0] + dAf[ 9]*tpx[1] + dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3] = (dAf[12]*tpx[0] + dAf[13]*tpx[1] + dAf[14]*tpx[2] + dAf[15]*tpx[3]);
  d2a[0] = (d2Af[ 0]*tpx[0] + d2Af[ 1]*tpx[1] + d2Af[ 2]*tpx[2] + d2Af[ 3]*tpx[3]);
  d2a[1] = (d2Af[ 4]*tpx[0] + d2Af[ 5]*tpx[1] + d2Af[ 6]*tpx[2] + d2Af[ 7]*tpx[3]);
  d2a[2] = (d2Af[ 8]*tpx[0] + d2Af[ 9]*tpx[1] + d2Af[10]*tpx[2] + d2Af[11]*tpx[3]);
  d2a[3] = (d2Af[12]*tpx[0] + d2Af[13]*tpx[1] + d2Af[14]*tpx[2] + d2Af[15]*tpx[3]);

  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 0]*tpy[0] + dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 4]*tpy[0] + dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 8]*tpy[0] + dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[12]*tpy[0] + dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  d2b[0] = (d2Af[ 0]*tpy[0] + d2Af[ 1]*tpy[1] + d2Af[ 2]*tpy[2] + d2Af[ 3]*tpy[3]);
  d2b[1] = (d2Af[ 4]*tpy[0] + d2Af[ 5]*tpy[1] + d2Af[ 6]*tpy[2] + d2Af[ 7]*tpy[3]);
  d2b[2] = (d2Af[ 8]*tpy[0] + d2Af[ 9]*tpy[1] + d2Af[10]*tpy[2] + d2Af[11]*tpy[3]);
  d2b[3] = (d2Af[12]*tpy[0] + d2Af[13]*tpy[1] + d2Af[14]*tpy[2] + d2Af[15]*tpy[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;

  for (n=0; n<spline->num_splines; n++) {
    vals[n] = 0.0;
    grads[2*n+0] = grads[2*n+1] = 0.0;
    for (i=0; i<4; i++)
      hess[4*n+i] = 0.0;
  }

  for (i=0; i<4; i++)
    for (j=0; j<4; j++){
      float ab = a[i]*b[j];
      float dab[2], d2ab[3];
      dab[0] = da[i]* b[j];
      dab[1] =  a[i]*db[j];
      d2ab[0] = d2a[i] *   b[j];
      d2ab[1] =  da[i] *  db[j];
      d2ab[2] =   a[i] * d2b[j];

      coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
      for (n=0; n<spline->num_splines; n++) {
        vals[n]      +=   ab   *coefs[n];
        grads[2*n+0] +=  dab[0]*coefs[n];
        grads[2*n+1] +=  dab[1]*coefs[n];
        hess [4*n+0] += d2ab[0]*coefs[n];
        hess [4*n+1] += d2ab[1]*coefs[n];
        hess [4*n+3] += d2ab[2]*coefs[n];
      }
    }

  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  for (n=0; n<spline->num_splines; n++) {
    grads[2*n+0] *= dxInv;
    grads[2*n+1] *= dyInv;
    hess[4*n+0] *= dxInv*dxInv;
    hess[4*n+1] *= dxInv*dyInv;
    hess[4*n+3] *= dyInv*dyInv;
    /* Copy hessian elements into lower half of 3x3 matrix */
    hess[4*n+2] = hess[4*n+1];
  }
}

/************************************************************/
/* 3D single precision, real evaulation functions           */
/************************************************************/
void eval_multi_UBspline_3d_s (multi_UBspline_3d_s *spline,
                               double x, double y, double z,
                               float* restrict vals)
{
  int n, i, j, k;
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  ty = modff (uy, &iparty);  int iy = (int) iparty;
  tz = modff (uz, &ipartz);  int iz = (int) ipartz;
  
  float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0] = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1] = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2] = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3] = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);

  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);

  c[0] = (Af[ 0]*tpz[0] + Af[ 1]*tpz[1] + Af[ 2]*tpz[2] + Af[ 3]*tpz[3]);
  c[1] = (Af[ 4]*tpz[0] + Af[ 5]*tpz[1] + Af[ 6]*tpz[2] + Af[ 7]*tpz[3]);
  c[2] = (Af[ 8]*tpz[0] + Af[ 9]*tpz[1] + Af[10]*tpz[2] + Af[11]*tpz[3]);
  c[3] = (Af[12]*tpz[0] + Af[13]*tpz[1] + Af[14]*tpz[2] + Af[15]*tpz[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;

  for (n=0; n<spline->num_splines; n++)
    vals[n] = 0.0;

  for (i=0; i<4; i++)
    for (j=0; j<4; j++) 
      for (k=0; k<4; k++) {
        float abc = a[i]*b[j]*c[k];
        coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
        for (n=0; n<spline->num_splines; n++) 
          vals[n] += abc*coefs[n];
      }
}

void eval_multi_UBspline_3d_s_vg (multi_UBspline_3d_s *spline,
                                  double x, double y, double z,
                                  float* restrict vals,
                                  float* restrict grads)
{
  int n, i, j, k;
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  ty = modff (uy, &iparty);  int iy = (int) iparty;
  tz = modff (uz, &ipartz);  int iz = (int) ipartz;
  
  float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], 
    da[4], db[4], dc[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0] = (dAf[ 0]*tpx[0] + dAf[ 1]*tpx[1] + dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1] = (dAf[ 4]*tpx[0] + dAf[ 5]*tpx[1] + dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2] = (dAf[ 8]*tpx[0] + dAf[ 9]*tpx[1] + dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3] = (dAf[12]*tpx[0] + dAf[13]*tpx[1] + dAf[14]*tpx[2] + dAf[15]*tpx[3]);

  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 0]*tpy[0] + dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 4]*tpy[0] + dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 8]*tpy[0] + dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[12]*tpy[0] + dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);

  c[0] = (Af[ 0]*tpz[0] + Af[ 1]*tpz[1] + Af[ 2]*tpz[2] + Af[ 3]*tpz[3]);
  c[1] = (Af[ 4]*tpz[0] + Af[ 5]*tpz[1] + Af[ 6]*tpz[2] + Af[ 7]*tpz[3]);
  c[2] = (Af[ 8]*tpz[0] + Af[ 9]*tpz[1] + Af[10]*tpz[2] + Af[11]*tpz[3]);
  c[3] = (Af[12]*tpz[0] + Af[13]*tpz[1] + Af[14]*tpz[2] + Af[15]*tpz[3]);
  dc[0] = (dAf[ 0]*tpz[0] + dAf[ 1]*tpz[1] + dAf[ 2]*tpz[2] + dAf[ 3]*tpz[3]);
  dc[1] = (dAf[ 4]*tpz[0] + dAf[ 5]*tpz[1] + dAf[ 6]*tpz[2] + dAf[ 7]*tpz[3]);
  dc[2] = (dAf[ 8]*tpz[0] + dAf[ 9]*tpz[1] + dAf[10]*tpz[2] + dAf[11]*tpz[3]);
  dc[3] = (dAf[12]*tpz[0] + dAf[13]*tpz[1] + dAf[14]*tpz[2] + dAf[15]*tpz[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;

  for (n=0; n<spline->num_splines; n++) {
    vals[n] = 0.0;
    grads[3*n+0] = grads[3*n+1] = grads[3*n+2] = 0.0;
  }

  for (i=0; i<4; i++)
    for (j=0; j<4; j++) 
      for (k=0; k<4; k++) {
        float abc = a[i]*b[j]*c[k];
        float dabc[3];
        dabc[0] = da[i]* b[j]* c[k];
        dabc[1] =  a[i]*db[j]* c[k];
        dabc[2] =  a[i]* b[j]*dc[k];

        coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
        for (n=0; n<spline->num_splines; n++) {
          vals[n]      +=   abc   *coefs[n];
          grads[3*n+0] +=  dabc[0]*coefs[n];
          grads[3*n+1] +=  dabc[1]*coefs[n];
          grads[3*n+2] +=  dabc[2]*coefs[n];
        }
      }

  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv; 
  for (n=0; n<spline->num_splines; n++) {
    grads[3*n+0] *= dxInv;
    grads[3*n+1] *= dyInv;
    grads[3*n+2] *= dzInv;
  }
}

void eval_multi_UBspline_3d_s_vgl (multi_UBspline_3d_s *spline,
                                   double x, double y, double z,
                                   float* restrict vals,
                                   float* restrict grads,
                                   float* restrict lapl)          
{
  int n, i, j, k;
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  ty = modff (uy, &iparty);  int iy = (int) iparty;
  tz = modff (uz, &ipartz);  int iz = (int) ipartz;
  
  float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], 
    da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0] = (dAf[ 0]*tpx[0] + dAf[ 1]*tpx[1] + dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1] = (dAf[ 4]*tpx[0] + dAf[ 5]*tpx[1] + dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2] = (dAf[ 8]*tpx[0] + dAf[ 9]*tpx[1] + dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3] = (dAf[12]*tpx[0] + dAf[13]*tpx[1] + dAf[14]*tpx[2] + dAf[15]*tpx[3]);
  d2a[0] = (d2Af[ 0]*tpx[0] + d2Af[ 1]*tpx[1] + d2Af[ 2]*tpx[2] + d2Af[ 3]*tpx[3]);
  d2a[1] = (d2Af[ 4]*tpx[0] + d2Af[ 5]*tpx[1] + d2Af[ 6]*tpx[2] + d2Af[ 7]*tpx[3]);
  d2a[2] = (d2Af[ 8]*tpx[0] + d2Af[ 9]*tpx[1] + d2Af[10]*tpx[2] + d2Af[11]*tpx[3]);
  d2a[3] = (d2Af[12]*tpx[0] + d2Af[13]*tpx[1] + d2Af[14]*tpx[2] + d2Af[15]*tpx[3]);

  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 0]*tpy[0] + dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 4]*tpy[0] + dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 8]*tpy[0] + dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[12]*tpy[0] + dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  d2b[0] = (d2Af[ 0]*tpy[0] + d2Af[ 1]*tpy[1] + d2Af[ 2]*tpy[2] + d2Af[ 3]*tpy[3]);
  d2b[1] = (d2Af[ 4]*tpy[0] + d2Af[ 5]*tpy[1] + d2Af[ 6]*tpy[2] + d2Af[ 7]*tpy[3]);
  d2b[2] = (d2Af[ 8]*tpy[0] + d2Af[ 9]*tpy[1] + d2Af[10]*tpy[2] + d2Af[11]*tpy[3]);
  d2b[3] = (d2Af[12]*tpy[0] + d2Af[13]*tpy[1] + d2Af[14]*tpy[2] + d2Af[15]*tpy[3]);

  c[0] = (Af[ 0]*tpz[0] + Af[ 1]*tpz[1] + Af[ 2]*tpz[2] + Af[ 3]*tpz[3]);
  c[1] = (Af[ 4]*tpz[0] + Af[ 5]*tpz[1] + Af[ 6]*tpz[2] + Af[ 7]*tpz[3]);
  c[2] = (Af[ 8]*tpz[0] + Af[ 9]*tpz[1] + Af[10]*tpz[2] + Af[11]*tpz[3]);
  c[3] = (Af[12]*tpz[0] + Af[13]*tpz[1] + Af[14]*tpz[2] + Af[15]*tpz[3]);
  dc[0] = (dAf[ 0]*tpz[0] + dAf[ 1]*tpz[1] + dAf[ 2]*tpz[2] + dAf[ 3]*tpz[3]);
  dc[1] = (dAf[ 4]*tpz[0] + dAf[ 5]*tpz[1] + dAf[ 6]*tpz[2] + dAf[ 7]*tpz[3]);
  dc[2] = (dAf[ 8]*tpz[0] + dAf[ 9]*tpz[1] + dAf[10]*tpz[2] + dAf[11]*tpz[3]);
  dc[3] = (dAf[12]*tpz[0] + dAf[13]*tpz[1] + dAf[14]*tpz[2] + dAf[15]*tpz[3]);
  d2c[0] = (d2Af[ 0]*tpz[0] + d2Af[ 1]*tpz[1] + d2Af[ 2]*tpz[2] + d2Af[ 3]*tpz[3]);
  d2c[1] = (d2Af[ 4]*tpz[0] + d2Af[ 5]*tpz[1] + d2Af[ 6]*tpz[2] + d2Af[ 7]*tpz[3]);
  d2c[2] = (d2Af[ 8]*tpz[0] + d2Af[ 9]*tpz[1] + d2Af[10]*tpz[2] + d2Af[11]*tpz[3]);
  d2c[3] = (d2Af[12]*tpz[0] + d2Af[13]*tpz[1] + d2Af[14]*tpz[2] + d2Af[15]*tpz[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;

  float lapl3[3*spline->num_splines];
  for (n=0; n<spline->num_splines; n++) {
    vals[n] = 0.0;
    grads[3*n+0] = grads[3*n+1] = grads[3*n+2] = 0.0;
    lapl3[3*n+0] = lapl3[3*n+1] = lapl3[3*n+2] = 0.0;
  }

  for (i=0; i<4; i++)
    for (j=0; j<4; j++) 
      for (k=0; k<4; k++) {
        float abc = a[i]*b[j]*c[k];
        float dabc[3], d2abc[3];
        dabc[0] = da[i]* b[j]* c[k];
        dabc[1] =  a[i]*db[j]* c[k];
        dabc[2] =  a[i]* b[j]*dc[k];
        d2abc[0] = d2a[i]*  b[j]*  c[k];
        d2abc[1] =   a[i]*d2b[j]*  c[k];
        d2abc[2] =   a[i]*  b[j]*d2c[k];

        coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
        for (n=0; n<spline->num_splines; n++) {
          vals[n]      +=   abc   *coefs[n];
          grads[3*n+0] +=  dabc[0]*coefs[n];
          grads[3*n+1] +=  dabc[1]*coefs[n];
          grads[3*n+2] +=  dabc[2]*coefs[n];
          lapl3[3*n+0] += d2abc[0]*coefs[n];
          lapl3[3*n+1] += d2abc[1]*coefs[n];
          lapl3[3*n+2] += d2abc[2]*coefs[n];
        }
      }

  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv; 
  for (n=0; n<spline->num_splines; n++) {
    grads[3*n+0] *= dxInv;
    grads[3*n+1] *= dyInv;
    grads[3*n+2] *= dzInv;
    lapl3[3*n+0] *= dxInv*dxInv;
    lapl3[3*n+1] *= dyInv*dyInv;
    lapl3[3*n+2] *= dzInv*dzInv;
    lapl[n] = lapl3[3*n+0] + lapl3[3*n+1] + lapl3[3*n+2];
  }
}

void eval_multi_UBspline_3d_s_vgh (multi_UBspline_3d_s *spline,
                                   double x, double y, double z,
                                   float* restrict vals,
                                   float* restrict grads,
                                   float* restrict hess)          
{
  int n, i, j, k;
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  float ux = x*spline->x_grid.delta_inv;
  float uy = y*spline->y_grid.delta_inv;
  float uz = z*spline->z_grid.delta_inv;
  float ipartx, iparty, ipartz, tx, ty, tz;
  tx = modff (ux, &ipartx);  int ix = (int) ipartx;
  ty = modff (uy, &iparty);  int iy = (int) iparty;
  tz = modff (uz, &ipartz);  int iz = (int) ipartz;
  
  float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], 
    da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  float* restrict coefs = spline->coefs;

  a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
  da[0] = (dAf[ 0]*tpx[0] + dAf[ 1]*tpx[1] + dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
  da[1] = (dAf[ 4]*tpx[0] + dAf[ 5]*tpx[1] + dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
  da[2] = (dAf[ 8]*tpx[0] + dAf[ 9]*tpx[1] + dAf[10]*tpx[2] + dAf[11]*tpx[3]);
  da[3] = (dAf[12]*tpx[0] + dAf[13]*tpx[1] + dAf[14]*tpx[2] + dAf[15]*tpx[3]);
  d2a[0] = (d2Af[ 0]*tpx[0] + d2Af[ 1]*tpx[1] + d2Af[ 2]*tpx[2] + d2Af[ 3]*tpx[3]);
  d2a[1] = (d2Af[ 4]*tpx[0] + d2Af[ 5]*tpx[1] + d2Af[ 6]*tpx[2] + d2Af[ 7]*tpx[3]);
  d2a[2] = (d2Af[ 8]*tpx[0] + d2Af[ 9]*tpx[1] + d2Af[10]*tpx[2] + d2Af[11]*tpx[3]);
  d2a[3] = (d2Af[12]*tpx[0] + d2Af[13]*tpx[1] + d2Af[14]*tpx[2] + d2Af[15]*tpx[3]);

  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);
  db[0] = (dAf[ 0]*tpy[0] + dAf[ 1]*tpy[1] + dAf[ 2]*tpy[2] + dAf[ 3]*tpy[3]);
  db[1] = (dAf[ 4]*tpy[0] + dAf[ 5]*tpy[1] + dAf[ 6]*tpy[2] + dAf[ 7]*tpy[3]);
  db[2] = (dAf[ 8]*tpy[0] + dAf[ 9]*tpy[1] + dAf[10]*tpy[2] + dAf[11]*tpy[3]);
  db[3] = (dAf[12]*tpy[0] + dAf[13]*tpy[1] + dAf[14]*tpy[2] + dAf[15]*tpy[3]);
  d2b[0] = (d2Af[ 0]*tpy[0] + d2Af[ 1]*tpy[1] + d2Af[ 2]*tpy[2] + d2Af[ 3]*tpy[3]);
  d2b[1] = (d2Af[ 4]*tpy[0] + d2Af[ 5]*tpy[1] + d2Af[ 6]*tpy[2] + d2Af[ 7]*tpy[3]);
  d2b[2] = (d2Af[ 8]*tpy[0] + d2Af[ 9]*tpy[1] + d2Af[10]*tpy[2] + d2Af[11]*tpy[3]);
  d2b[3] = (d2Af[12]*tpy[0] + d2Af[13]*tpy[1] + d2Af[14]*tpy[2] + d2Af[15]*tpy[3]);

  c[0] = (Af[ 0]*tpz[0] + Af[ 1]*tpz[1] + Af[ 2]*tpz[2] + Af[ 3]*tpz[3]);
  c[1] = (Af[ 4]*tpz[0] + Af[ 5]*tpz[1] + Af[ 6]*tpz[2] + Af[ 7]*tpz[3]);
  c[2] = (Af[ 8]*tpz[0] + Af[ 9]*tpz[1] + Af[10]*tpz[2] + Af[11]*tpz[3]);
  c[3] = (Af[12]*tpz[0] + Af[13]*tpz[1] + Af[14]*tpz[2] + Af[15]*tpz[3]);
  dc[0] = (dAf[ 0]*tpz[0] + dAf[ 1]*tpz[1] + dAf[ 2]*tpz[2] + dAf[ 3]*tpz[3]);
  dc[1] = (dAf[ 4]*tpz[0] + dAf[ 5]*tpz[1] + dAf[ 6]*tpz[2] + dAf[ 7]*tpz[3]);
  dc[2] = (dAf[ 8]*tpz[0] + dAf[ 9]*tpz[1] + dAf[10]*tpz[2] + dAf[11]*tpz[3]);
  dc[3] = (dAf[12]*tpz[0] + dAf[13]*tpz[1] + dAf[14]*tpz[2] + dAf[15]*tpz[3]);
  d2c[0] = (d2Af[ 0]*tpz[0] + d2Af[ 1]*tpz[1] + d2Af[ 2]*tpz[2] + d2Af[ 3]*tpz[3]);
  d2c[1] = (d2Af[ 4]*tpz[0] + d2Af[ 5]*tpz[1] + d2Af[ 6]*tpz[2] + d2Af[ 7]*tpz[3]);
  d2c[2] = (d2Af[ 8]*tpz[0] + d2Af[ 9]*tpz[1] + d2Af[10]*tpz[2] + d2Af[11]*tpz[3]);
  d2c[3] = (d2Af[12]*tpz[0] + d2Af[13]*tpz[1] + d2Af[14]*tpz[2] + d2Af[15]*tpz[3]);

  intptr_t xs = spline->x_stride;
  intptr_t ys = spline->y_stride;
  intptr_t zs = spline->z_stride;

  for (n=0; n<spline->num_splines; n++) {
    vals[n] = 0.0;
    grads[3*n+0] = grads[3*n+1] = grads[3*n+2] = 0.0;
    for (i=0; i<9; i++)
      hess[9*n+i] = 0.0;
  }

  for (i=0; i<4; i++)
    for (j=0; j<4; j++) 
      for (k=0; k<4; k++) {
        float abc = a[i]*b[j]*c[k];
        float dabc[3], d2abc[6];
        dabc[0] = da[i]* b[j]* c[k];
        dabc[1] =  a[i]*db[j]* c[k];
        dabc[2] =  a[i]* b[j]*dc[k];
        d2abc[0] = d2a[i]*  b[j]*  c[k];
        d2abc[1] =  da[i]* db[j]*  c[k];
        d2abc[2] =  da[i]*  b[j]* dc[k];
        d2abc[3] =   a[i]*d2b[j]*  c[k];
        d2abc[4] =   a[i]* db[j]* dc[k];
        d2abc[5] =   a[i]*  b[j]*d2c[k];

        coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
        for (n=0; n<spline->num_splines; n++) {
          vals[n]      +=   abc   *coefs[n];
          grads[3*n+0] +=  dabc[0]*coefs[n];
          grads[3*n+1] +=  dabc[1]*coefs[n];
          grads[3*n+2] +=  dabc[2]*coefs[n];
          hess [9*n+0] += d2abc[0]*coefs[n];
          hess [9*n+1] += d2abc[1]*coefs[n];
          hess [9*n+2] += d2abc[2]*coefs[n];
          hess [9*n+4] += d2abc[3]*coefs[n];
          hess [9*n+5] += d2abc[4]*coefs[n];
          hess [9*n+8] += d2abc[5]*coefs[n];
        }
      }

  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv; 
  for (n=0; n<spline->num_splines; n++) {
    grads[3*n+0] *= dxInv;
    grads[3*n+1] *= dyInv;
    grads[3*n+2] *= dzInv;
    hess [9*n+0] *= dxInv*dxInv;
    hess [9*n+4] *= dyInv*dyInv;
    hess [9*n+8] *= dzInv*dzInv;
    hess [9*n+1] *= dxInv*dyInv;
    hess [9*n+2] *= dxInv*dzInv;
    hess [9*n+5] *= dyInv*dzInv;
    /* Copy hessian elements into lower half of 3x3 matrix */
    hess [9*n+3] = hess[9*n+1];
    hess [9*n+6] = hess[9*n+2];
    hess [9*n+7] = hess[9*n+5];
  }
}


void destroy_Bspline (void *spline) {
  Bspline *sp = (Bspline *)spline;
  if (sp->sp_code <= U3D) 
    destroy_UBspline (sp);
/*  else if (sp->sp_code <= NU3D) 
    destroy_NUBspline (sp);*/
  else if (sp->sp_code <= MULTI_U3D)
    destroy_multi_UBspline (sp);
  else
    fprintf (stderr, "Error in destroy_Bspline:  invalide spline code %d.\n",
             sp->sp_code);
}

