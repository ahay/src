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

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>

#include "einspline.h"
#include "alloc.h"

#ifndef __USE_ISOC99
static double einspl_fmin (double x, double y) {
    if (x < y)
        return x;
    else
        return y;
}
#define fmin einspl_fmin
#endif

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
const float *Af = A44f;

const float dA44f[16] =
  {  0.0, -0.5,  1.0, -0.5,
     0.0,  1.5, -2.0,  0.0,
     0.0, -1.5,  1.0,  0.5,
     0.0,  0.5,  0.0,  0.0 };
const float *dAf = dA44f;

const float d2A44f[16] = 
  {  0.0, 0.0, -1.0,  1.0,
     0.0, 0.0,  3.0, -2.0,
     0.0, 0.0, -3.0,  1.0,
     0.0, 0.0,  1.0,  0.0 };
const float *d2Af = d2A44f;

/********************
 * Double precision *
 ********************/
const double A44d[16] = 
  { -1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
     3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
    -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
     1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0 };
const double *Ad = A44d;

const double dA44d[16] =
  {  0.0, -0.5,  1.0, -0.5,
     0.0,  1.5, -2.0,  0.0,
     0.0, -1.5,  1.0,  0.5,
     0.0,  0.5,  0.0,  0.0 };
const double *dAd = dA44d;

const double d2A44d[16] = 
  {  0.0, 0.0, -1.0,  1.0,
     0.0, 0.0,  3.0, -2.0,
     0.0, 0.0, -3.0,  1.0,
     0.0, 0.0,  1.0,  0.0 };
const double *d2Ad = d2A44d;

/*
 * bspline_create.c
 */

void find_coefs_1d_d (Ugrid grid, BCtype_d bc, 
                      double *data,  intptr_t dstride,
                      double *coefs, intptr_t cstride);

void  solve_deriv_interp_1d_s (float bands[], float coefs[],
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
void  solve_periodic_interp_1d_s (float bands[], float coefs[],
                                  int M, int cstride) {
  float *lastCol;
  int row;

  lastCol = sf_floatalloc(M);

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

void find_coefs_1d_s (Ugrid grid, BCtype_s bc, 
                      float *data,  intptr_t dstride,
                      float *coefs, intptr_t cstride) {
  int M = grid.num, i, j;
  float basis[4] = {1.0/6.0, 2.0/3.0, 1.0/6.0, 0.0};
  float *bands;

  if (bc.lCode == PERIODIC) {
      bands = sf_floatalloc(4*M);
      for (i=0; i<M; i++) {
	  bands[4*i+0] = basis[0];
	  bands[4*i+1] = basis[1];
	  bands[4*i+2] = basis[2];
	  bands[4*i+3] = data[i*dstride];
      }
      solve_periodic_interp_1d_s (bands, coefs, M, cstride);
      free (bands);
  }
  else {
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
    bands = sf_floatalloc ((M+2)*4);
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
  UBspline_1d_s* spline = malloc (sizeof(UBspline_1d_s));
  spline->spcode = U1D;
  spline->tcode  = SINGLE_REAL;
  spline->xBC = xBC; spline->x_grid = x_grid;

  /* Setup internal variables */
  M = x_grid.num;

  if (xBC.lCode == PERIODIC) {
    x_grid.delta     = (x_grid.end-x_grid.start)/(double)(x_grid.num);
    N = M+3;
  }
  else {
    x_grid.delta     = (x_grid.end-x_grid.start)/(double)(x_grid.num-1);
    N = M+2;
  }

  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;
  spline->coefs = malloc (sizeof(float)*N);
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

  UBspline_2d_s* spline = malloc (sizeof(UBspline_2d_s));
  spline->spcode = U2D;
  spline->tcode  = SINGLE_REAL;
  spline->xBC = xBC; 
  spline->yBC = yBC; 
  /* Setup internal variables */
  Mx = x_grid.num;
  My = y_grid.num;


  if (xBC.lCode == PERIODIC)     Nx = Mx+3;
  else                           Nx = Mx+2;
  x_grid.delta = (x_grid.end - x_grid.start)/(double)(Nx-3);
  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC)     Ny = My+3;
  else                           Ny = My+2;
  y_grid.delta = (y_grid.end - y_grid.start)/(double)(Ny-3);
  y_grid.delta_inv = 1.0/y_grid.delta;
  spline->y_grid   = y_grid;
  spline->x_stride = Ny;
  spline->coefs = malloc (sizeof(float)*Nx*Ny);

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

  if (spline->xBC.lCode == PERIODIC)     Nx = Mx+3;
  else                           Nx = Mx+2;
  if (spline->yBC.lCode == PERIODIC)     Ny = My+3;
  else                           Ny = My+2;

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
    
  UBspline_3d_s* spline = malloc (sizeof(UBspline_3d_s));
  spline->spcode = U3D;
  spline->tcode  = SINGLE_REAL;
  spline->xBC = xBC; 
  spline->yBC = yBC; 
  spline->zBC = zBC; 
  /* Setup internal variables */

  if (xBC.lCode == PERIODIC)     Nx = Mx+3;
  else                           Nx = Mx+2;
  x_grid.delta = (x_grid.end - x_grid.start)/(double)(Nx-3);
  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC)     Ny = My+3;
  else                           Ny = My+2;
  y_grid.delta = (y_grid.end - y_grid.start)/(double)(Ny-3);
  y_grid.delta_inv = 1.0/y_grid.delta;
  spline->y_grid   = y_grid;

  if (zBC.lCode == PERIODIC)     Nz = Mz+3;
  else                           Nz = Mz+2;
  z_grid.delta = (z_grid.end - z_grid.start)/(double)(Nz-3);
  z_grid.delta_inv = 1.0/z_grid.delta;
  spline->z_grid   = z_grid;

  spline->x_stride = Ny*Nz;
  spline->y_stride = Nz;

  spline->coefs      = malloc (sizeof(float)*Nx*Ny*Nz);

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

  if (spline->xBC.lCode == PERIODIC)     Nx = Mx+3;
  else                           Nx = Mx+2;
  if (spline->yBC.lCode == PERIODIC)     Ny = My+3;
  else                           Ny = My+2;
  if (spline->zBC.lCode == PERIODIC)     Nz = Mz+3;
  else                           Nz = Mz+2;

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

#ifndef NO_COMPLEX
/***********************************************************

        Single-Precision, Complex Creation Routines

 ***********************************************************

  On input, bands should be filled with:
  row 0   :  abcdInitial from boundary conditions
  rows 1:M:  basis functions in first 3 cols, data in last
  row M+1 :  abcdFinal   from boundary conditions
  cstride gives the stride between values in coefs.
  On exit, coefs with contain interpolating B-spline coefs */

UBspline_1d_c* create_UBspline_1d_c (Ugrid x_grid, BCtype_c xBC, sf_complex *data) {
  /* Create new spline */
  UBspline_1d_c* spline = malloc (sizeof(UBspline_1d_c));
  spline->spcode = U1D;
  spline->tcode  = SINGLE_COMPLEX;
  spline->xBC = xBC; 
  /* Setup internal variables */
  int M = x_grid.num;
  int N;

  if (xBC.lCode == PERIODIC) {
    x_grid.delta     = (x_grid.end-x_grid.start)/(double)(x_grid.num);
    N = M+3;
  }
  else {
    x_grid.delta     = (x_grid.end-x_grid.start)/(double)(x_grid.num-1);
    N = M+2;
  }

  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;
  spline->coefs = malloc (2*sizeof(float)*N);

  BCtype_s xBC_r, xBC_i;
  xBC_r.lCode = xBC.lCode;  xBC_r.rCode = xBC.rCode;
  xBC_r.lVal  = xBC.lVal_r; xBC_r.rVal  = xBC.rVal_r;
  xBC_i.lCode = xBC.lCode;  xBC_i.rCode = xBC.rCode;
  xBC_i.lVal  = xBC.lVal_i; xBC_i.rVal  = xBC.rVal_i;
  /* Real part */
  find_coefs_1d_s (spline->x_grid, xBC_r, 
                   (float*)data, 2, (float*)spline->coefs, 2);
  /* Imaginarty part */
  find_coefs_1d_s (spline->x_grid, xBC_i, 
                   ((float*)data)+1, 2, ((float*)spline->coefs+1), 2);

  return spline;
}

void recompute_UBspline_1d_c (UBspline_1d_c* spline, sf_complex *data) {

  BCtype_s xBC_r, xBC_i;
  xBC_r.lCode = spline->xBC.lCode;  xBC_r.rCode = spline->xBC.rCode;
  xBC_r.lVal  = spline->xBC.lVal_r; xBC_r.rVal  = spline->xBC.rVal_r;
  xBC_i.lCode = spline->xBC.lCode;  xBC_i.rCode = spline->xBC.rCode;
  xBC_i.lVal  = spline->xBC.lVal_i; xBC_i.rVal  = spline->xBC.rVal_i;

  /* Real part */
  find_coefs_1d_s (spline->x_grid, xBC_r, 
                   (float*)data, 2, (float*)spline->coefs, 2);
  /* Imaginarty part */
  find_coefs_1d_s (spline->x_grid, xBC_i, 
                   ((float*)data)+1, 2, ((float*)spline->coefs+1), 2);
}

UBspline_2d_c* create_UBspline_2d_c (Ugrid x_grid, Ugrid y_grid,
                                     BCtype_c xBC, BCtype_c yBC, sf_complex *data) {
  /* Create new spline */
  UBspline_2d_c* spline = malloc (sizeof(UBspline_2d_c));
  spline->spcode = U2D;
  spline->tcode  = SINGLE_COMPLEX;
  spline->xBC = xBC; 
  spline->yBC = yBC; 

  /* Setup internal variables */
  int Mx = x_grid.num;
  int My = y_grid.num;
  int Nx, Ny, ix, iy;

  if (xBC.lCode == PERIODIC)     Nx = Mx+3;
  else                           Nx = Mx+2;
  x_grid.delta = (x_grid.end - x_grid.start)/(double)(Nx-3);
  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC)     Ny = My+3;
  else                           Ny = My+2;
  y_grid.delta = (y_grid.end - y_grid.start)/(double)(Ny-3);
  y_grid.delta_inv = 1.0/y_grid.delta;
  spline->y_grid   = y_grid;
  spline->x_stride = Ny;

  spline->coefs = malloc (2*sizeof(float)*Nx*Ny);

  BCtype_s xBC_r, xBC_i, yBC_r, yBC_i;
  xBC_r.lCode = xBC.lCode;  xBC_r.rCode = xBC.rCode;
  xBC_r.lVal  = xBC.lVal_r; xBC_r.rVal  = xBC.rVal_r;
  xBC_i.lCode = xBC.lCode;  xBC_i.rCode = xBC.rCode;
  xBC_i.lVal  = xBC.lVal_i; xBC_i.rVal  = xBC.rVal_i;
  yBC_r.lCode = yBC.lCode;  yBC_r.rCode = yBC.rCode;
  yBC_r.lVal  = yBC.lVal_r; yBC_r.rVal  = yBC.rVal_r;
  yBC_i.lCode = yBC.lCode;  yBC_i.rCode = yBC.rCode;
  yBC_i.lVal  = yBC.lVal_i; yBC_i.rVal  = yBC.rVal_i;
  /* First, solve in the X-direction */
  for (iy=0; iy<My; iy++) {
    intptr_t doffset = 2*iy;
    intptr_t coffset = 2*iy;
    /* Real part */
    find_coefs_1d_s (spline->x_grid, xBC_r, ((float*)data)+doffset, 2*My,
                     (float*)spline->coefs+coffset, 2*Ny);
    /* Imag part */
    find_coefs_1d_s (spline->x_grid, xBC_i, ((float*)data)+doffset+1, 2*My,
                     ((float*)spline->coefs)+coffset+1, 2*Ny);
  }

  /* Now, solve in the Y-direction */
  for (ix=0; ix<Nx; ix++) {
    intptr_t doffset = 2*ix*Ny;
    intptr_t coffset = 2*ix*Ny;
    /* Real part */
    find_coefs_1d_s (spline->y_grid, yBC_r, ((float*)spline->coefs)+doffset, 2, 
                     ((float*)spline->coefs)+coffset, 2);
    /* Imag part */
    find_coefs_1d_s (spline->y_grid, yBC_i, ((float*)spline->coefs)+doffset+1, 2, 
                     ((float*)spline->coefs)+coffset+1, 2);
  }
  return spline;
}

void recompute_UBspline_2d_c (UBspline_2d_c* spline, sf_complex *data) {
  /* Setup internal variables */
  int Mx = spline->x_grid.num;
  int My = spline->y_grid.num;
  int Nx, Ny, ix, iy;

  if (spline->xBC.lCode == PERIODIC)     Nx = Mx+3;
  else                           Nx = Mx+2;
  if (spline->yBC.lCode == PERIODIC)     Ny = My+3;
  else                           Ny = My+2;

  BCtype_s xBC_r, xBC_i, yBC_r, yBC_i;
  xBC_r.lCode = spline->xBC.lCode;  xBC_r.rCode = spline->xBC.rCode;
  xBC_r.lVal  = spline->xBC.lVal_r; xBC_r.rVal  = spline->xBC.rVal_r;
  xBC_i.lCode = spline->xBC.lCode;  xBC_i.rCode = spline->xBC.rCode;
  xBC_i.lVal  = spline->xBC.lVal_i; xBC_i.rVal  = spline->xBC.rVal_i;
  yBC_r.lCode = spline->yBC.lCode;  yBC_r.rCode = spline->yBC.rCode;
  yBC_r.lVal  = spline->yBC.lVal_r; yBC_r.rVal  = spline->yBC.rVal_r;
  yBC_i.lCode = spline->yBC.lCode;  yBC_i.rCode = spline->yBC.rCode;
  yBC_i.lVal  = spline->yBC.lVal_i; yBC_i.rVal  = spline->yBC.rVal_i;

  /* First, solve in the X-direction */
  for (iy=0; iy<My; iy++) {
    intptr_t doffset = 2*iy;
    intptr_t coffset = 2*iy;
    /* Real part */
    find_coefs_1d_s (spline->x_grid, xBC_r, ((float*)data)+doffset, 2*My,
                     (float*)spline->coefs+coffset, 2*Ny);
    /* Imag part */
    find_coefs_1d_s (spline->x_grid, xBC_i, ((float*)data)+doffset+1, 2*My,
                     ((float*)spline->coefs)+coffset+1, 2*Ny);
  }

  /* Now, solve in the Y-direction */
  for (ix=0; ix<Nx; ix++) {
    intptr_t doffset = 2*ix*Ny;
    intptr_t coffset = 2*ix*Ny;
    /* Real part */
    find_coefs_1d_s (spline->y_grid, yBC_r, ((float*)spline->coefs)+doffset, 2, 
                     ((float*)spline->coefs)+coffset, 2);
    /* Imag part */
    find_coefs_1d_s (spline->y_grid, yBC_i, ((float*)spline->coefs)+doffset+1, 2, 
                     ((float*)spline->coefs)+coffset+1, 2);
  }  
}

UBspline_3d_c* create_UBspline_3d_c (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
                                     BCtype_c xBC, BCtype_c yBC, BCtype_c zBC,
                                     sf_complex *data) {
  /* Create new spline */
  UBspline_3d_c* spline = malloc (sizeof(UBspline_3d_c));
  spline->spcode = U3D;
  spline->tcode  = SINGLE_COMPLEX;
  spline->xBC = xBC; 
  spline->yBC = yBC; 
  spline->zBC = zBC; 

  /* Setup internal variables */
  int Mx = x_grid.num;  int My = y_grid.num; int Mz = z_grid.num;
  int Nx, Ny, Nz, ix, iy, iz;

  if (xBC.lCode == PERIODIC)     Nx = Mx+3;
  else                           Nx = Mx+2;
  x_grid.delta = (x_grid.end - x_grid.start)/(double)(Nx-3);
  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC)     Ny = My+3;
  else                           Ny = My+2;
  y_grid.delta = (y_grid.end - y_grid.start)/(double)(Ny-3);
  y_grid.delta_inv = 1.0/y_grid.delta;
  spline->y_grid   = y_grid;

  if (zBC.lCode == PERIODIC)     Nz = Mz+3;
  else                           Nz = Mz+2;
  z_grid.delta = (z_grid.end - z_grid.start)/(double)(Nz-3);
  z_grid.delta_inv = 1.0/z_grid.delta;
  spline->z_grid   = z_grid;

  spline->x_stride = Ny*Nz;
  spline->y_stride = Nz;

  spline->coefs      = malloc (2*sizeof(float)*Nx*Ny*Nz);

  BCtype_s xBC_r, xBC_i, yBC_r, yBC_i, zBC_r, zBC_i;
  xBC_r.lCode = xBC.lCode;  xBC_r.rCode = xBC.rCode;
  xBC_r.lVal  = xBC.lVal_r; xBC_r.rVal  = xBC.rVal_r;
  xBC_i.lCode = xBC.lCode;  xBC_i.rCode = xBC.rCode;
  xBC_i.lVal  = xBC.lVal_i; xBC_i.rVal  = xBC.rVal_i;
  yBC_r.lCode = yBC.lCode;  yBC_r.rCode = yBC.rCode;
  yBC_r.lVal  = yBC.lVal_r; yBC_r.rVal  = yBC.rVal_r;
  yBC_i.lCode = yBC.lCode;  yBC_i.rCode = yBC.rCode;
  yBC_i.lVal  = yBC.lVal_i; yBC_i.rVal  = yBC.rVal_i;
  zBC_r.lCode = zBC.lCode;  zBC_r.rCode = zBC.rCode;
  zBC_r.lVal  = zBC.lVal_r; zBC_r.rVal  = zBC.rVal_r;
  zBC_i.lCode = zBC.lCode;  zBC_i.rCode = zBC.rCode;
  zBC_i.lVal  = zBC.lVal_i; zBC_i.rVal  = zBC.rVal_i;
  /* First, solve in the X-direction */
  for (iy=0; iy<My; iy++) 
    for (iz=0; iz<Mz; iz++) {
      intptr_t doffset = 2*(iy*Mz+iz);
      intptr_t coffset = 2*(iy*Nz+iz);
      /* Real part */
      find_coefs_1d_s (spline->x_grid, xBC_r, ((float*)data)+doffset, 2*My*Mz,
                       ((float*)spline->coefs)+coffset, 2*Ny*Nz);
      /* Imag part */
      find_coefs_1d_s (spline->x_grid, xBC_i, ((float*)data)+doffset+1, 2*My*Mz,
                       ((float*)spline->coefs)+coffset+1, 2*Ny*Nz);
    }

  /* Now, solve in the Y-direction */
  for (ix=0; ix<Nx; ix++) 
    for (iz=0; iz<Nz; iz++) {
      intptr_t doffset = 2*(ix*Ny*Nz + iz);
      intptr_t coffset = 2*(ix*Ny*Nz + iz);
      /* Real part */
      find_coefs_1d_s (spline->y_grid, yBC_r, ((float*)spline->coefs)+doffset, 2*Nz, 
                       ((float*)spline->coefs)+coffset, 2*Nz);
      /* Imag part */
      find_coefs_1d_s (spline->y_grid, yBC_i, ((float*)spline->coefs)+doffset+1, 2*Nz, 
                       ((float*)spline->coefs)+coffset+1, 2*Nz);
    }

  /* Now, solve in the Z-direction */
  for (ix=0; ix<Nx; ix++) 
    for (iy=0; iy<Ny; iy++) {
      intptr_t doffset = 2*((ix*Ny+iy)*Nz);
      intptr_t coffset = 2*((ix*Ny+iy)*Nz);
      /* Real part */
      find_coefs_1d_s (spline->z_grid, zBC_r, ((float*)spline->coefs)+doffset, 2, 
                       ((float*)spline->coefs)+coffset, 2);
      /* Imag part */
      find_coefs_1d_s (spline->z_grid, zBC_i, ((float*)spline->coefs)+doffset+1, 2, 
                       ((float*)spline->coefs)+coffset+1, 2);
    }

  return spline;
}

void recompute_UBspline_3d_c (UBspline_3d_c* spline, sf_complex *data) {
  /* Setup internal variables */
  int Mx = spline->x_grid.num;
  int My = spline->y_grid.num;
  int Mz = spline->z_grid.num;
  int Nx, Ny, Nz, ix, iy, iz;

  if (spline->xBC.lCode == PERIODIC)     Nx = Mx+3;
  else                           Nx = Mx+2;
  if (spline->yBC.lCode == PERIODIC)     Ny = My+3;
  else                           Ny = My+2;
  if (spline->zBC.lCode == PERIODIC)     Nz = Mz+3;
  else                           Nz = Mz+2;

  BCtype_s xBC_r, xBC_i, yBC_r, yBC_i, zBC_r, zBC_i;
  xBC_r.lCode = spline->xBC.lCode;  xBC_r.rCode = spline->xBC.rCode;
  xBC_r.lVal  = spline->xBC.lVal_r; xBC_r.rVal  = spline->xBC.rVal_r;
  xBC_i.lCode = spline->xBC.lCode;  xBC_i.rCode = spline->xBC.rCode;
  xBC_i.lVal  = spline->xBC.lVal_i; xBC_i.rVal  = spline->xBC.rVal_i;
  yBC_r.lCode = spline->yBC.lCode;  yBC_r.rCode = spline->yBC.rCode;
  yBC_r.lVal  = spline->yBC.lVal_r; yBC_r.rVal  = spline->yBC.rVal_r;
  yBC_i.lCode = spline->yBC.lCode;  yBC_i.rCode = spline->yBC.rCode;
  yBC_i.lVal  = spline->yBC.lVal_i; yBC_i.rVal  = spline->yBC.rVal_i;
  zBC_r.lCode = spline->zBC.lCode;  zBC_r.rCode = spline->zBC.rCode;
  zBC_r.lVal  = spline->zBC.lVal_r; zBC_r.rVal  = spline->zBC.rVal_r;
  zBC_i.lCode = spline->zBC.lCode;  zBC_i.rCode = spline->zBC.rCode;
  zBC_i.lVal  = spline->zBC.lVal_i; zBC_i.rVal  = spline->zBC.rVal_i;
  /* First, solve in the X-direction */
  for (iy=0; iy<My; iy++) 
    for (iz=0; iz<Mz; iz++) {
      intptr_t doffset = 2*(iy*Mz+iz);
      intptr_t coffset = 2*(iy*Nz+iz);
      /* Real part */
      find_coefs_1d_s (spline->x_grid, xBC_r, ((float*)data)+doffset, 2*My*Mz,
                       ((float*)spline->coefs)+coffset, 2*Ny*Nz);
      /* Imag part */
      find_coefs_1d_s (spline->x_grid, xBC_i, ((float*)data)+doffset+1, 2*My*Mz,
                       ((float*)spline->coefs)+coffset+1, 2*Ny*Nz);
    }

  /* Now, solve in the Y-direction */
  for (ix=0; ix<Nx; ix++) 
    for (iz=0; iz<Nz; iz++) {
      intptr_t doffset = 2*(ix*Ny*Nz + iz);
      intptr_t coffset = 2*(ix*Ny*Nz + iz);
      /* Real part */
      find_coefs_1d_s (spline->y_grid, yBC_r, ((float*)spline->coefs)+doffset, 2*Nz, 
                       ((float*)spline->coefs)+coffset, 2*Nz);
      /* Imag part */
      find_coefs_1d_s (spline->y_grid, yBC_i, ((float*)spline->coefs)+doffset+1, 2*Nz, 
                       ((float*)spline->coefs)+coffset+1, 2*Nz);
    }

  /* Now, solve in the Z-direction */
  for (ix=0; ix<Nx; ix++) 
    for (iy=0; iy<Ny; iy++) {
      intptr_t doffset = 2*((ix*Ny+iy)*Nz);
      intptr_t coffset = 2*((ix*Ny+iy)*Nz);
      /* Real part */
      find_coefs_1d_s (spline->z_grid, zBC_r, ((float*)spline->coefs)+doffset, 2, 
                       ((float*)spline->coefs)+coffset, 2);
      /* Imag part */
      find_coefs_1d_s (spline->z_grid, zBC_i, ((float*)spline->coefs)+doffset+1, 2, 
                       ((float*)spline->coefs)+coffset+1, 2);
    }
}
#endif /* NO_COMPLEX */

/***********************************************************

        Double-Precision, Real Creation Routines

 ***********************************************************

  On input, bands should be filled with:
  row 0   :  abcdInitial from boundary conditions
  rows 1:M:  basis functions in first 3 cols, data in last
  row M+1 :  abcdFinal   from boundary conditions
  cstride gives the stride between values in coefs.
  On exit, coefs with contain interpolating B-spline coefs */

void solve_deriv_interp_1d_d (double bands[], double coefs[],
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
void  solve_periodic_interp_1d_d (double bands[], double coefs[],
                                  int M, int cstride)
{
  double *lastCol;
  int row;

  lastCol = (double*) sf_alloc(M,sizeof(double));

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
     The [2] element and [0] element are now on top of each other  */
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

void find_coefs_1d_d (Ugrid grid, BCtype_d bc, 
                      double *data,  intptr_t dstride,
                      double *coefs, intptr_t cstride) {
  int M = grid.num, i, j;
  double basis[4] = {1.0/6.0, 2.0/3.0, 1.0/6.0, 0.0};
  double *bands;

  if (bc.lCode == PERIODIC) {
    bands = (double*) sf_alloc (4*M,sizeof(double));
    for (i=0; i<M; i++) {
      bands[4*i+0] = basis[0];
      bands[4*i+1] = basis[1];
      bands[4*i+2] = basis[2];
      bands[4*i+3] = data[i*dstride];
    }
    solve_periodic_interp_1d_d (bands, coefs, M, cstride);
    free (bands);
  }
  else {
    /* Setup boundary conditions */
    double abcd_left[4], abcd_right[4];
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
    bands = (double*) sf_alloc ((M+2)*4,sizeof(double));
    for (i=0; i<4; i++) {
      bands[4*( 0 )+i] = abcd_left[i];
      bands[4*(M+1)+i] = abcd_right[i];
    }
    for (i=0; i<M; i++) {
      for (j=0; j<3; j++)
        bands[4*(i+1)+j] = basis[j];
      bands[4*(i+1)+3] = data[i*dstride];
    }   
    /* Now, solve for coefficients */
    solve_deriv_interp_1d_d (bands, coefs, M, cstride);
    free (bands);
  }
}

UBspline_1d_d* create_UBspline_1d_d (Ugrid x_grid, BCtype_d xBC, double *data) {
  /* Create new spline */
    int M, N;

  UBspline_1d_d* spline = malloc (sizeof(UBspline_1d_d));
  spline->spcode = U1D;
  spline->tcode  = DOUBLE_REAL;
  spline->xBC = xBC; 

  /* Setup internal variables */
  M = x_grid.num;

  if (xBC.lCode == PERIODIC) {
    x_grid.delta     = (x_grid.end-x_grid.start)/(double)(x_grid.num);
    N = M+3;
  }
  else {
    x_grid.delta     = (x_grid.end-x_grid.start)/(double)(x_grid.num-1);
    N = M+2;
  }

  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;

  spline->coefs = malloc (sizeof(double)*N);
  find_coefs_1d_d (spline->x_grid, xBC, data, 1, spline->coefs, 1);

  return spline;
}

void recompute_UBspline_1d_d (UBspline_1d_d* spline, double *data) {
  find_coefs_1d_d (spline->x_grid, spline->xBC, data, 1, spline->coefs, 1);
}

UBspline_2d_d* create_UBspline_2d_d (Ugrid x_grid, Ugrid y_grid,
                                     BCtype_d xBC, BCtype_d yBC, double *data) {
  /* Create new spline */
    int Mx, My;
    int Nx, Ny, ix, iy;

  UBspline_2d_d* spline = malloc (sizeof(UBspline_2d_d));
  spline->spcode = U2D;
  spline->tcode  = DOUBLE_REAL;
  spline->xBC = xBC; 
  spline->yBC = yBC; 

  /* Setup internal variables */
  Mx = x_grid.num;
  My = y_grid.num;

  if (xBC.lCode == PERIODIC)     Nx = Mx+3;
  else                           Nx = Mx+2;
  x_grid.delta = (x_grid.end - x_grid.start)/(double)(Nx-3);
  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC)     Ny = My+3;
  else                           Ny = My+2;
  y_grid.delta = (y_grid.end - y_grid.start)/(double)(Ny-3);
  y_grid.delta_inv = 1.0/y_grid.delta;
  spline->y_grid   = y_grid;
  spline->x_stride = Ny;

  spline->coefs = malloc (sizeof(double)*Nx*Ny);

  /* First, solve in the X-direction */
  for (iy=0; iy<My; iy++) {
    intptr_t doffset = iy;
    intptr_t coffset = iy;
    find_coefs_1d_d (spline->x_grid, xBC, data+doffset, My,
                     spline->coefs+coffset, Ny);
  }

  /* Now, solve in the Y-direction */
  for (ix=0; ix<Nx; ix++) {
    intptr_t doffset = ix*Ny;
    intptr_t coffset = ix*Ny;
    find_coefs_1d_d (spline->y_grid, yBC, spline->coefs+doffset, 1, 
                     spline->coefs+coffset, 1);
  }

  return spline;
}

void recompute_UBspline_2d_d (UBspline_2d_d* spline, double *data) {
  int Mx = spline->x_grid.num;
  int My = spline->y_grid.num;
  int Nx, Ny, ix, iy;

  if (spline->xBC.lCode == PERIODIC)     Nx = Mx+3;
  else                                   Nx = Mx+2;

  if (spline->yBC.lCode == PERIODIC)     Ny = My+3;
  else                                   Ny = My+2;

  /* First, solve in the X-direction */
  for (iy=0; iy<My; iy++) {
    intptr_t doffset = iy;
    intptr_t coffset = iy;
    find_coefs_1d_d (spline->x_grid, spline->xBC, data+doffset, My,
                     spline->coefs+coffset, Ny);
  }

  /* Now, solve in the Y-direction */
  for (ix=0; ix<Nx; ix++) {
    intptr_t doffset = ix*Ny;
    intptr_t coffset = ix*Ny;
    find_coefs_1d_d (spline->y_grid, spline->yBC, spline->coefs+doffset, 1, 
                     spline->coefs+coffset, 1);
  }
}

UBspline_3d_d* create_UBspline_3d_d (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
                                     BCtype_d xBC, BCtype_d yBC, BCtype_d zBC,
                                     double *data) {
  /* Create new spline */
    int Mx, My, Mz;
    int Nx, Ny, Nz, ix, iy, iz;

  UBspline_3d_d* spline = malloc (sizeof(UBspline_3d_d));
  spline->spcode = U3D;
  spline->tcode  = DOUBLE_REAL;
  spline->xBC = xBC; 
  spline->yBC = yBC; 
  spline->zBC = zBC; 

  /* Setup internal variables */
  Mx = x_grid.num;  My = y_grid.num; Mz = z_grid.num;

  if (xBC.lCode == PERIODIC)     Nx = Mx+3;
  else                           Nx = Mx+2;
  x_grid.delta = (x_grid.end - x_grid.start)/(double)(Nx-3);
  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC)     Ny = My+3;
  else                           Ny = My+2;
  y_grid.delta = (y_grid.end - y_grid.start)/(double)(Ny-3);
  y_grid.delta_inv = 1.0/y_grid.delta;
  spline->y_grid   = y_grid;

  if (zBC.lCode == PERIODIC)     Nz = Mz+3;
  else                           Nz = Mz+2;
  z_grid.delta = (z_grid.end - z_grid.start)/(double)(Nz-3);
  z_grid.delta_inv = 1.0/z_grid.delta;
  spline->z_grid   = z_grid;

  spline->x_stride = Ny*Nz;
  spline->y_stride = Nz;

  spline->coefs      = malloc (sizeof(double)*Nx*Ny*Nz);

  /* First, solve in the X-direction */
  for (iy=0; iy<My; iy++) 
    for (iz=0; iz<Mz; iz++) {
      intptr_t doffset = iy*Mz+iz;
      intptr_t coffset = iy*Nz+iz;
      find_coefs_1d_d (spline->x_grid, xBC, data+doffset, My*Mz,
                       spline->coefs+coffset, Ny*Nz);
    }

  /* Now, solve in the Y-direction */
  for (ix=0; ix<Nx; ix++) 
    for (iz=0; iz<Nz; iz++) {
      intptr_t doffset = ix*Ny*Nz + iz;
      intptr_t coffset = ix*Ny*Nz + iz;
      find_coefs_1d_d (spline->y_grid, yBC, spline->coefs+doffset, Nz, 
                       spline->coefs+coffset, Nz);
    }

  /* Now, solve in the Z-direction */
  for (ix=0; ix<Nx; ix++) 
    for (iy=0; iy<Ny; iy++) {
      intptr_t doffset = (ix*Ny+iy)*Nz;
      intptr_t coffset = (ix*Ny+iy)*Nz;
      find_coefs_1d_d (spline->z_grid, zBC, spline->coefs+doffset, 1, 
                       spline->coefs+coffset, 1);
    }

  return spline;
}

void recompute_UBspline_3d_d (UBspline_3d_d* spline, double *data) {
  int Mx = spline->x_grid.num;  
  int My = spline->y_grid.num; 
  int Mz = spline->z_grid.num;
  int Nx, Ny, Nz, ix, iy, iz;

  if (spline->xBC.lCode == PERIODIC)     Nx = Mx+3;
  else                                   Nx = Mx+2;
  if (spline->yBC.lCode == PERIODIC)     Ny = My+3;
  else                                   Ny = My+2;
  if (spline->zBC.lCode == PERIODIC)     Nz = Mz+3;
  else                                   Nz = Mz+2;

  /* First, solve in the X-direction */
  for (iy=0; iy<My; iy++) 
    for (iz=0; iz<Mz; iz++) {
      intptr_t doffset = iy*Mz+iz;
      intptr_t coffset = iy*Nz+iz;
      find_coefs_1d_d (spline->x_grid, spline->xBC, data+doffset, My*Mz,
                       spline->coefs+coffset, Ny*Nz);
    }

  /* Now, solve in the Y-direction */
  for (ix=0; ix<Nx; ix++) 
    for (iz=0; iz<Nz; iz++) {
      intptr_t doffset = ix*Ny*Nz + iz;
      intptr_t coffset = ix*Ny*Nz + iz;
      find_coefs_1d_d (spline->y_grid, spline->yBC, spline->coefs+doffset, Nz, 
                       spline->coefs+coffset, Nz);
    }

  /* Now, solve in the Z-direction */
  for (ix=0; ix<Nx; ix++) 
    for (iy=0; iy<Ny; iy++) {
      intptr_t doffset = (ix*Ny+iy)*Nz;
      intptr_t coffset = (ix*Ny+iy)*Nz;
      find_coefs_1d_d (spline->z_grid, spline->zBC, spline->coefs+doffset, 1, 
                       spline->coefs+coffset, 1);
    }
}

#ifndef NO_COMPLEX
/***********************************************************

        Double-Precision, Complex Creation Routines

 ***********************************************************

  On input, bands should be filled with:
  row 0   :  abcdInitial from boundary conditions
  rows 1:M:  basis functions in first 3 cols, data in last
  row M+1 :  abcdFinal   from boundary conditions
  cstride gives the stride between values in coefs.
  On exit, coefs with contain interpolating B-spline coefs */

UBspline_1d_z* create_UBspline_1d_z (Ugrid x_grid, BCtype_z xBC, sf_double_complex *data) {
    int M, N;

  /* Create new spline */
  UBspline_1d_z* spline = malloc (sizeof(UBspline_1d_z));
  spline->spcode = U1D;
  spline->tcode  = DOUBLE_COMPLEX;
  spline->xBC = xBC; 

  /* Setup internal variables */
  M = x_grid.num;

  if (xBC.lCode == PERIODIC) {
    x_grid.delta     = (x_grid.end-x_grid.start)/(double)(x_grid.num);
    N = M+3;
  }
  else {
    x_grid.delta     = (x_grid.end-x_grid.start)/(double)(x_grid.num-1);
    N = M+2;
  }

  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;
  spline->coefs = malloc (2*sizeof(double)*N);

  BCtype_d xBC_r, xBC_i;
  xBC_r.lCode = xBC.lCode;  xBC_r.rCode = xBC.rCode;
  xBC_r.lVal  = xBC.lVal_r; xBC_r.rVal  = xBC.rVal_r;
  xBC_i.lCode = xBC.lCode;  xBC_i.rCode = xBC.rCode;
  xBC_i.lVal  = xBC.lVal_i; xBC_i.rVal  = xBC.rVal_i;
  /* Real part */
  find_coefs_1d_d (spline->x_grid, xBC_r, (double*)data, 2, 
                   (double*)spline->coefs, 2);
  /* Imaginarty part */
  find_coefs_1d_d (spline->x_grid, xBC_i, ((double*)data)+1, 2, 
                   ((double*)spline->coefs)+1, 2);

  return spline;
}

void recompute_UBspline_1d_z (UBspline_1d_z* spline, sf_double_complex *data) {
  int M = spline->x_grid.num;
  int N;

  if (spline->xBC.lCode == PERIODIC)   N = M+3;
  else                                 N = M+2;

  BCtype_d xBC_r, xBC_i;
  xBC_r.lCode = spline->xBC.lCode;  xBC_r.rCode = spline->xBC.rCode;
  xBC_r.lVal  = spline->xBC.lVal_r; xBC_r.rVal  = spline->xBC.rVal_r;
  xBC_i.lCode = spline->xBC.lCode;  xBC_i.rCode = spline->xBC.rCode;
  xBC_i.lVal  = spline->xBC.lVal_i; xBC_i.rVal  = spline->xBC.rVal_i;
  /* Real part */
  find_coefs_1d_d (spline->x_grid, xBC_r, (double*)data, 2, 
                   (double*)spline->coefs, 2);
  /* Imaginarty part */
  find_coefs_1d_d (spline->x_grid, xBC_i, ((double*)data)+1, 2, 
                   ((double*)spline->coefs)+1, 2);
}

UBspline_2d_z* create_UBspline_2d_z (Ugrid x_grid, Ugrid y_grid,
                                     BCtype_z xBC, BCtype_z yBC, sf_double_complex *data) {
    int Mx, My;
    int Nx, Ny, ix, iy;

  /* Create new spline */
  UBspline_2d_z* spline = malloc (sizeof(UBspline_2d_z));
  spline->spcode = U2D;
  spline->tcode  = DOUBLE_COMPLEX;
  spline->xBC = xBC; 
  spline->yBC = yBC; 

  /* Setup internal variables */
  Mx = x_grid.num;
  My = y_grid.num;

  if (xBC.lCode == PERIODIC)     Nx = Mx+3;
  else                           Nx = Mx+2;
  x_grid.delta = (x_grid.end - x_grid.start)/(double)(Nx-3);
  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC)     Ny = My+3;
  else                           Ny = My+2;
  y_grid.delta = (y_grid.end - y_grid.start)/(double)(Ny-3);
  y_grid.delta_inv = 1.0/y_grid.delta;
  spline->y_grid   = y_grid;
  spline->x_stride = Ny;

  spline->coefs = malloc (2*sizeof(double)*Nx*Ny);

  BCtype_d xBC_r, xBC_i, yBC_r, yBC_i;
  xBC_r.lCode = xBC.lCode;  xBC_r.rCode = xBC.rCode;
  xBC_r.lVal  = xBC.lVal_r; xBC_r.rVal  = xBC.rVal_r;
  xBC_i.lCode = xBC.lCode;  xBC_i.rCode = xBC.rCode;
  xBC_i.lVal  = xBC.lVal_i; xBC_i.rVal  = xBC.rVal_i;
  yBC_r.lCode = yBC.lCode;  yBC_r.rCode = yBC.rCode;
  yBC_r.lVal  = yBC.lVal_r; yBC_r.rVal  = yBC.rVal_r;
  yBC_i.lCode = yBC.lCode;  yBC_i.rCode = yBC.rCode;
  yBC_i.lVal  = yBC.lVal_i; yBC_i.rVal  = yBC.rVal_i;
  /* First, solve in the X-direction */
  for (iy=0; iy<My; iy++) {
    intptr_t doffset = 2*iy;
    intptr_t coffset = 2*iy;
    /* Real part */
    find_coefs_1d_d (spline->x_grid, xBC_r, ((double*)data+doffset), 2*My,
                     (double*)spline->coefs+coffset, 2*Ny);
    /* Imag part */
    find_coefs_1d_d (spline->x_grid, xBC_i, ((double*)data)+doffset+1, 2*My,
                     ((double*)spline->coefs)+coffset+1, 2*Ny);
  }

  /* Now, solve in the Y-direction */
  for (ix=0; ix<Nx; ix++) {
    intptr_t doffset = 2*ix*Ny;
    intptr_t coffset = 2*ix*Ny;
    /* Real part */
    find_coefs_1d_d (spline->y_grid, yBC_r, ((double*)spline->coefs)+doffset, 2, 
                     (double*)spline->coefs+coffset, 2);
    /* Imag part */
    find_coefs_1d_d (spline->y_grid, yBC_i, (double*)spline->coefs+doffset+1, 2, 
                     ((double*)spline->coefs)+coffset+1, 2);
  }

  return spline;
}


void recompute_UBspline_2d_z (UBspline_2d_z* spline, sf_double_complex *data) {
  int Mx = spline->x_grid.num;
  int My = spline->y_grid.num;
  int Nx, Ny, ix, iy;

  if (spline->xBC.lCode == PERIODIC)     Nx = Mx+3;
  else                           Nx = Mx+2;
  if (spline->yBC.lCode == PERIODIC)     Ny = My+3;
  else                           Ny = My+2;

  BCtype_d xBC_r, xBC_i, yBC_r, yBC_i;
  xBC_r.lCode = spline->xBC.lCode;  xBC_r.rCode = spline->xBC.rCode;
  xBC_r.lVal  = spline->xBC.lVal_r; xBC_r.rVal  = spline->xBC.rVal_r;
  xBC_i.lCode = spline->xBC.lCode;  xBC_i.rCode = spline->xBC.rCode;
  xBC_i.lVal  = spline->xBC.lVal_i; xBC_i.rVal  = spline->xBC.rVal_i;
  yBC_r.lCode = spline->yBC.lCode;  yBC_r.rCode = spline->yBC.rCode;
  yBC_r.lVal  = spline->yBC.lVal_r; yBC_r.rVal  = spline->yBC.rVal_r;
  yBC_i.lCode = spline->yBC.lCode;  yBC_i.rCode = spline->yBC.rCode;
  yBC_i.lVal  = spline->yBC.lVal_i; yBC_i.rVal  = spline->yBC.rVal_i;

  /* First, solve in the X-direction */
  for (iy=0; iy<My; iy++) {
    intptr_t doffset = 2*iy;
    intptr_t coffset = 2*iy;
    /* Real part */
    find_coefs_1d_d (spline->x_grid, xBC_r, ((double*)data+doffset), 2*My,
                     (double*)spline->coefs+coffset, 2*Ny);
    /* Imag part */
    find_coefs_1d_d (spline->x_grid, xBC_i, ((double*)data)+doffset+1, 2*My,
                     ((double*)spline->coefs)+coffset+1, 2*Ny);
  }

  /* Now, solve in the Y-direction */
  for (ix=0; ix<Nx; ix++) {
    intptr_t doffset = 2*ix*Ny;
    intptr_t coffset = 2*ix*Ny;
    /* Real part */
    find_coefs_1d_d (spline->y_grid, yBC_r, ((double*)spline->coefs)+doffset, 2, 
                     (double*)spline->coefs+coffset, 2);
    /* Imag part */
    find_coefs_1d_d (spline->y_grid, yBC_i, (double*)spline->coefs+doffset+1, 2, 
                     ((double*)spline->coefs)+coffset+1, 2);
  }
}

UBspline_3d_z* create_UBspline_3d_z (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
                                     BCtype_z xBC, BCtype_z yBC, BCtype_z zBC,
                                     sf_double_complex *data) {
  /* Create new spline */
    int Mx, My, Mz;
    int Nx, Ny, Nz, ix, iy, iz;

  UBspline_3d_z* spline = malloc (sizeof(UBspline_3d_z));
  spline->spcode = U3D;
  spline->tcode  = DOUBLE_COMPLEX;
  spline->xBC = xBC; 
  spline->yBC = yBC; 
  spline->zBC = zBC;

  /* Setup internal variables */
  Mx = x_grid.num;  My = y_grid.num; Mz = z_grid.num;

  if (xBC.lCode == PERIODIC)     Nx = Mx+3;
  else                           Nx = Mx+2;
  x_grid.delta = (x_grid.end - x_grid.start)/(double)(Nx-3);
  x_grid.delta_inv = 1.0/x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC)     Ny = My+3;
  else                           Ny = My+2;
  y_grid.delta = (y_grid.end - y_grid.start)/(double)(Ny-3);
  y_grid.delta_inv = 1.0/y_grid.delta;
  spline->y_grid   = y_grid;

  if (zBC.lCode == PERIODIC)     Nz = Mz+3;
  else                           Nz = Mz+2;
  z_grid.delta = (z_grid.end - z_grid.start)/(double)(Nz-3);
  z_grid.delta_inv = 1.0/z_grid.delta;
  spline->z_grid   = z_grid;

  spline->x_stride = Ny*Nz;
  spline->y_stride = Nz;

  spline->coefs      = malloc (2*sizeof(double)*Nx*Ny*Nz);

  BCtype_d xBC_r, xBC_i, yBC_r, yBC_i, zBC_r, zBC_i;
  xBC_r.lCode = xBC.lCode;  xBC_r.rCode = xBC.rCode;
  xBC_r.lVal  = xBC.lVal_r; xBC_r.rVal  = xBC.rVal_r;
  xBC_i.lCode = xBC.lCode;  xBC_i.rCode = xBC.rCode;
  xBC_i.lVal  = xBC.lVal_i; xBC_i.rVal  = xBC.rVal_i;
  yBC_r.lCode = yBC.lCode;  yBC_r.rCode = yBC.rCode;
  yBC_r.lVal  = yBC.lVal_r; yBC_r.rVal  = yBC.rVal_r;
  yBC_i.lCode = yBC.lCode;  yBC_i.rCode = yBC.rCode;
  yBC_i.lVal  = yBC.lVal_i; yBC_i.rVal  = yBC.rVal_i;
  zBC_r.lCode = zBC.lCode;  zBC_r.rCode = zBC.rCode;
  zBC_r.lVal  = zBC.lVal_r; zBC_r.rVal  = zBC.rVal_r;
  zBC_i.lCode = zBC.lCode;  zBC_i.rCode = zBC.rCode;
  zBC_i.lVal  = zBC.lVal_i; zBC_i.rVal  = zBC.rVal_i;
  /* First, solve in the X-direction */
  for (iy=0; iy<My; iy++) 
    for (iz=0; iz<Mz; iz++) {
      intptr_t doffset = 2*(iy*Mz+iz);
      intptr_t coffset = 2*(iy*Nz+iz);
      /* Real part */
      find_coefs_1d_d (spline->x_grid, xBC_r, ((double*)data)+doffset, 2*My*Mz,
                       ((double*)spline->coefs)+coffset, 2*Ny*Nz);
      /* Imag part */
      find_coefs_1d_d (spline->x_grid, xBC_i, ((double*)data)+doffset+1, 2*My*Mz,
                       ((double*)spline->coefs)+coffset+1, 2*Ny*Nz);
    }

  /* Now, solve in the Y-direction */
  for (ix=0; ix<Nx; ix++) 
    for (iz=0; iz<Nz; iz++) {
      intptr_t doffset = 2*(ix*Ny*Nz + iz);
      intptr_t coffset = 2*(ix*Ny*Nz + iz);
      /* Real part */
      find_coefs_1d_d (spline->y_grid, yBC_r, ((double*)spline->coefs)+doffset, 2*Nz, 
                       ((double*)spline->coefs)+coffset, 2*Nz);
      /* Imag part */
      find_coefs_1d_d (spline->y_grid, yBC_i, ((double*)spline->coefs)+doffset+1, 2*Nz, 
                       ((double*)spline->coefs)+coffset+1, 2*Nz);
    }

  /* Now, solve in the Z-direction */
  for (ix=0; ix<Nx; ix++) 
    for (iy=0; iy<Ny; iy++) {
      intptr_t doffset = 2*((ix*Ny+iy)*Nz);
      intptr_t coffset = 2*((ix*Ny+iy)*Nz);
      /* Real part */
      find_coefs_1d_d (spline->z_grid, zBC_r, ((double*)spline->coefs)+doffset, 2, 
                       ((double*)spline->coefs)+coffset, 2);
      /* Imag part */
      find_coefs_1d_d (spline->z_grid, zBC_i, ((double*)spline->coefs)+doffset+1, 2, 
                       ((double*)spline->coefs)+coffset+1, 2);
    }
  return spline;
}

void recompute_UBspline_3d_z (UBspline_3d_z* spline, sf_double_complex *data) {
  /* Setup internal variables */
  int Mx = spline->x_grid.num;  
  int My = spline->y_grid.num; 
  int Mz = spline->z_grid.num;
  int Nx, Ny, Nz, ix, iy, iz;

  if (spline->xBC.lCode == PERIODIC)     Nx = Mx+3;
  else                                   Nx = Mx+2;
  if (spline->yBC.lCode == PERIODIC)     Ny = My+3;
  else                                   Ny = My+2;
  if (spline->zBC.lCode == PERIODIC)     Nz = Mz+3;
  else                                   Nz = Mz+2;

  BCtype_d xBC_r, xBC_i, yBC_r, yBC_i, zBC_r, zBC_i;
  xBC_r.lCode = spline->xBC.lCode;  xBC_r.rCode = spline->xBC.rCode;
  xBC_r.lVal  = spline->xBC.lVal_r; xBC_r.rVal  = spline->xBC.rVal_r;
  xBC_i.lCode = spline->xBC.lCode;  xBC_i.rCode = spline->xBC.rCode;
  xBC_i.lVal  = spline->xBC.lVal_i; xBC_i.rVal  = spline->xBC.rVal_i;
  yBC_r.lCode = spline->yBC.lCode;  yBC_r.rCode = spline->yBC.rCode;
  yBC_r.lVal  = spline->yBC.lVal_r; yBC_r.rVal  = spline->yBC.rVal_r;
  yBC_i.lCode = spline->yBC.lCode;  yBC_i.rCode = spline->yBC.rCode;
  yBC_i.lVal  = spline->yBC.lVal_i; yBC_i.rVal  = spline->yBC.rVal_i;
  zBC_r.lCode = spline->zBC.lCode;  zBC_r.rCode = spline->zBC.rCode;
  zBC_r.lVal  = spline->zBC.lVal_r; zBC_r.rVal  = spline->zBC.rVal_r;
  zBC_i.lCode = spline->zBC.lCode;  zBC_i.rCode = spline->zBC.rCode;
  zBC_i.lVal  = spline->zBC.lVal_i; zBC_i.rVal  = spline->zBC.rVal_i;
  /* First, solve in the X-direction */
  for (iy=0; iy<My; iy++) 
    for (iz=0; iz<Mz; iz++) {
      intptr_t doffset = 2*(iy*Mz+iz);
      intptr_t coffset = 2*(iy*Nz+iz);
      /* Real part */
      find_coefs_1d_d (spline->x_grid, xBC_r, ((double*)data)+doffset, 2*My*Mz,
                       ((double*)spline->coefs)+coffset, 2*Ny*Nz);
      /* Imag part */
      find_coefs_1d_d (spline->x_grid, xBC_i, ((double*)data)+doffset+1, 2*My*Mz,
                       ((double*)spline->coefs)+coffset+1, 2*Ny*Nz);
    }

  /* Now, solve in the Y-direction */
  for (ix=0; ix<Nx; ix++) 
    for (iz=0; iz<Nz; iz++) {
      intptr_t doffset = 2*(ix*Ny*Nz + iz);
      intptr_t coffset = 2*(ix*Ny*Nz + iz);
      /* Real part */
      find_coefs_1d_d (spline->y_grid, yBC_r, ((double*)spline->coefs)+doffset, 2*Nz, 
                       ((double*)spline->coefs)+coffset, 2*Nz);
      /* Imag part */
      find_coefs_1d_d (spline->y_grid, yBC_i, ((double*)spline->coefs)+doffset+1, 2*Nz, 
                       ((double*)spline->coefs)+coffset+1, 2*Nz);
    }

  /* Now, solve in the Z-direction */
  for (ix=0; ix<Nx; ix++) 
    for (iy=0; iy<Ny; iy++) {
      intptr_t doffset = 2*((ix*Ny+iy)*Nz);
      intptr_t coffset = 2*((ix*Ny+iy)*Nz);
      /* Real part */
      find_coefs_1d_d (spline->z_grid, zBC_r, ((double*)spline->coefs)+doffset, 2, 
                       ((double*)spline->coefs)+coffset, 2);
      /* Imag part */
      find_coefs_1d_d (spline->z_grid, zBC_i, ((double*)spline->coefs)+doffset+1, 2, 
                       ((double*)spline->coefs)+coffset+1, 2);
    }
}
#endif /* NO_COMPLEX */

void
destroy_UBspline (Bspline *spline) {
  free (spline->coefs);
  free (spline);
}

void destroy_Bspline (void *spline) {
  Bspline *sp = (Bspline *)spline;
  if (sp->sp_code <= U3D)
    destroy_UBspline (sp);
  else
    fprintf (stderr, "Error in destroy_Bspline:  invalide spline code %d.\n",
             sp->sp_code);
}

/*
 * bspline_eval_std_s|d|c|z.h
 */

/************************************************************/
/* 1D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_1d_s (UBspline_1d_s * spline, 
                         double x, float* val) {
  float u;
  float ipart, t;
  int i = (int) ipart;

  float tp[4];
  float* coefs = spline->coefs;

  x -= spline->x_grid.start;

  u = x*spline->x_grid.delta_inv;
  t = modff (u, &ipart);
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;

  *val = 
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
}

/* Value and first derivative */
void eval_UBspline_1d_s_vg (UBspline_1d_s * spline, double x, 
                            float* val, float* grad) {
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;

  float tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  float* coefs = spline->coefs;

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
void eval_UBspline_1d_s_vgl (UBspline_1d_s * spline, double x, 
                             float* val, float* grad,
                             float* lapl) {
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;

  float tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  float* coefs = spline->coefs;

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

void eval_UBspline_1d_s_vgh (UBspline_1d_s * spline, double x, 
                             float* val, float* grad,
                             float* hess) {
  eval_UBspline_1d_s_vgl (spline, x, val, grad, hess);
}

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_2d_s (UBspline_2d_s * spline, 
                         double x, double y, float* val) {
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
  float* coefs = spline->coefs;

  a[0] = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1] = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2] = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3] = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);

  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);

  int xs = spline->x_stride;
#define C(i,j) coefs[(ix+(i))*xs+iy+(j)]
  *val = (a[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
	  a[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
	  a[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
	  a[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));
#undef C

}

/* Value and gradient */
void eval_UBspline_2d_s_vg (UBspline_2d_s * spline, 
                            double x, double y, 
                            float* val, float* grad) {
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
  float* coefs = spline->coefs;

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

  int xs = spline->x_stride;
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
void eval_UBspline_2d_s_vgl (UBspline_2d_s * spline, 
                             double x, double y, float* val, 
                             float* grad, float* lapl) {
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
  float* coefs = spline->coefs;

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

  int xs = spline->x_stride;
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
void eval_UBspline_2d_s_vgh (UBspline_2d_s * spline, 
                             double x, double y, float* val, 
                             float* grad, float* hess) {
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
  float* coefs = spline->coefs;

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

  int xs = spline->x_stride;
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
void eval_UBspline_3d_s (UBspline_3d_s * spline, 
                         double x, double y, double z,
                         float* val) {
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
  float* coefs = spline->coefs;

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

  int xs = spline->x_stride;
  int ys = spline->y_stride;
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
void eval_UBspline_3d_s_vg (UBspline_3d_s * spline, 
                            double x, double y, double z,
                            float* val, float* grad) {
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
  float* coefs = spline->coefs;

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

  int xs = spline->x_stride;
  int ys = spline->y_stride;
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
void eval_UBspline_3d_s_vgl (UBspline_3d_s * spline, 
                             double x, double y, double z,
                             float* val, float* grad, float* lapl) {
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
  float* coefs = spline->coefs;

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

  int xs = spline->x_stride;
  int ys = spline->y_stride;
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
void eval_UBspline_3d_s_vgh (UBspline_3d_s * spline, 
                             double x, double y, double z,
                             float* val, float* grad, float* hess) {
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

  float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], da[4], db[4], dc[4], 
    d2a[4], d2b[4], d2c[4], cP[16], dcP[16], d2cP[16], bcP[4], dbcP[4],
    d2bcP[4], dbdcP[4], bd2cP[4], bdcP[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  float* coefs = spline->coefs;

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

  int xs = spline->x_stride;
  int ys = spline->y_stride;

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

/************************************************************/
/* 1D double-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_1d_d (UBspline_1d_d * spline, 
                         double x, double* val) {
  x -= spline->x_grid.start;
  double u = x*spline->x_grid.delta_inv;
  double ipart, t;
  t = modf (u, &ipart);
  int i = (int) ipart;

  double tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  double* coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Ad[ 0]*tp[0] + Ad[ 1]*tp[1] + Ad[ 2]*tp[2] + Ad[ 3]*tp[3])+
     coefs[i+1]*(Ad[ 4]*tp[0] + Ad[ 5]*tp[1] + Ad[ 6]*tp[2] + Ad[ 7]*tp[3])+
     coefs[i+2]*(Ad[ 8]*tp[0] + Ad[ 9]*tp[1] + Ad[10]*tp[2] + Ad[11]*tp[3])+
     coefs[i+3]*(Ad[12]*tp[0] + Ad[13]*tp[1] + Ad[14]*tp[2] + Ad[15]*tp[3]));
}

/* Value and first derivative */
void eval_UBspline_1d_d_vg (UBspline_1d_d * spline, double x, 
                            double* val, double* grad) {
  x -= spline->x_grid.start;
  double u = x*spline->x_grid.delta_inv;
  double ipart, t;
  t = modf (u, &ipart);
  int i = (int) ipart;

  double tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  double* coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Ad[ 0]*tp[0] + Ad[ 1]*tp[1] + Ad[ 2]*tp[2] + Ad[ 3]*tp[3])+
     coefs[i+1]*(Ad[ 4]*tp[0] + Ad[ 5]*tp[1] + Ad[ 6]*tp[2] + Ad[ 7]*tp[3])+
     coefs[i+2]*(Ad[ 8]*tp[0] + Ad[ 9]*tp[1] + Ad[10]*tp[2] + Ad[11]*tp[3])+
     coefs[i+3]*(Ad[12]*tp[0] + Ad[13]*tp[1] + Ad[14]*tp[2] + Ad[15]*tp[3]));
  *grad = spline->x_grid.delta_inv * 
    (coefs[i+0]*(dAd[ 1]*tp[1] + dAd[ 2]*tp[2] + dAd[ 3]*tp[3])+
     coefs[i+1]*(dAd[ 5]*tp[1] + dAd[ 6]*tp[2] + dAd[ 7]*tp[3])+
     coefs[i+2]*(dAd[ 9]*tp[1] + dAd[10]*tp[2] + dAd[11]*tp[3])+
     coefs[i+3]*(dAd[13]*tp[1] + dAd[14]*tp[2] + dAd[15]*tp[3]));
}

/* Value, first derivative, and second derivative */
void eval_UBspline_1d_d_vgl (UBspline_1d_d * spline, double x, 
                             double* val, double* grad,
                             double* lapl) {
  x -= spline->x_grid.start;
  double u = x*spline->x_grid.delta_inv;
  double ipart, t;
  t = modf (u, &ipart);
  int i = (int) ipart;

  double tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  double* coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Ad[ 0]*tp[0] + Ad[ 1]*tp[1] + Ad[ 2]*tp[2] + Ad[ 3]*tp[3])+
     coefs[i+1]*(Ad[ 4]*tp[0] + Ad[ 5]*tp[1] + Ad[ 6]*tp[2] + Ad[ 7]*tp[3])+
     coefs[i+2]*(Ad[ 8]*tp[0] + Ad[ 9]*tp[1] + Ad[10]*tp[2] + Ad[11]*tp[3])+
     coefs[i+3]*(Ad[12]*tp[0] + Ad[13]*tp[1] + Ad[14]*tp[2] + Ad[15]*tp[3]));
  *grad = spline->x_grid.delta_inv * 
    (coefs[i+0]*(dAd[ 1]*tp[1] + dAd[ 2]*tp[2] + dAd[ 3]*tp[3])+
     coefs[i+1]*(dAd[ 5]*tp[1] + dAd[ 6]*tp[2] + dAd[ 7]*tp[3])+
     coefs[i+2]*(dAd[ 9]*tp[1] + dAd[10]*tp[2] + dAd[11]*tp[3])+
     coefs[i+3]*(dAd[13]*tp[1] + dAd[14]*tp[2] + dAd[15]*tp[3]));
  *lapl = spline->x_grid.delta_inv * spline->x_grid.delta_inv * 
    (coefs[i+0]*(d2Ad[ 2]*tp[2] + d2Ad[ 3]*tp[3])+
     coefs[i+1]*(d2Ad[ 6]*tp[2] + d2Ad[ 7]*tp[3])+
     coefs[i+2]*(d2Ad[10]*tp[2] + d2Ad[11]*tp[3])+
     coefs[i+3]*(d2Ad[14]*tp[2] + d2Ad[15]*tp[3]));
}

void eval_UBspline_1d_d_vgh (UBspline_1d_d * spline, double x, 
                             double* val, double* grad,
                             double* hess) {
  eval_UBspline_1d_d_vgl (spline, x, val, grad, hess);
}

/************************************************************/
/* 2D double-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_2d_d (UBspline_2d_d * spline, 
                         double x, double y, double* val) {
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double ipartx, iparty, tx, ty;
  tx = modf (ux, &ipartx);
  ty = modf (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;

  double tpx[4], tpy[4], a[4], b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  double* coefs = spline->coefs;

  a[0] = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1] = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2] = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3] = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);

  b[0] = (Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1] = (Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2] = (Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3] = (Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);

  int xs = spline->x_stride;
#define C(i,j) coefs[(ix+(i))*xs+iy+(j)]
  *val = (a[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
	  a[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
	  a[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
	  a[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));
#undef C

}

/* Value and gradient */
void eval_UBspline_2d_d_vg (UBspline_2d_d * spline, 
                            double x, double y, 
                            double* val, double* grad) {
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double ipartx, iparty, tx, ty;
  tx = modf (ux, &ipartx);
  ty = modf (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;

  double tpx[4], tpy[4], a[4], b[4], da[4], db[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  double* coefs = spline->coefs;

  a[0]  = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]  = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]  = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]  = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  da[0] = (dAd[ 1]*tpx[1] + dAd[ 2]*tpx[2] + dAd[ 3]*tpx[3]);
  da[1] = (dAd[ 5]*tpx[1] + dAd[ 6]*tpx[2] + dAd[ 7]*tpx[3]);
  da[2] = (dAd[ 9]*tpx[1] + dAd[10]*tpx[2] + dAd[11]*tpx[3]);
  da[3] = (dAd[13]*tpx[1] + dAd[14]*tpx[2] + dAd[15]*tpx[3]);

  b[0]  = (Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1]  = (Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2]  = (Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3]  = (Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  db[0] = (dAd[ 1]*tpy[1] + dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
  db[1] = (dAd[ 5]*tpy[1] + dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
  db[2] = (dAd[ 9]*tpy[1] + dAd[10]*tpy[2] + dAd[11]*tpy[3]);
  db[3] = (dAd[13]*tpy[1] + dAd[14]*tpy[2] + dAd[15]*tpy[3]);

  int xs = spline->x_stride;
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
void eval_UBspline_2d_d_vgl (UBspline_2d_d * spline, 
                             double x, double y, double* val, 
                             double* grad, double* lapl) {
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double ipartx, iparty, tx, ty;
  tx = modf (ux, &ipartx);
  ty = modf (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;

  double tpx[4], tpy[4], a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  double* coefs = spline->coefs;

  a[0]   = (  Ad[ 0]*tpx[0] +   Ad[ 1]*tpx[1] +  Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]   = (  Ad[ 4]*tpx[0] +   Ad[ 5]*tpx[1] +  Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]   = (  Ad[ 8]*tpx[0] +   Ad[ 9]*tpx[1] +  Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]   = (  Ad[12]*tpx[0] +   Ad[13]*tpx[1] +  Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  da[0]  = ( dAd[ 1]*tpx[1] +  dAd[ 2]*tpx[2] + dAd[ 3]*tpx[3]);
  da[1]  = ( dAd[ 5]*tpx[1] +  dAd[ 6]*tpx[2] + dAd[ 7]*tpx[3]);
  da[2]  = ( dAd[ 9]*tpx[1] +  dAd[10]*tpx[2] + dAd[11]*tpx[3]);
  da[3]  = ( dAd[13]*tpx[1] +  dAd[14]*tpx[2] + dAd[15]*tpx[3]);
  d2a[0] = (d2Ad[ 2]*tpx[2] + d2Ad[ 3]*tpx[3]);
  d2a[1] = (d2Ad[ 6]*tpx[2] + d2Ad[ 7]*tpx[3]);
  d2a[2] = (d2Ad[10]*tpx[2] + d2Ad[11]*tpx[3]);
  d2a[3] = (d2Ad[14]*tpx[2] + d2Ad[15]*tpx[3]);

  b[0]  = ( Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1]  = ( Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2]  = ( Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3]  = ( Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  db[0] = (dAd[ 1]*tpy[1] + dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
  db[1] = (dAd[ 5]*tpy[1] + dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
  db[2] = (dAd[ 9]*tpy[1] + dAd[10]*tpy[2] + dAd[11]*tpy[3]);
  db[3] = (dAd[13]*tpy[1] + dAd[14]*tpy[2] + dAd[15]*tpy[3]);
  d2b[0] = (d2Ad[ 2]*tpy[2] + d2Ad[ 3]*tpy[3]);
  d2b[1] = (d2Ad[ 6]*tpy[2] + d2Ad[ 7]*tpy[3]);
  d2b[2] = (d2Ad[10]*tpy[2] + d2Ad[11]*tpy[3]);
  d2b[3] = (d2Ad[14]*tpy[2] + d2Ad[15]*tpy[3]);

  int xs = spline->x_stride;
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
void eval_UBspline_2d_d_vgh (UBspline_2d_d * spline, 
                             double x, double y, double* val, 
                             double* grad, double* hess) {
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double ipartx, iparty, tx, ty;
  tx = modf (ux, &ipartx);
  ty = modf (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;

  double tpx[4], tpy[4], a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  double* coefs = spline->coefs;

  a[0]   = (  Ad[ 0]*tpx[0] +   Ad[ 1]*tpx[1] +  Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]   = (  Ad[ 4]*tpx[0] +   Ad[ 5]*tpx[1] +  Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]   = (  Ad[ 8]*tpx[0] +   Ad[ 9]*tpx[1] +  Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]   = (  Ad[12]*tpx[0] +   Ad[13]*tpx[1] +  Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  da[0]  = ( dAd[ 1]*tpx[1] +  dAd[ 2]*tpx[2] + dAd[ 3]*tpx[3]);
  da[1]  = ( dAd[ 5]*tpx[1] +  dAd[ 6]*tpx[2] + dAd[ 7]*tpx[3]);
  da[2]  = ( dAd[ 9]*tpx[1] +  dAd[10]*tpx[2] + dAd[11]*tpx[3]);
  da[3]  = ( dAd[13]*tpx[1] +  dAd[14]*tpx[2] + dAd[15]*tpx[3]);
  d2a[0] = (d2Ad[ 2]*tpx[2] + d2Ad[ 3]*tpx[3]);
  d2a[1] = (d2Ad[ 6]*tpx[2] + d2Ad[ 7]*tpx[3]);
  d2a[2] = (d2Ad[10]*tpx[2] + d2Ad[11]*tpx[3]);
  d2a[3] = (d2Ad[14]*tpx[2] + d2Ad[15]*tpx[3]);

  b[0]   = (  Ad[ 0]*tpy[0] +   Ad[ 1]*tpy[1] +  Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1]   = (  Ad[ 4]*tpy[0] +   Ad[ 5]*tpy[1] +  Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2]   = (  Ad[ 8]*tpy[0] +   Ad[ 9]*tpy[1] +  Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3]   = (  Ad[12]*tpy[0] +   Ad[13]*tpy[1] +  Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  db[0]  = ( dAd[ 1]*tpy[1] +  dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
  db[1]  = ( dAd[ 5]*tpy[1] +  dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
  db[2]  = ( dAd[ 9]*tpy[1] +  dAd[10]*tpy[2] + dAd[11]*tpy[3]);
  db[3]  = ( dAd[13]*tpy[1] +  dAd[14]*tpy[2] + dAd[15]*tpy[3]);
  d2b[0] = (d2Ad[ 2]*tpy[2] + d2Ad[ 3]*tpy[3]);
  d2b[1] = (d2Ad[ 6]*tpy[2] + d2Ad[ 7]*tpy[3]);
  d2b[2] = (d2Ad[10]*tpy[2] + d2Ad[11]*tpy[3]);
  d2b[3] = (d2Ad[14]*tpy[2] + d2Ad[15]*tpy[3]);

  int xs = spline->x_stride;
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
/* 3D double-precision, real evaulation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_3d_d (UBspline_3d_d * spline, 
                         double x, double y, double z,
                         double* val) {
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  tz = modf (uz, &ipartz);  int iz = (int) ipartz;




  double tpx[4], tpy[4], tpz[4], a[4], b[4], c[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  double* coefs = spline->coefs;

  a[0] = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1] = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2] = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3] = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);

  b[0] = (Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1] = (Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2] = (Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3] = (Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);

  c[0] = (Ad[ 0]*tpz[0] + Ad[ 1]*tpz[1] + Ad[ 2]*tpz[2] + Ad[ 3]*tpz[3]);
  c[1] = (Ad[ 4]*tpz[0] + Ad[ 5]*tpz[1] + Ad[ 6]*tpz[2] + Ad[ 7]*tpz[3]);
  c[2] = (Ad[ 8]*tpz[0] + Ad[ 9]*tpz[1] + Ad[10]*tpz[2] + Ad[11]*tpz[3]);
  c[3] = (Ad[12]*tpz[0] + Ad[13]*tpz[1] + Ad[14]*tpz[2] + Ad[15]*tpz[3]);

  int xs = spline->x_stride;
  int ys = spline->y_stride;
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
void eval_UBspline_3d_d_vg (UBspline_3d_d * spline, 
                            double x, double y, double z,
                            double* val, double* grad) {
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;  
  ty = modf (uy, &iparty);  int iy = (int) iparty; 
  tz = modf (uz, &ipartz);  int iz = (int) ipartz; 

  double tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], da[4], db[4], dc[4], 
    cP[16], bcP[4], dbcP[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  double* coefs = spline->coefs;

  a[0]   = (  Ad[ 0]*tpx[0] +   Ad[ 1]*tpx[1] +  Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]   = (  Ad[ 4]*tpx[0] +   Ad[ 5]*tpx[1] +  Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]   = (  Ad[ 8]*tpx[0] +   Ad[ 9]*tpx[1] +  Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]   = (  Ad[12]*tpx[0] +   Ad[13]*tpx[1] +  Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  da[0]  = ( dAd[ 1]*tpx[1] +  dAd[ 2]*tpx[2] + dAd[ 3]*tpx[3]);
  da[1]  = ( dAd[ 5]*tpx[1] +  dAd[ 6]*tpx[2] + dAd[ 7]*tpx[3]);
  da[2]  = ( dAd[ 9]*tpx[1] +  dAd[10]*tpx[2] + dAd[11]*tpx[3]);
  da[3]  = ( dAd[13]*tpx[1] +  dAd[14]*tpx[2] + dAd[15]*tpx[3]);

  b[0]  = ( Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1]  = ( Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2]  = ( Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3]  = ( Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  db[0] = (dAd[ 1]*tpy[1] + dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
  db[1] = (dAd[ 5]*tpy[1] + dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
  db[2] = (dAd[ 9]*tpy[1] + dAd[10]*tpy[2] + dAd[11]*tpy[3]);
  db[3] = (dAd[13]*tpy[1] + dAd[14]*tpy[2] + dAd[15]*tpy[3]);

  c[0]  = ( Ad[ 0]*tpz[0] + Ad[ 1]*tpz[1] + Ad[ 2]*tpz[2] + Ad[ 3]*tpz[3]);
  c[1]  = ( Ad[ 4]*tpz[0] + Ad[ 5]*tpz[1] + Ad[ 6]*tpz[2] + Ad[ 7]*tpz[3]);
  c[2]  = ( Ad[ 8]*tpz[0] + Ad[ 9]*tpz[1] + Ad[10]*tpz[2] + Ad[11]*tpz[3]);
  c[3]  = ( Ad[12]*tpz[0] + Ad[13]*tpz[1] + Ad[14]*tpz[2] + Ad[15]*tpz[3]);
  dc[0] = (dAd[ 1]*tpz[1] + dAd[ 2]*tpz[2] + dAd[ 3]*tpz[3]);
  dc[1] = (dAd[ 5]*tpz[1] + dAd[ 6]*tpz[2] + dAd[ 7]*tpz[3]);
  dc[2] = (dAd[ 9]*tpz[1] + dAd[10]*tpz[2] + dAd[11]*tpz[3]);
  dc[3] = (dAd[13]*tpz[1] + dAd[14]*tpz[2] + dAd[15]*tpz[3]);

  int xs = spline->x_stride;
  int ys = spline->y_stride;
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
void eval_UBspline_3d_d_vgl (UBspline_3d_d * spline, 
                             double x, double y, double z,
                             double* val, double* grad, double* lapl) {
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;  
  ty = modf (uy, &iparty);  int iy = (int) iparty; 
  tz = modf (uz, &ipartz);  int iz = (int) ipartz; 

  double tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], da[4], db[4], dc[4], 
    d2a[4], d2b[4], d2c[4], cP[16], dcP[16], bcP[4], dbcP[4], d2bcP[4], bdcP[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  double* coefs = spline->coefs;

  a[0]   = (  Ad[ 0]*tpx[0] +   Ad[ 1]*tpx[1] +  Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]   = (  Ad[ 4]*tpx[0] +   Ad[ 5]*tpx[1] +  Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]   = (  Ad[ 8]*tpx[0] +   Ad[ 9]*tpx[1] +  Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]   = (  Ad[12]*tpx[0] +   Ad[13]*tpx[1] +  Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  da[0]  = ( dAd[ 1]*tpx[1] +  dAd[ 2]*tpx[2] + dAd[ 3]*tpx[3]);
  da[1]  = ( dAd[ 5]*tpx[1] +  dAd[ 6]*tpx[2] + dAd[ 7]*tpx[3]);
  da[2]  = ( dAd[ 9]*tpx[1] +  dAd[10]*tpx[2] + dAd[11]*tpx[3]);
  da[3]  = ( dAd[13]*tpx[1] +  dAd[14]*tpx[2] + dAd[15]*tpx[3]);
  d2a[0] = (d2Ad[ 2]*tpx[2] + d2Ad[ 3]*tpx[3]);
  d2a[1] = (d2Ad[ 6]*tpx[2] + d2Ad[ 7]*tpx[3]);
  d2a[2] = (d2Ad[10]*tpx[2] + d2Ad[11]*tpx[3]);
  d2a[3] = (d2Ad[14]*tpx[2] + d2Ad[15]*tpx[3]);

  b[0]  = ( Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1]  = ( Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2]  = ( Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3]  = ( Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  db[0] = (dAd[ 1]*tpy[1] + dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
  db[1] = (dAd[ 5]*tpy[1] + dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
  db[2] = (dAd[ 9]*tpy[1] + dAd[10]*tpy[2] + dAd[11]*tpy[3]);
  db[3] = (dAd[13]*tpy[1] + dAd[14]*tpy[2] + dAd[15]*tpy[3]);
  d2b[0] = (d2Ad[ 2]*tpy[2] + d2Ad[ 3]*tpy[3]);
  d2b[1] = (d2Ad[ 6]*tpy[2] + d2Ad[ 7]*tpy[3]);
  d2b[2] = (d2Ad[10]*tpy[2] + d2Ad[11]*tpy[3]);
  d2b[3] = (d2Ad[14]*tpy[2] + d2Ad[15]*tpy[3]);

  c[0]  = ( Ad[ 0]*tpz[0] + Ad[ 1]*tpz[1] + Ad[ 2]*tpz[2] + Ad[ 3]*tpz[3]);
  c[1]  = ( Ad[ 4]*tpz[0] + Ad[ 5]*tpz[1] + Ad[ 6]*tpz[2] + Ad[ 7]*tpz[3]);
  c[2]  = ( Ad[ 8]*tpz[0] + Ad[ 9]*tpz[1] + Ad[10]*tpz[2] + Ad[11]*tpz[3]);
  c[3]  = ( Ad[12]*tpz[0] + Ad[13]*tpz[1] + Ad[14]*tpz[2] + Ad[15]*tpz[3]);
  dc[0] = (dAd[ 1]*tpz[1] + dAd[ 2]*tpz[2] + dAd[ 3]*tpz[3]);
  dc[1] = (dAd[ 5]*tpz[1] + dAd[ 6]*tpz[2] + dAd[ 7]*tpz[3]);
  dc[2] = (dAd[ 9]*tpz[1] + dAd[10]*tpz[2] + dAd[11]*tpz[3]);
  dc[3] = (dAd[13]*tpz[1] + dAd[14]*tpz[2] + dAd[15]*tpz[3]);
  d2c[0] = (d2Ad[ 2]*tpz[2] + d2Ad[ 3]*tpz[3]);
  d2c[1] = (d2Ad[ 6]*tpz[2] + d2Ad[ 7]*tpz[3]);
  d2c[2] = (d2Ad[10]*tpz[2] + d2Ad[11]*tpz[3]);
  d2c[3] = (d2Ad[14]*tpz[2] + d2Ad[15]*tpz[3]);

  int xs = spline->x_stride;
  int ys = spline->y_stride;
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
void eval_UBspline_3d_d_vgh (UBspline_3d_d * spline, 
                             double x, double y, double z,
                             double* val, double* grad, double* hess) {
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;  
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
  uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
  uz = fmin (uz, (double)(spline->z_grid.num)-1.0e-5);
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  tz = modf (uz, &ipartz);  int iz = (int) ipartz;

  double tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], da[4], db[4], dc[4], 
    d2a[4], d2b[4], d2c[4], cP[16], dcP[16], d2cP[16], bcP[4], dbcP[4],
    d2bcP[4], dbdcP[4], bd2cP[4], bdcP[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  double* coefs = spline->coefs;

  a[0]   = (  Ad[ 0]*tpx[0] +   Ad[ 1]*tpx[1] +  Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]   = (  Ad[ 4]*tpx[0] +   Ad[ 5]*tpx[1] +  Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]   = (  Ad[ 8]*tpx[0] +   Ad[ 9]*tpx[1] +  Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]   = (  Ad[12]*tpx[0] +   Ad[13]*tpx[1] +  Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  da[0]  = ( dAd[ 1]*tpx[1] +  dAd[ 2]*tpx[2] + dAd[ 3]*tpx[3]);
  da[1]  = ( dAd[ 5]*tpx[1] +  dAd[ 6]*tpx[2] + dAd[ 7]*tpx[3]);
  da[2]  = ( dAd[ 9]*tpx[1] +  dAd[10]*tpx[2] + dAd[11]*tpx[3]);
  da[3]  = ( dAd[13]*tpx[1] +  dAd[14]*tpx[2] + dAd[15]*tpx[3]);
  d2a[0] = (d2Ad[ 2]*tpx[2] + d2Ad[ 3]*tpx[3]);
  d2a[1] = (d2Ad[ 6]*tpx[2] + d2Ad[ 7]*tpx[3]);
  d2a[2] = (d2Ad[10]*tpx[2] + d2Ad[11]*tpx[3]);
  d2a[3] = (d2Ad[14]*tpx[2] + d2Ad[15]*tpx[3]);

  b[0]  = ( Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1]  = ( Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2]  = ( Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3]  = ( Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  db[0] = (dAd[ 1]*tpy[1] + dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
  db[1] = (dAd[ 5]*tpy[1] + dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
  db[2] = (dAd[ 9]*tpy[1] + dAd[10]*tpy[2] + dAd[11]*tpy[3]);
  db[3] = (dAd[13]*tpy[1] + dAd[14]*tpy[2] + dAd[15]*tpy[3]);
  d2b[0] = (d2Ad[ 2]*tpy[2] + d2Ad[ 3]*tpy[3]);
  d2b[1] = (d2Ad[ 6]*tpy[2] + d2Ad[ 7]*tpy[3]);
  d2b[2] = (d2Ad[10]*tpy[2] + d2Ad[11]*tpy[3]);
  d2b[3] = (d2Ad[14]*tpy[2] + d2Ad[15]*tpy[3]);

  c[0]  = ( Ad[ 0]*tpz[0] + Ad[ 1]*tpz[1] + Ad[ 2]*tpz[2] + Ad[ 3]*tpz[3]);
  c[1]  = ( Ad[ 4]*tpz[0] + Ad[ 5]*tpz[1] + Ad[ 6]*tpz[2] + Ad[ 7]*tpz[3]);
  c[2]  = ( Ad[ 8]*tpz[0] + Ad[ 9]*tpz[1] + Ad[10]*tpz[2] + Ad[11]*tpz[3]);
  c[3]  = ( Ad[12]*tpz[0] + Ad[13]*tpz[1] + Ad[14]*tpz[2] + Ad[15]*tpz[3]);
  dc[0] = (dAd[ 1]*tpz[1] + dAd[ 2]*tpz[2] + dAd[ 3]*tpz[3]);
  dc[1] = (dAd[ 5]*tpz[1] + dAd[ 6]*tpz[2] + dAd[ 7]*tpz[3]);
  dc[2] = (dAd[ 9]*tpz[1] + dAd[10]*tpz[2] + dAd[11]*tpz[3]);
  dc[3] = (dAd[13]*tpz[1] + dAd[14]*tpz[2] + dAd[15]*tpz[3]);
  d2c[0] = (d2Ad[ 2]*tpz[2] + d2Ad[ 3]*tpz[3]);
  d2c[1] = (d2Ad[ 6]*tpz[2] + d2Ad[ 7]*tpz[3]);
  d2c[2] = (d2Ad[10]*tpz[2] + d2Ad[11]*tpz[3]);
  d2c[3] = (d2Ad[14]*tpz[2] + d2Ad[15]*tpz[3]);

  int xs = spline->x_stride;
  int ys = spline->y_stride;

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
  /* dx dz; */
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

#ifndef NO_COMPLEX
/************************************************************/
/* 1D single-precision, complex evaulation functions        */
/************************************************************/

/* Value only */
void eval_UBspline_1d_c (UBspline_1d_c * spline, 
                         double x, sf_complex* val) {
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;

  float tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  sf_complex* coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
}

/* Value and gradient */
void eval_UBspline_1d_c_vg (UBspline_1d_c * spline, double x, 
                            sf_complex* val, sf_complex* grad) {
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;

  float tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  sf_complex* coefs = spline->coefs;

  float dxInv = spline->x_grid.delta_inv;

  *val = 
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
  *grad = dxInv * 
    (coefs[i+0]*(dAf[ 1]*tp[1] + dAf[ 2]*tp[2] + dAf[ 3]*tp[3])+
     coefs[i+1]*(dAf[ 5]*tp[1] + dAf[ 6]*tp[2] + dAf[ 7]*tp[3])+
     coefs[i+2]*(dAf[ 9]*tp[1] + dAf[10]*tp[2] + dAf[11]*tp[3])+
     coefs[i+3]*(dAf[13]*tp[1] + dAf[14]*tp[2] + dAf[15]*tp[3]));
}

/* Value, gradient, and laplacian */
void eval_UBspline_1d_c_vgl (UBspline_1d_c * spline, double x, 
                             sf_complex* val, sf_complex* grad,
                             sf_complex* lapl) {
  x -= spline->x_grid.start;
  float u = x*spline->x_grid.delta_inv;
  float ipart, t;
  t = modff (u, &ipart);
  int i = (int) ipart;

  float tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  sf_complex* coefs = spline->coefs;

  float dxInv = spline->x_grid.delta_inv;

  *val = 
    (coefs[i+0]*(Af[ 0]*tp[0] + Af[ 1]*tp[1] + Af[ 2]*tp[2] + Af[ 3]*tp[3])+
     coefs[i+1]*(Af[ 4]*tp[0] + Af[ 5]*tp[1] + Af[ 6]*tp[2] + Af[ 7]*tp[3])+
     coefs[i+2]*(Af[ 8]*tp[0] + Af[ 9]*tp[1] + Af[10]*tp[2] + Af[11]*tp[3])+
     coefs[i+3]*(Af[12]*tp[0] + Af[13]*tp[1] + Af[14]*tp[2] + Af[15]*tp[3]));
  *grad = dxInv * 
    (coefs[i+0]*(dAf[ 1]*tp[1] + dAf[ 2]*tp[2] + dAf[ 3]*tp[3])+
     coefs[i+1]*(dAf[ 5]*tp[1] + dAf[ 6]*tp[2] + dAf[ 7]*tp[3])+
     coefs[i+2]*(dAf[ 9]*tp[1] + dAf[10]*tp[2] + dAf[11]*tp[3])+
     coefs[i+3]*(dAf[13]*tp[1] + dAf[14]*tp[2] + dAf[15]*tp[3]));
  *lapl = dxInv * dxInv * 
    (coefs[i+0]*(d2Af[ 2]*tp[2] + d2Af[ 3]*tp[3])+
     coefs[i+1]*(d2Af[ 6]*tp[2] + d2Af[ 7]*tp[3])+
     coefs[i+2]*(d2Af[10]*tp[2] + d2Af[11]*tp[3])+
     coefs[i+3]*(d2Af[14]*tp[2] + d2Af[15]*tp[3]));
}

/* Value, gradient, and Hessian */
void eval_UBspline_1d_c_vgh (UBspline_1d_c * spline, double x, 
                             sf_complex* val, 
                             sf_complex* grad,
                             sf_complex* hess) {
  eval_UBspline_1d_c_vgl (spline, x, val, grad, hess);
}

/************************************************************/
/* 2D single-precision, complex evaulation functions        */
/************************************************************/

/* Value only */
void eval_UBspline_2d_c (UBspline_2d_c * spline, 
                         double x, double y, sf_complex* val) {
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
  sf_complex* coefs = spline->coefs;

  a[0] = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
  a[1] = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
  a[2] = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
  a[3] = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);

  b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
  b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
  b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
  b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);

  int xs = spline->x_stride;
#define C(i,j) coefs[(ix+(i))*xs+iy+(j)]
  *val = (a[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
	  a[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
	  a[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
	  a[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));
#undef C

}

/* Value and gradient */
void eval_UBspline_2d_c_vg (UBspline_2d_c * spline, 
                            double x, double y, 
                            sf_complex* val, sf_complex* grad) {
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
  sf_complex* coefs = spline->coefs;

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

  int xs = spline->x_stride;
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
#define C(i,j) coefs[(ix+(i))*xs+iy+(j)]
  *val =    
    (a[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
     a[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
     a[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
     a[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));
  grad[0] = dxInv *
    (da[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
     da[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
     da[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
     da[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));
  grad[1] = dyInv * 
    (a[0]*(C(0,0)*db[0]+C(0,1)*db[1]+C(0,2)*db[2]+C(0,3)*db[3])+
     a[1]*(C(1,0)*db[0]+C(1,1)*db[1]+C(1,2)*db[2]+C(1,3)*db[3])+
     a[2]*(C(2,0)*db[0]+C(2,1)*db[1]+C(2,2)*db[2]+C(2,3)*db[3])+
     a[3]*(C(3,0)*db[0]+C(3,1)*db[1]+C(3,2)*db[2]+C(3,3)*db[3]));
#undef C

}

/* Value, gradient, and laplacian */
void eval_UBspline_2d_c_vgl (UBspline_2d_c * spline, 
                             double x, double y, sf_complex* val, 
                             sf_complex* grad, sf_complex* lapl) {
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
  sf_complex* coefs = spline->coefs;

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

  int xs = spline->x_stride;

  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;

#define C(i,j) coefs[(ix+(i))*xs+iy+(j)]
  *val =    
    (a[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
     a[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
     a[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
     a[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));
  grad[0] = dxInv *
    (da[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
     da[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
     da[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
     da[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));
  grad[1] = dyInv*
    (a[0]*(C(0,0)*db[0]+C(0,1)*db[1]+C(0,2)*db[2]+C(0,3)*db[3])+
     a[1]*(C(1,0)*db[0]+C(1,1)*db[1]+C(1,2)*db[2]+C(1,3)*db[3])+
     a[2]*(C(2,0)*db[0]+C(2,1)*db[1]+C(2,2)*db[2]+C(2,3)*db[3])+
     a[3]*(C(3,0)*db[0]+C(3,1)*db[1]+C(3,2)*db[2]+C(3,3)*db[3]));
  *lapl   = 
    dyInv * dyInv *
    (a[0]*(C(0,0)*d2b[0]+C(0,1)*d2b[1]+C(0,2)*d2b[2]+C(0,3)*d2b[3])+
      a[1]*(C(1,0)*d2b[0]+C(1,1)*d2b[1]+C(1,2)*d2b[2]+C(1,3)*d2b[3])+
      a[2]*(C(2,0)*d2b[0]+C(2,1)*d2b[1]+C(2,2)*d2b[2]+C(2,3)*d2b[3])+
     a[3]*(C(3,0)*d2b[0]+C(3,1)*d2b[1]+C(3,2)*d2b[2]+C(3,3)*d2b[3])) + 
    dxInv * dxInv *
     (d2a[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
      d2a[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
      d2a[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
      d2a[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));

#undef C

}

/* Value, gradient, and Hessian */
void eval_UBspline_2d_c_vgh (UBspline_2d_c * spline, 
                             double x, double y, sf_complex* val, 
                             sf_complex* grad, sf_complex* hess) {
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
  sf_complex* coefs = spline->coefs;

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

  int xs = spline->x_stride;
  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
#define C(i,j) coefs[(ix+(i))*xs+iy+(j)]
  *val =    
    (  a[0]*(C(0,0)*  b[0]+C(0,1)*  b[1]+C(0,2)*  b[2]+C(0,3)*  b[3])+
       a[1]*(C(1,0)*  b[0]+C(1,1)*  b[1]+C(1,2)*  b[2]+C(1,3)*  b[3])+
       a[2]*(C(2,0)*  b[0]+C(2,1)*  b[1]+C(2,2)*  b[2]+C(2,3)*  b[3])+
       a[3]*(C(3,0)*  b[0]+C(3,1)*  b[1]+C(3,2)*  b[2]+C(3,3)*  b[3]));
  grad[0] = dxInv * 
    ( da[0]*(C(0,0)*  b[0]+C(0,1)*  b[1]+C(0,2)*  b[2]+C(0,3)*  b[3])+
      da[1]*(C(1,0)*  b[0]+C(1,1)*  b[1]+C(1,2)*  b[2]+C(1,3)*  b[3])+
      da[2]*(C(2,0)*  b[0]+C(2,1)*  b[1]+C(2,2)*  b[2]+C(2,3)*  b[3])+
      da[3]*(C(3,0)*  b[0]+C(3,1)*  b[1]+C(3,2)*  b[2]+C(3,3)*  b[3]));
  grad[1] = dyInv *
    (  a[0]*(C(0,0)* db[0]+C(0,1)* db[1]+C(0,2)* db[2]+C(0,3)* db[3])+
       a[1]*(C(1,0)* db[0]+C(1,1)* db[1]+C(1,2)* db[2]+C(1,3)* db[3])+
       a[2]*(C(2,0)* db[0]+C(2,1)* db[1]+C(2,2)* db[2]+C(2,3)* db[3])+
       a[3]*(C(3,0)* db[0]+C(3,1)* db[1]+C(3,2)* db[2]+C(3,3)* db[3]));
  hess[0] = dxInv * dxInv *
    (d2a[0]*(C(0,0)*  b[0]+C(0,1)*  b[1]+C(0,2)*  b[2]+C(0,3)*  b[3])+
     d2a[1]*(C(1,0)*  b[0]+C(1,1)*  b[1]+C(1,2)*  b[2]+C(1,3)*  b[3])+
     d2a[2]*(C(2,0)*  b[0]+C(2,1)*  b[1]+C(2,2)*  b[2]+C(2,3)*  b[3])+
     d2a[3]*(C(3,0)*  b[0]+C(3,1)*  b[1]+C(3,2)*  b[2]+C(3,3)*  b[3]));
  hess[1] = dxInv * dyInv * 
    ( da[0]*(C(0,0)* db[0]+C(0,1)* db[1]+C(0,2)* db[2]+C(0,3)* db[3])+
      da[1]*(C(1,0)* db[0]+C(1,1)* db[1]+C(1,2)* db[2]+C(1,3)* db[3])+
      da[2]*(C(2,0)* db[0]+C(2,1)* db[1]+C(2,2)* db[2]+C(2,3)* db[3])+
      da[3]*(C(3,0)* db[0]+C(3,1)* db[1]+C(3,2)* db[2]+C(3,3)* db[3]));
  hess[3] = dyInv * dyInv * 
    (  a[0]*(C(0,0)*d2b[0]+C(0,1)*d2b[1]+C(0,2)*d2b[2]+C(0,3)*d2b[3])+
       a[1]*(C(1,0)*d2b[0]+C(1,1)*d2b[1]+C(1,2)*d2b[2]+C(1,3)*d2b[3])+
       a[2]*(C(2,0)*d2b[0]+C(2,1)*d2b[1]+C(2,2)*d2b[2]+C(2,3)*d2b[3])+
       a[3]*(C(3,0)*d2b[0]+C(3,1)*d2b[1]+C(3,2)*d2b[2]+C(3,3)*d2b[3]));
  hess[2] = hess[1];

#undef C

}

/************************************************************/
/* 3D single-precision, complex evaulation functions        */
/************************************************************/

/* Value only */
void eval_UBspline_3d_c (UBspline_3d_c * spline, 
                         double x, double y, double z,
                         sf_complex* val) {
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
  sf_complex* coefs = spline->coefs;

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

  int xs = spline->x_stride;
  int ys = spline->y_stride;
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
void eval_UBspline_3d_c_vg (UBspline_3d_c * spline, 
                            double x, double y, double z,
                            sf_complex* val, sf_complex* grad) {
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

  float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], da[4], db[4], dc[4];
  sf_complex cP[16], bcP[4], dbcP[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  sf_complex* coefs = spline->coefs;

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

  int xs = spline->x_stride;
  int ys = spline->y_stride;

  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv;

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
  grad[0] = dxInv * 
    (da[0]*bcP[0] + da[1]*bcP[1] + da[2]*bcP[2] + da[3]*bcP[3]);
  grad[1] = dyInv *
    (a[0]*dbcP[0] + a[1]*dbcP[1] + a[2]*dbcP[2] + a[3]*dbcP[3]);
  grad[2] = dzInv * 
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
void eval_UBspline_3d_c_vgl (UBspline_3d_c * spline, 
                             double x, double y, double z,
                             sf_complex* val, sf_complex* grad, 
                             sf_complex* lapl) {
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
    d2a[4], d2b[4], d2c[4];
  sf_complex cP[16], dcP[16], bcP[4], dbcP[4], d2bcP[4], bdcP[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  sf_complex* coefs = spline->coefs;

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

  int xs = spline->x_stride;
  int ys = spline->y_stride;

  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv;

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

  grad[0] = dxInv *
    (da[0]*bcP[0] + da[1]*bcP[1] + da[2]*bcP[2] + da[3]*bcP[3]);
  grad[1] = dyInv * 
    (a[0]*dbcP[0] + a[1]*dbcP[1] + a[2]*dbcP[2] + a[3]*dbcP[3]);
  grad[2] = dzInv * 
    (a[0]*bdcP[0] + a[1]*bdcP[1] + a[2]*bdcP[2] + a[3]*bdcP[3]);

  *lapl = 
    dxInv * dxInv * 
    (d2a[0]*bcP[0] + d2a[1]*bcP[1] + d2a[2]*bcP[2] + d2a[3]*bcP[3])

    + dyInv * dyInv * 
    (a[0]*d2bcP[0] + a[1]*d2bcP[1] + a[2]*d2bcP[2] + a[3]*d2bcP[3]) +

    + dzInv * dzInv * 
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
void eval_UBspline_3d_c_vgh (UBspline_3d_c * spline, 
                             double x, double y, double z,
                             sf_complex* val, sf_complex* grad, 
                             sf_complex* hess) {
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

  float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], da[4], db[4], dc[4], 
    d2a[4], d2b[4], d2c[4];
  sf_complex cP[16], dcP[16], d2cP[16], bcP[4], dbcP[4],
    d2bcP[4], dbdcP[4], bd2cP[4], bdcP[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  sf_complex* coefs = spline->coefs;

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

  int xs = spline->x_stride;
  int ys = spline->y_stride;

  float dxInv = spline->x_grid.delta_inv;
  float dyInv = spline->y_grid.delta_inv;
  float dzInv = spline->z_grid.delta_inv;

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
  grad[0] = dxInv *
    (da[0] *bcP[0] + da[1]*bcP[1] + da[2]*bcP[2] + da[3]*bcP[3]);
  grad[1] = dyInv *
    (a[0]*dbcP[0] + a[1]*dbcP[1] + a[2]*dbcP[2] + a[3]*dbcP[3]);
  grad[2] = dzInv *
    (a[0]*bdcP[0] + a[1]*bdcP[1] + a[2]*bdcP[2] + a[3]*bdcP[3]);
  /* d2x */
  hess[0] = dxInv * dxInv *
    (d2a[0]*bcP[0] + d2a[1]*bcP[1] + d2a[2]*bcP[2] + d2a[3]*bcP[3]);
  /* dx dy */
  hess[1] = dxInv * dyInv *
    (da[0]*dbcP[0] + da[1]*dbcP[1] + da[2]*dbcP[2] + da[3]*dbcP[3]);
  hess[3] = hess[1];
  /* dx dz */
  hess[2] = dxInv * dzInv *
    (da[0]*bdcP[0] + da[1]*bdcP[1] + da[2]*bdcP[2] + da[3]*bdcP[3]);
  hess[6] = hess[2];
  /* d2y */
  hess[4] = dyInv * dyInv *
    (a[0]*d2bcP[0] + a[1]*d2bcP[1] + a[2]*d2bcP[2] + a[3]*d2bcP[3]);
  /* dy dz */
  hess[5] = dyInv * dzInv *
    (a[0]*dbdcP[0] + a[1]*dbdcP[1] + a[2]*dbdcP[2] + a[3]*dbdcP[3]);
  hess[7] = hess[5];
  /* d2z */
  hess[8] = dzInv * dzInv *
    (a[0]*bd2cP[0] + a[1]*bd2cP[1] + a[2]*bd2cP[2] + a[3]*bd2cP[3]);
#undef P

}

/************************************************************/
/* 1D double-precision, complex evaulation functions        */
/************************************************************/

/* Value only */
void eval_UBspline_1d_z (UBspline_1d_z * spline, 
                         double x, sf_double_complex* val) {
  x -= spline->x_grid.start;
  double u = x*spline->x_grid.delta_inv;
  double ipart, t;
  t = modf (u, &ipart);
  int i = (int) ipart;

  double tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  sf_double_complex* coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Ad[ 0]*tp[0] + Ad[ 1]*tp[1] + Ad[ 2]*tp[2] + Ad[ 3]*tp[3])+
     coefs[i+1]*(Ad[ 4]*tp[0] + Ad[ 5]*tp[1] + Ad[ 6]*tp[2] + Ad[ 7]*tp[3])+
     coefs[i+2]*(Ad[ 8]*tp[0] + Ad[ 9]*tp[1] + Ad[10]*tp[2] + Ad[11]*tp[3])+
     coefs[i+3]*(Ad[12]*tp[0] + Ad[13]*tp[1] + Ad[14]*tp[2] + Ad[15]*tp[3]));
}

/* Value and gradient */
void eval_UBspline_1d_z_vg (UBspline_1d_z * spline, double x, 
                            sf_double_complex* val, 
                            sf_double_complex* grad) {
  x -= spline->x_grid.start;
  double u = x*spline->x_grid.delta_inv;
  double ipart, t;
  t = modf (u, &ipart);
  int i = (int) ipart;

  double tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  sf_double_complex* coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Ad[ 0]*tp[0] + Ad[ 1]*tp[1] + Ad[ 2]*tp[2] + Ad[ 3]*tp[3])+
     coefs[i+1]*(Ad[ 4]*tp[0] + Ad[ 5]*tp[1] + Ad[ 6]*tp[2] + Ad[ 7]*tp[3])+
     coefs[i+2]*(Ad[ 8]*tp[0] + Ad[ 9]*tp[1] + Ad[10]*tp[2] + Ad[11]*tp[3])+
     coefs[i+3]*(Ad[12]*tp[0] + Ad[13]*tp[1] + Ad[14]*tp[2] + Ad[15]*tp[3]));
  *grad = spline->x_grid.delta_inv * 
    (coefs[i+0]*(dAd[ 1]*tp[1] + dAd[ 2]*tp[2] + dAd[ 3]*tp[3])+
     coefs[i+1]*(dAd[ 5]*tp[1] + dAd[ 6]*tp[2] + dAd[ 7]*tp[3])+
     coefs[i+2]*(dAd[ 9]*tp[1] + dAd[10]*tp[2] + dAd[11]*tp[3])+
     coefs[i+3]*(dAd[13]*tp[1] + dAd[14]*tp[2] + dAd[15]*tp[3]));
}

/* Value, gradient, and laplacian */
void eval_UBspline_1d_z_vgl (UBspline_1d_z * spline, double x, 
                             sf_double_complex* val, sf_double_complex* grad,
                             sf_double_complex* lapl) {
  x -= spline->x_grid.start;
  double u = x*spline->x_grid.delta_inv;
  double ipart, t;
  t = modf (u, &ipart);
  int i = (int) ipart;

  double tp[4];
  tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
  sf_double_complex* coefs = spline->coefs;

  *val = 
    (coefs[i+0]*(Ad[ 0]*tp[0] + Ad[ 1]*tp[1] + Ad[ 2]*tp[2] + Ad[ 3]*tp[3])+
     coefs[i+1]*(Ad[ 4]*tp[0] + Ad[ 5]*tp[1] + Ad[ 6]*tp[2] + Ad[ 7]*tp[3])+
     coefs[i+2]*(Ad[ 8]*tp[0] + Ad[ 9]*tp[1] + Ad[10]*tp[2] + Ad[11]*tp[3])+
     coefs[i+3]*(Ad[12]*tp[0] + Ad[13]*tp[1] + Ad[14]*tp[2] + Ad[15]*tp[3]));
  *grad = spline->x_grid.delta_inv * 
    (coefs[i+0]*(dAd[ 1]*tp[1] + dAd[ 2]*tp[2] + dAd[ 3]*tp[3])+
     coefs[i+1]*(dAd[ 5]*tp[1] + dAd[ 6]*tp[2] + dAd[ 7]*tp[3])+
     coefs[i+2]*(dAd[ 9]*tp[1] + dAd[10]*tp[2] + dAd[11]*tp[3])+
     coefs[i+3]*(dAd[13]*tp[1] + dAd[14]*tp[2] + dAd[15]*tp[3]));
  *lapl = spline->x_grid.delta_inv * spline->x_grid.delta_inv * 
    (coefs[i+0]*(d2Ad[ 2]*tp[2] + d2Ad[ 3]*tp[3])+
     coefs[i+1]*(d2Ad[ 6]*tp[2] + d2Ad[ 7]*tp[3])+
     coefs[i+2]*(d2Ad[10]*tp[2] + d2Ad[11]*tp[3])+
     coefs[i+3]*(d2Ad[14]*tp[2] + d2Ad[15]*tp[3]));
}

/* Value, gradient, and Hessian */
void eval_UBspline_1d_z_vgh (UBspline_1d_z * spline, double x, 
                             sf_double_complex* val, 
                             sf_double_complex* grad,
                             sf_double_complex* hess) {
  eval_UBspline_1d_z_vgh (spline, x, val, grad, hess);
}

/************************************************************/
/* 2D double-precision, complex evaulation functions        */
/************************************************************/

/* Value only */
void eval_UBspline_2d_z (UBspline_2d_z * spline, 
                         double x, double y, sf_double_complex* val) {
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double ipartx, iparty, tx, ty;
  tx = modf (ux, &ipartx);
  ty = modf (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;

  double tpx[4], tpy[4], a[4], b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  sf_double_complex* coefs = spline->coefs;

  a[0] = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1] = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2] = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3] = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);

  b[0] = (Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1] = (Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2] = (Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3] = (Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);

  int xs = spline->x_stride;
#define C(i,j) coefs[(ix+(i))*xs+iy+(j)]
  *val = (a[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
	  a[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
	  a[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
	  a[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));
#undef C

}

/* Value and gradient */
void eval_UBspline_2d_z_vg (UBspline_2d_z * spline, 
                            double x, double y, 
                            sf_double_complex* val, 
                            sf_double_complex* grad) {
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double ipartx, iparty, tx, ty;
  tx = modf (ux, &ipartx);
  ty = modf (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;

  double tpx[4], tpy[4], a[4], b[4], da[4], db[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  sf_double_complex* coefs = spline->coefs;

  a[0]  = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]  = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]  = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]  = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  da[0] = (dAd[ 1]*tpx[1] + dAd[ 2]*tpx[2] + dAd[ 3]*tpx[3]);
  da[1] = (dAd[ 5]*tpx[1] + dAd[ 6]*tpx[2] + dAd[ 7]*tpx[3]);
  da[2] = (dAd[ 9]*tpx[1] + dAd[10]*tpx[2] + dAd[11]*tpx[3]);
  da[3] = (dAd[13]*tpx[1] + dAd[14]*tpx[2] + dAd[15]*tpx[3]);

  b[0]  = (Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1]  = (Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2]  = (Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3]  = (Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  db[0] = (dAd[ 1]*tpy[1] + dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
  db[1] = (dAd[ 5]*tpy[1] + dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
  db[2] = (dAd[ 9]*tpy[1] + dAd[10]*tpy[2] + dAd[11]*tpy[3]);
  db[3] = (dAd[13]*tpy[1] + dAd[14]*tpy[2] + dAd[15]*tpy[3]);

  int xs = spline->x_stride;
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
void eval_UBspline_2d_z_vgl (UBspline_2d_z * spline, 
                             double x, double y, sf_double_complex* val, 
                             sf_double_complex* grad, 
                             sf_double_complex* lapl) {
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double ipartx, iparty, tx, ty;
  tx = modf (ux, &ipartx);
  ty = modf (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;

  double tpx[4], tpy[4], a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  sf_double_complex* coefs = spline->coefs;

  a[0]   = (  Ad[ 0]*tpx[0] +   Ad[ 1]*tpx[1] +  Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]   = (  Ad[ 4]*tpx[0] +   Ad[ 5]*tpx[1] +  Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]   = (  Ad[ 8]*tpx[0] +   Ad[ 9]*tpx[1] +  Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]   = (  Ad[12]*tpx[0] +   Ad[13]*tpx[1] +  Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  da[0]  = ( dAd[ 1]*tpx[1] +  dAd[ 2]*tpx[2] + dAd[ 3]*tpx[3]);
  da[1]  = ( dAd[ 5]*tpx[1] +  dAd[ 6]*tpx[2] + dAd[ 7]*tpx[3]);
  da[2]  = ( dAd[ 9]*tpx[1] +  dAd[10]*tpx[2] + dAd[11]*tpx[3]);
  da[3]  = ( dAd[13]*tpx[1] +  dAd[14]*tpx[2] + dAd[15]*tpx[3]);
  d2a[0] = (d2Ad[ 2]*tpx[2] + d2Ad[ 3]*tpx[3]);
  d2a[1] = (d2Ad[ 6]*tpx[2] + d2Ad[ 7]*tpx[3]);
  d2a[2] = (d2Ad[10]*tpx[2] + d2Ad[11]*tpx[3]);
  d2a[3] = (d2Ad[14]*tpx[2] + d2Ad[15]*tpx[3]);

  b[0]  = ( Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1]  = ( Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2]  = ( Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3]  = ( Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  db[0] = (dAd[ 1]*tpy[1] + dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
  db[1] = (dAd[ 5]*tpy[1] + dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
  db[2] = (dAd[ 9]*tpy[1] + dAd[10]*tpy[2] + dAd[11]*tpy[3]);
  db[3] = (dAd[13]*tpy[1] + dAd[14]*tpy[2] + dAd[15]*tpy[3]);
  d2b[0] = (d2Ad[ 2]*tpy[2] + d2Ad[ 3]*tpy[3]);
  d2b[1] = (d2Ad[ 6]*tpy[2] + d2Ad[ 7]*tpy[3]);
  d2b[2] = (d2Ad[10]*tpy[2] + d2Ad[11]*tpy[3]);
  d2b[3] = (d2Ad[14]*tpy[2] + d2Ad[15]*tpy[3]);

  int xs = spline->x_stride;
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
void eval_UBspline_2d_z_vgh (UBspline_2d_z * spline, 
                             double x, double y, sf_double_complex* val, 
                             sf_double_complex* grad, 
                             sf_double_complex* hess) {
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double ipartx, iparty, tx, ty;
  tx = modf (ux, &ipartx);
  ty = modf (uy, &iparty);
  int ix = (int) ipartx;
  int iy = (int) iparty;

  double tpx[4], tpy[4], a[4], b[4], da[4], db[4], d2a[4], d2b[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  sf_double_complex* coefs = spline->coefs;

  a[0]   = (  Ad[ 0]*tpx[0] +   Ad[ 1]*tpx[1] +  Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]   = (  Ad[ 4]*tpx[0] +   Ad[ 5]*tpx[1] +  Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]   = (  Ad[ 8]*tpx[0] +   Ad[ 9]*tpx[1] +  Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]   = (  Ad[12]*tpx[0] +   Ad[13]*tpx[1] +  Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  da[0]  = ( dAd[ 1]*tpx[1] +  dAd[ 2]*tpx[2] + dAd[ 3]*tpx[3]);
  da[1]  = ( dAd[ 5]*tpx[1] +  dAd[ 6]*tpx[2] + dAd[ 7]*tpx[3]);
  da[2]  = ( dAd[ 9]*tpx[1] +  dAd[10]*tpx[2] + dAd[11]*tpx[3]);
  da[3]  = ( dAd[13]*tpx[1] +  dAd[14]*tpx[2] + dAd[15]*tpx[3]);
  d2a[0] = (d2Ad[ 2]*tpx[2] + d2Ad[ 3]*tpx[3]);
  d2a[1] = (d2Ad[ 6]*tpx[2] + d2Ad[ 7]*tpx[3]);
  d2a[2] = (d2Ad[10]*tpx[2] + d2Ad[11]*tpx[3]);
  d2a[3] = (d2Ad[14]*tpx[2] + d2Ad[15]*tpx[3]);

  b[0]   = (  Ad[ 0]*tpy[0] +   Ad[ 1]*tpy[1] +  Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1]   = (  Ad[ 4]*tpy[0] +   Ad[ 5]*tpy[1] +  Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2]   = (  Ad[ 8]*tpy[0] +   Ad[ 9]*tpy[1] +  Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3]   = (  Ad[12]*tpy[0] +   Ad[13]*tpy[1] +  Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  db[0]  = ( dAd[ 1]*tpy[1] +  dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
  db[1]  = ( dAd[ 5]*tpy[1] +  dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
  db[2]  = ( dAd[ 9]*tpy[1] +  dAd[10]*tpy[2] + dAd[11]*tpy[3]);
  db[3]  = ( dAd[13]*tpy[1] +  dAd[14]*tpy[2] + dAd[15]*tpy[3]);
  d2b[0] = (d2Ad[ 2]*tpy[2] + d2Ad[ 3]*tpy[3]);
  d2b[1] = (d2Ad[ 6]*tpy[2] + d2Ad[ 7]*tpy[3]);
  d2b[2] = (d2Ad[10]*tpy[2] + d2Ad[11]*tpy[3]);
  d2b[3] = (d2Ad[14]*tpy[2] + d2Ad[15]*tpy[3]);

  int xs = spline->x_stride;
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
/* 3D double-precision, complex evaulation functions           */
/************************************************************/

/* Value only */
void eval_UBspline_3d_z (UBspline_3d_z * spline, 
                         double x, double y, double z,
                         sf_double_complex* val) {
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  tz = modf (uz, &ipartz);  int iz = (int) ipartz;

  double tpx[4], tpy[4], tpz[4], a[4], b[4], c[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  sf_double_complex* coefs = spline->coefs;

  a[0] = (Ad[ 0]*tpx[0] + Ad[ 1]*tpx[1] + Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1] = (Ad[ 4]*tpx[0] + Ad[ 5]*tpx[1] + Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2] = (Ad[ 8]*tpx[0] + Ad[ 9]*tpx[1] + Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3] = (Ad[12]*tpx[0] + Ad[13]*tpx[1] + Ad[14]*tpx[2] + Ad[15]*tpx[3]);

  b[0] = (Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1] = (Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2] = (Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3] = (Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);

  c[0] = (Ad[ 0]*tpz[0] + Ad[ 1]*tpz[1] + Ad[ 2]*tpz[2] + Ad[ 3]*tpz[3]);
  c[1] = (Ad[ 4]*tpz[0] + Ad[ 5]*tpz[1] + Ad[ 6]*tpz[2] + Ad[ 7]*tpz[3]);
  c[2] = (Ad[ 8]*tpz[0] + Ad[ 9]*tpz[1] + Ad[10]*tpz[2] + Ad[11]*tpz[3]);
  c[3] = (Ad[12]*tpz[0] + Ad[13]*tpz[1] + Ad[14]*tpz[2] + Ad[15]*tpz[3]);

  int xs = spline->x_stride;
  int ys = spline->y_stride;
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
void eval_UBspline_3d_z_vg (UBspline_3d_z * spline, 
                            double x, double y, double z,
                            sf_double_complex* val, 
                            sf_double_complex* grad) {
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;  
  ty = modf (uy, &iparty);  int iy = (int) iparty; 
  tz = modf (uz, &ipartz);  int iz = (int) ipartz; 

  double tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], da[4], db[4], dc[4];
  sf_double_complex cP[16], bcP[4], dbcP[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  sf_double_complex* coefs = spline->coefs;

  a[0]   = (  Ad[ 0]*tpx[0] +   Ad[ 1]*tpx[1] +  Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]   = (  Ad[ 4]*tpx[0] +   Ad[ 5]*tpx[1] +  Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]   = (  Ad[ 8]*tpx[0] +   Ad[ 9]*tpx[1] +  Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]   = (  Ad[12]*tpx[0] +   Ad[13]*tpx[1] +  Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  da[0]  = ( dAd[ 1]*tpx[1] +  dAd[ 2]*tpx[2] + dAd[ 3]*tpx[3]);
  da[1]  = ( dAd[ 5]*tpx[1] +  dAd[ 6]*tpx[2] + dAd[ 7]*tpx[3]);
  da[2]  = ( dAd[ 9]*tpx[1] +  dAd[10]*tpx[2] + dAd[11]*tpx[3]);
  da[3]  = ( dAd[13]*tpx[1] +  dAd[14]*tpx[2] + dAd[15]*tpx[3]);

  b[0]  = ( Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1]  = ( Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2]  = ( Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3]  = ( Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  db[0] = (dAd[ 1]*tpy[1] + dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
  db[1] = (dAd[ 5]*tpy[1] + dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
  db[2] = (dAd[ 9]*tpy[1] + dAd[10]*tpy[2] + dAd[11]*tpy[3]);
  db[3] = (dAd[13]*tpy[1] + dAd[14]*tpy[2] + dAd[15]*tpy[3]);

  c[0]  = ( Ad[ 0]*tpz[0] + Ad[ 1]*tpz[1] + Ad[ 2]*tpz[2] + Ad[ 3]*tpz[3]);
  c[1]  = ( Ad[ 4]*tpz[0] + Ad[ 5]*tpz[1] + Ad[ 6]*tpz[2] + Ad[ 7]*tpz[3]);
  c[2]  = ( Ad[ 8]*tpz[0] + Ad[ 9]*tpz[1] + Ad[10]*tpz[2] + Ad[11]*tpz[3]);
  c[3]  = ( Ad[12]*tpz[0] + Ad[13]*tpz[1] + Ad[14]*tpz[2] + Ad[15]*tpz[3]);
  dc[0] = (dAd[ 1]*tpz[1] + dAd[ 2]*tpz[2] + dAd[ 3]*tpz[3]);
  dc[1] = (dAd[ 5]*tpz[1] + dAd[ 6]*tpz[2] + dAd[ 7]*tpz[3]);
  dc[2] = (dAd[ 9]*tpz[1] + dAd[10]*tpz[2] + dAd[11]*tpz[3]);
  dc[3] = (dAd[13]*tpz[1] + dAd[14]*tpz[2] + dAd[15]*tpz[3]);

  int xs = spline->x_stride;
  int ys = spline->y_stride;
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
void eval_UBspline_3d_z_vgl (UBspline_3d_z * spline, 
                             double x, double y, double z,
                             sf_double_complex* val, 
                             sf_double_complex* grad, 
                             sf_double_complex* lapl) {
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;  
  ty = modf (uy, &iparty);  int iy = (int) iparty; 
  tz = modf (uz, &ipartz);  int iz = (int) ipartz; 

  double tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], da[4], db[4], dc[4], 
    d2a[4], d2b[4], d2c[4];
  sf_double_complex cP[16], dcP[16], bcP[4], dbcP[4], d2bcP[4], bdcP[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  sf_double_complex* coefs = spline->coefs;

  a[0]   = (  Ad[ 0]*tpx[0] +   Ad[ 1]*tpx[1] +  Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]   = (  Ad[ 4]*tpx[0] +   Ad[ 5]*tpx[1] +  Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]   = (  Ad[ 8]*tpx[0] +   Ad[ 9]*tpx[1] +  Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]   = (  Ad[12]*tpx[0] +   Ad[13]*tpx[1] +  Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  da[0]  = ( dAd[ 1]*tpx[1] +  dAd[ 2]*tpx[2] + dAd[ 3]*tpx[3]);
  da[1]  = ( dAd[ 5]*tpx[1] +  dAd[ 6]*tpx[2] + dAd[ 7]*tpx[3]);
  da[2]  = ( dAd[ 9]*tpx[1] +  dAd[10]*tpx[2] + dAd[11]*tpx[3]);
  da[3]  = ( dAd[13]*tpx[1] +  dAd[14]*tpx[2] + dAd[15]*tpx[3]);
  d2a[0] = (d2Ad[ 2]*tpx[2] + d2Ad[ 3]*tpx[3]);
  d2a[1] = (d2Ad[ 6]*tpx[2] + d2Ad[ 7]*tpx[3]);
  d2a[2] = (d2Ad[10]*tpx[2] + d2Ad[11]*tpx[3]);
  d2a[3] = (d2Ad[14]*tpx[2] + d2Ad[15]*tpx[3]);

  b[0]  = ( Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1]  = ( Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2]  = ( Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3]  = ( Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  db[0] = (dAd[ 1]*tpy[1] + dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
  db[1] = (dAd[ 5]*tpy[1] + dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
  db[2] = (dAd[ 9]*tpy[1] + dAd[10]*tpy[2] + dAd[11]*tpy[3]);
  db[3] = (dAd[13]*tpy[1] + dAd[14]*tpy[2] + dAd[15]*tpy[3]);
  d2b[0] = (d2Ad[ 2]*tpy[2] + d2Ad[ 3]*tpy[3]);
  d2b[1] = (d2Ad[ 6]*tpy[2] + d2Ad[ 7]*tpy[3]);
  d2b[2] = (d2Ad[10]*tpy[2] + d2Ad[11]*tpy[3]);
  d2b[3] = (d2Ad[14]*tpy[2] + d2Ad[15]*tpy[3]);

  c[0]  = ( Ad[ 0]*tpz[0] + Ad[ 1]*tpz[1] + Ad[ 2]*tpz[2] + Ad[ 3]*tpz[3]);
  c[1]  = ( Ad[ 4]*tpz[0] + Ad[ 5]*tpz[1] + Ad[ 6]*tpz[2] + Ad[ 7]*tpz[3]);
  c[2]  = ( Ad[ 8]*tpz[0] + Ad[ 9]*tpz[1] + Ad[10]*tpz[2] + Ad[11]*tpz[3]);
  c[3]  = ( Ad[12]*tpz[0] + Ad[13]*tpz[1] + Ad[14]*tpz[2] + Ad[15]*tpz[3]);
  dc[0] = (dAd[ 1]*tpz[1] + dAd[ 2]*tpz[2] + dAd[ 3]*tpz[3]);
  dc[1] = (dAd[ 5]*tpz[1] + dAd[ 6]*tpz[2] + dAd[ 7]*tpz[3]);
  dc[2] = (dAd[ 9]*tpz[1] + dAd[10]*tpz[2] + dAd[11]*tpz[3]);
  dc[3] = (dAd[13]*tpz[1] + dAd[14]*tpz[2] + dAd[15]*tpz[3]);
  d2c[0] = (d2Ad[ 2]*tpz[2] + d2Ad[ 3]*tpz[3]);
  d2c[1] = (d2Ad[ 6]*tpz[2] + d2Ad[ 7]*tpz[3]);
  d2c[2] = (d2Ad[10]*tpz[2] + d2Ad[11]*tpz[3]);
  d2c[3] = (d2Ad[14]*tpz[2] + d2Ad[15]*tpz[3]);

  int xs = spline->x_stride;
  int ys = spline->y_stride;
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
void eval_UBspline_3d_z_vgh (UBspline_3d_z * spline, 
                             double x, double y, double z,
                             sf_double_complex* val, 
                             sf_double_complex* grad, 
                             sf_double_complex* hess) {
  x -= spline->x_grid.start;
  y -= spline->y_grid.start;  
  z -= spline->z_grid.start;
  double ux = x*spline->x_grid.delta_inv;
  double uy = y*spline->y_grid.delta_inv;
  double uz = z*spline->z_grid.delta_inv;
  ux = fmin (ux, (double)(spline->x_grid.num)-1.0e-5);
  uy = fmin (uy, (double)(spline->y_grid.num)-1.0e-5);
  uz = fmin (uz, (double)(spline->z_grid.num)-1.0e-5);
  double ipartx, iparty, ipartz, tx, ty, tz;
  tx = modf (ux, &ipartx);  int ix = (int) ipartx;
  ty = modf (uy, &iparty);  int iy = (int) iparty;
  tz = modf (uz, &ipartz);  int iz = (int) ipartz;

  double tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], da[4], db[4], dc[4], 
    d2a[4], d2b[4], d2c[4];
  sf_double_complex cP[16], dcP[16], d2cP[16], bcP[4], dbcP[4],
    d2bcP[4], dbdcP[4], bd2cP[4], bdcP[4];
  tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
  tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
  tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;
  sf_double_complex* coefs = spline->coefs;

  a[0]   = (  Ad[ 0]*tpx[0] +   Ad[ 1]*tpx[1] +  Ad[ 2]*tpx[2] + Ad[ 3]*tpx[3]);
  a[1]   = (  Ad[ 4]*tpx[0] +   Ad[ 5]*tpx[1] +  Ad[ 6]*tpx[2] + Ad[ 7]*tpx[3]);
  a[2]   = (  Ad[ 8]*tpx[0] +   Ad[ 9]*tpx[1] +  Ad[10]*tpx[2] + Ad[11]*tpx[3]);
  a[3]   = (  Ad[12]*tpx[0] +   Ad[13]*tpx[1] +  Ad[14]*tpx[2] + Ad[15]*tpx[3]);
  da[0]  = ( dAd[ 1]*tpx[1] +  dAd[ 2]*tpx[2] + dAd[ 3]*tpx[3]);
  da[1]  = ( dAd[ 5]*tpx[1] +  dAd[ 6]*tpx[2] + dAd[ 7]*tpx[3]);
  da[2]  = ( dAd[ 9]*tpx[1] +  dAd[10]*tpx[2] + dAd[11]*tpx[3]);
  da[3]  = ( dAd[13]*tpx[1] +  dAd[14]*tpx[2] + dAd[15]*tpx[3]);
  d2a[0] = (d2Ad[ 2]*tpx[2] + d2Ad[ 3]*tpx[3]);
  d2a[1] = (d2Ad[ 6]*tpx[2] + d2Ad[ 7]*tpx[3]);
  d2a[2] = (d2Ad[10]*tpx[2] + d2Ad[11]*tpx[3]);
  d2a[3] = (d2Ad[14]*tpx[2] + d2Ad[15]*tpx[3]);

  b[0]  = ( Ad[ 0]*tpy[0] + Ad[ 1]*tpy[1] + Ad[ 2]*tpy[2] + Ad[ 3]*tpy[3]);
  b[1]  = ( Ad[ 4]*tpy[0] + Ad[ 5]*tpy[1] + Ad[ 6]*tpy[2] + Ad[ 7]*tpy[3]);
  b[2]  = ( Ad[ 8]*tpy[0] + Ad[ 9]*tpy[1] + Ad[10]*tpy[2] + Ad[11]*tpy[3]);
  b[3]  = ( Ad[12]*tpy[0] + Ad[13]*tpy[1] + Ad[14]*tpy[2] + Ad[15]*tpy[3]);
  db[0] = (dAd[ 1]*tpy[1] + dAd[ 2]*tpy[2] + dAd[ 3]*tpy[3]);
  db[1] = (dAd[ 5]*tpy[1] + dAd[ 6]*tpy[2] + dAd[ 7]*tpy[3]);
  db[2] = (dAd[ 9]*tpy[1] + dAd[10]*tpy[2] + dAd[11]*tpy[3]);
  db[3] = (dAd[13]*tpy[1] + dAd[14]*tpy[2] + dAd[15]*tpy[3]);
  d2b[0] = (d2Ad[ 2]*tpy[2] + d2Ad[ 3]*tpy[3]);
  d2b[1] = (d2Ad[ 6]*tpy[2] + d2Ad[ 7]*tpy[3]);
  d2b[2] = (d2Ad[10]*tpy[2] + d2Ad[11]*tpy[3]);
  d2b[3] = (d2Ad[14]*tpy[2] + d2Ad[15]*tpy[3]);

  c[0]  = ( Ad[ 0]*tpz[0] + Ad[ 1]*tpz[1] + Ad[ 2]*tpz[2] + Ad[ 3]*tpz[3]);
  c[1]  = ( Ad[ 4]*tpz[0] + Ad[ 5]*tpz[1] + Ad[ 6]*tpz[2] + Ad[ 7]*tpz[3]);
  c[2]  = ( Ad[ 8]*tpz[0] + Ad[ 9]*tpz[1] + Ad[10]*tpz[2] + Ad[11]*tpz[3]);
  c[3]  = ( Ad[12]*tpz[0] + Ad[13]*tpz[1] + Ad[14]*tpz[2] + Ad[15]*tpz[3]);
  dc[0] = (dAd[ 1]*tpz[1] + dAd[ 2]*tpz[2] + dAd[ 3]*tpz[3]);
  dc[1] = (dAd[ 5]*tpz[1] + dAd[ 6]*tpz[2] + dAd[ 7]*tpz[3]);
  dc[2] = (dAd[ 9]*tpz[1] + dAd[10]*tpz[2] + dAd[11]*tpz[3]);
  dc[3] = (dAd[13]*tpz[1] + dAd[14]*tpz[2] + dAd[15]*tpz[3]);
  d2c[0] = (d2Ad[ 2]*tpz[2] + d2Ad[ 3]*tpz[3]);
  d2c[1] = (d2Ad[ 6]*tpz[2] + d2Ad[ 7]*tpz[3]);
  d2c[2] = (d2Ad[10]*tpz[2] + d2Ad[11]*tpz[3]);
  d2c[3] = (d2Ad[14]*tpz[2] + d2Ad[15]*tpz[3]);

  int xs = spline->x_stride;
  int ys = spline->y_stride;

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
#endif /* NO_COMPLEX */
