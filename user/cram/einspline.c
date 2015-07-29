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

#ifdef HAVE_SSE
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>

__m128 *restrict A_s = (__m128 *)0;

static void init_sse_data (void)
{
    if (A_s == 0) {
	posix_memalign ((void**)&A_s, 16, (sizeof(__m128)*12));
	A_s[0]  = _mm_setr_ps ( 1.0/6.0, -3.0/6.0,  3.0/6.0, -1.0/6.0 );
	A_s[0]  = _mm_setr_ps ( 1.0/6.0, -3.0/6.0,  3.0/6.0, -1.0/6.0 );          
	A_s[1]  = _mm_setr_ps ( 4.0/6.0,  0.0/6.0, -6.0/6.0,  3.0/6.0 );          
	A_s[2]  = _mm_setr_ps ( 1.0/6.0,  3.0/6.0,  3.0/6.0, -3.0/6.0 );          
	A_s[3]  = _mm_setr_ps ( 0.0/6.0,  0.0/6.0,  0.0/6.0,  1.0/6.0 );          
	A_s[4]  = _mm_setr_ps ( -0.5,  1.0, -0.5, 0.0  );                  
	A_s[5]  = _mm_setr_ps (  0.0, -2.0,  1.5, 0.0  );                  
	A_s[6]  = _mm_setr_ps (  0.5,  1.0, -1.5, 0.0  );                  
	A_s[7]  = _mm_setr_ps (  0.0,  0.0,  0.5, 0.0  );                  
	A_s[8]  = _mm_setr_ps (  1.0, -1.0,  0.0, 0.0  );                  
	A_s[9]  = _mm_setr_ps ( -2.0,  3.0,  0.0, 0.0  );                  
	A_s[10] = _mm_setr_ps (  1.0, -3.0,  0.0, 0.0  );                  
	A_s[11] = _mm_setr_ps (  0.0,  1.0,  0.0, 0.0  );                  
    }
}

#define _MM_MATVEC4_PS(M0, M1, M2, M3, v, r)				\
    do {                                                                \
	__m128 r0 = _mm_hadd_ps (_mm_mul_ps (M0, v), _mm_mul_ps (M1, v)); \
	__m128 r1 = _mm_hadd_ps (_mm_mul_ps (M2, v), _mm_mul_ps (M3, v)); \
	r = _mm_hadd_ps (r0, r1);					\
    } while (0);
#define _MM_DOT4_PS(A, B, _p)			\
    do {					\
	__m128 t  = _mm_mul_ps (A, B);		\
	__m128 t1 = _mm_hadd_ps (t,t);		\
	__m128 r  = _mm_hadd_ps (t1, t1);	\
	_mm_store_ss (&(_p), r);		\
    } while(0);

#endif

void init_einspline (void) {
#ifdef HAVE_SSE
    init_sse_data ();
#endif
}

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
                                     intptr_t M, intptr_t cstride) {
    intptr_t row;
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
                                        intptr_t M, intptr_t cstride) {
    float *lastCol;
    intptr_t row;

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
                                            intptr_t M, intptr_t cstride)
{
    intptr_t row;
    float *lastCol;

    lastCol = (float*) malloc(M*sizeof(float));

    bands[4*0+0]     *= -1.0;
    bands[4*(M-1)+2] *= -1.0;

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

    free(lastCol);
}

static void find_coefs_1d_s (Ugrid grid, BCtype_s bc, 
                             float *data,  intptr_t dstride,
                             float *coefs, intptr_t cstride) {
    intptr_t M = grid.num, i, j;
    float basis[4] = {1.0/6.0, 2.0/3.0, 1.0/6.0, 0.0};
    float *bands;
    float abcd_left[4], abcd_right[4];

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
    intptr_t M, N;
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
#ifndef HAVE_SSE
    spline->coefs = malloc ((intptr_t)sizeof(float)*(intptr_t)N);
#else
    posix_memalign ((void**)&spline->coefs, 16, (sizeof(float)*N));
#endif
    spline->nc = (intptr_t)sizeof(float)*(intptr_t)N;
    find_coefs_1d_s (spline->x_grid, xBC, data, 1, spline->coefs, 1);
#ifdef HAVE_SSE
    init_sse_data ();
#endif
    return spline;
}

void recompute_UBspline_1d_s (UBspline_1d_s* spline, float *data) {
    find_coefs_1d_s (spline->x_grid, spline->xBC, data, 1, spline->coefs, 1);
}

UBspline_2d_s* create_UBspline_2d_s (Ugrid x_grid, Ugrid y_grid,
                                     BCtype_s xBC, BCtype_s yBC, float *data) {
    /* Create new spline */
    intptr_t Mx, My;
    intptr_t Nx, Ny, iy, ix;
    intptr_t doffset;
    intptr_t coffset;

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
#ifndef HAVE_SSE
    spline->coefs = malloc ((intptr_t)sizeof(float)*(intptr_t)Nx*(intptr_t)Ny);
#else
    posix_memalign ((void**)&spline->coefs, 16, sizeof(float)*Nx*Ny);
#endif
    spline->nc = (intptr_t)sizeof(float)*(intptr_t)Nx*(intptr_t)Ny;

    /* First, solve in the X-direction */
    for (iy=0; iy<My; iy++) {
	doffset = iy;
	coffset = iy;
	find_coefs_1d_s (spline->x_grid, spline->xBC, data+doffset, My,
			 spline->coefs+coffset, Ny);
    }

    /* Now, solve in the Y-direction */
    for (ix=0; ix<Nx; ix++) {
	doffset = ix*Ny;
	coffset = ix*Ny;
	find_coefs_1d_s (spline->y_grid, spline->yBC, spline->coefs+doffset, 1, 
			 spline->coefs+coffset, 1);
    }
#ifdef HAVE_SSE
    init_sse_data ();
#endif
    return spline;
}

void recompute_UBspline_2d_s (UBspline_2d_s* spline, float *data) {
    intptr_t Mx = spline->x_grid.num;
    intptr_t My = spline->y_grid.num;
    intptr_t Nx, Ny, iy, ix;
    intptr_t doffset;
    intptr_t coffset;

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
	doffset = iy;
	coffset = iy;
	find_coefs_1d_s (spline->x_grid, spline->xBC, data+doffset, My,
			 spline->coefs+coffset, Ny);
    }

    /* Now, solve in the Y-direction */
    for (ix=0; ix<Nx; ix++) {
	doffset = ix*Ny;
	coffset = ix*Ny;
	find_coefs_1d_s (spline->y_grid, spline->yBC, spline->coefs+doffset, 1, 
			 spline->coefs+coffset, 1);
    }
}

UBspline_3d_s* create_UBspline_3d_s (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
                                     BCtype_s xBC, BCtype_s yBC, BCtype_s zBC,
                                     float *data) {
    /* Create new spline */
    intptr_t Mx, My, Mz;
    intptr_t Nx, Ny, Nz, iy, ix, iz;
    intptr_t doffset;
    intptr_t coffset;
    
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
#ifndef HAVE_SSE
    spline->coefs = malloc ((intptr_t)sizeof(float)*(intptr_t)Nx*(intptr_t)Ny*(intptr_t)Nz);
#else
    posix_memalign ((void**)&spline->coefs, 16, (sizeof(float)*Nx*Ny*Nz));
#endif
    spline->nc = (intptr_t)sizeof(float)*(intptr_t)Nx*(intptr_t)Ny*(intptr_t)Nz;

    /* First, solve in the X-direction */
    for (iy=0; iy<My; iy++) 
	for (iz=0; iz<Mz; iz++) {
	    doffset = iy*Mz+iz;
	    coffset = iy*Nz+iz;
	    find_coefs_1d_s (spline->x_grid, xBC, data+doffset, My*Mz,
			     spline->coefs+coffset, Ny*Nz);
	}

    /* Now, solve in the Y-direction */
    for (ix=0; ix<Nx; ix++) 
	for (iz=0; iz<Nz; iz++) {
	    doffset = ix*Ny*Nz + iz;
	    coffset = ix*Ny*Nz + iz;
	    find_coefs_1d_s (spline->y_grid, yBC, spline->coefs+doffset, Nz, 
			     spline->coefs+coffset, Nz);
	}

    /* Now, solve in the Z-direction */
    for (ix=0; ix<Nx; ix++) 
	for (iy=0; iy<Ny; iy++) {
	    doffset = (ix*Ny+iy)*Nz;
	    coffset = (ix*Ny+iy)*Nz;
	    find_coefs_1d_s (spline->z_grid, zBC, spline->coefs+doffset, 1, 
			     spline->coefs+coffset, 1);
	}
#ifdef HAVE_SSE
    init_sse_data ();
#endif
    return spline;
}

void recompute_UBspline_3d_s (UBspline_3d_s* spline, float *data) {
    intptr_t Mx = spline->x_grid.num;
    intptr_t My = spline->y_grid.num;
    intptr_t Mz = spline->z_grid.num;
    intptr_t Nx, Ny, Nz, ix, iy, iz;
    intptr_t doffset;
    intptr_t coffset;

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
	    doffset = iy*Mz+iz;
	    coffset = iy*Nz+iz;
	    find_coefs_1d_s (spline->x_grid, spline->xBC, data+doffset, My*Mz,
			     spline->coefs+coffset, Ny*Nz);
	}

    /* Now, solve in the Y-direction */
    for (ix=0; ix<Nx; ix++) 
	for (iz=0; iz<Nz; iz++) {
	    doffset = ix*Ny*Nz + iz;
	    coffset = ix*Ny*Nz + iz;
	    find_coefs_1d_s (spline->y_grid, spline->yBC, spline->coefs+doffset, Nz, 
			     spline->coefs+coffset, Nz);
	}

    /* Now, solve in the Z-direction */
    for (ix=0; ix<Nx; ix++) 
	for (iy=0; iy<Ny; iy++) {
	    doffset = (ix*Ny+iy)*Nz;
	    coffset = (ix*Ny+iy)*Nz;
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
    float u;
    float ipart, t;
    intptr_t i;
    float tp[4];
    float* restrict coefs = spline->coefs;

    x -= spline->x_grid.start;
    u = x*spline->x_grid.delta_inv;

    t = modff (u, &ipart);
    i = (intptr_t) ipart;
  
    tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;
 
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
    float u;
    float ipart, t;
    intptr_t i;
    float tp[4];
    float* restrict coefs = spline->coefs;

    x -= spline->x_grid.start;
    u = x*spline->x_grid.delta_inv;

    t = modff (u, &ipart);
    i = (intptr_t) ipart;
  
    tp[0] = t*t*t;  tp[1] = t*t;  tp[2] = t;  tp[3] = 1.0;

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

/************************************************************/
/* 2D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
#ifndef HAVE_SSE
void eval_UBspline_2d_s (UBspline_2d_s * restrict spline, 
                         double x, double y, float* restrict val)
{
    intptr_t xs;
    float ipartx, iparty, tx, ty, ux, uy;    
    intptr_t ix, iy;
    float tpx[4], tpy[4], a[4], b[4];
    float* restrict coefs = spline->coefs;

    x -= spline->x_grid.start;
    y -= spline->y_grid.start;
    ux = x*spline->x_grid.delta_inv;
    uy = y*spline->y_grid.delta_inv;

    tx = modff (ux, &ipartx);
    ty = modff (uy, &iparty);
    ix = (intptr_t) ipartx;
    iy = (intptr_t) iparty;
  

    tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
    tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;

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
#else
void eval_UBspline_2d_s (UBspline_2d_s * restrict spline, 
                         double x, double y, float* restrict val)
{
    _mm_prefetch ((const char*)  &A_s[0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[1],_MM_HINT_T0);  
    _mm_prefetch ((const char*)  &A_s[2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[3],_MM_HINT_T0);
    /* SSE mesh point determination */
    __m128 xy        = _mm_set_ps (x, y, 0.0, 0.0);
    __m128 x0y0      = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 0.0, 0.0);
    __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 0.0, 0.0);
    xy = _mm_sub_ps (xy, x0y0);
    /* ux = (x - x0)/delta_x and same for y */
    __m128 uxuy    = _mm_mul_ps (xy, delta_inv);
    /* intpart = trunc (ux, uy) */
    __m128i intpart  = _mm_cvttps_epi32(uxuy);
    __m128i ixiy;
    _mm_storeu_si128 (&ixiy, intpart);
    /* Store to memory for use in C expressions
       xmm registers are stored to memory in reverse order */
    int ix = ((int *)&ixiy)[3];
    int iy = ((int *)&ixiy)[2];

    int xs = spline->x_stride;
    /* This macro is used to give the pointer to coefficient data.
       i and j should be in the range [0,3].  Coefficients are read four
       at a time, so no j value is needed. */
#define P(i) (spline->coefs+(ix+(i))*xs+iy)
    /* Prefetch the data from main memory into cache so it's available
       when we need to use it. */
    _mm_prefetch ((const char*)P(0), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(1), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(2), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(3), _MM_HINT_T0);

    /* Now compute the vectors:
       tpx = [t_x^3 t_x^2 t_x 1]
       tpy = [t_y^3 t_y^2 t_y 1]
       tpz = [t_z^3 t_z^2 t_z 1] */
    __m128 ipart  = _mm_cvtepi32_ps (intpart);
    __m128 txty   = _mm_sub_ps (uxuy, ipart);
    __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
    __m128 t2     = _mm_mul_ps (txty, txty);
    __m128 t3     = _mm_mul_ps (t2, txty);
    __m128 tpx    = t3;
    __m128 tpy    = t2;
    __m128 tpz    = txty;
    __m128 zero   = one;
    _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

    /* a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
       da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
       A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3] */
    __m128 a, b, bP, tmp0, tmp1, tmp2, tmp3;
    /* x-dependent vectors */
    _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpx,   a);
    /* y-dependent vectors */
    _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpy,   b);
    /* Compute cP, dcP, and d2cP products 1/4 at a time to maximize
       register reuse and avoid rerereading from memory or cache.
       1st quarter */
    tmp0 = _mm_loadu_ps (P(0));
    tmp1 = _mm_loadu_ps (P(1));
    tmp2 = _mm_loadu_ps (P(2));
    tmp3 = _mm_loadu_ps (P(3));
    _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   b,   bP);
    /* Compute value */
    _MM_DOT4_PS (a, bP, *val);
#undef P
}
#endif

/* Value and gradient */
#ifndef HAVE_SSE
void eval_UBspline_2d_s_vg (UBspline_2d_s * restrict spline, 
                            double x, double y, 
                            float* restrict val, float* restrict grad)
{
    intptr_t xs;
    float ipartx, iparty, tx, ty, ux, uy;
    intptr_t ix, iy;
    float tpx[4], tpy[4], a[4], b[4], da[4], db[4];
    float* restrict coefs = spline->coefs;

    x -= spline->x_grid.start;
    y -= spline->y_grid.start;
    ux = x*spline->x_grid.delta_inv;
    uy = y*spline->y_grid.delta_inv;

    tx = modff (ux, &ipartx);
    ty = modff (uy, &iparty);
    ix = (intptr_t) ipartx;
    iy = (intptr_t) iparty;
  
    tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
    tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;

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
#else
void eval_UBspline_2d_s_vg (UBspline_2d_s * restrict spline, 
                            double x, double y, 
                            float* restrict val, float* restrict grad)
{
    _mm_prefetch ((const char*)  &A_s[0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[1],_MM_HINT_T0);  
    _mm_prefetch ((const char*)  &A_s[2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[3],_MM_HINT_T0);
    _mm_prefetch ((const char*)  &A_s[4],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[5],_MM_HINT_T0);  
    _mm_prefetch ((const char*)  &A_s[6],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[7],_MM_HINT_T0);
    /* SSE mesh point determination */
    __m128 xy        = _mm_set_ps (x, y, 0.0, 0.0);
    __m128 x0y0      = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 0.0, 0.0);
    __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 0.0, 0.0);
    xy = _mm_sub_ps (xy, x0y0);
    /* ux = (x - x0)/delta_x and same for y */
    __m128 uxuy    = _mm_mul_ps (xy, delta_inv);
    /* intpart = trunc (ux, uy) */
    __m128i intpart  = _mm_cvttps_epi32(uxuy);
    __m128i ixiy;
    _mm_storeu_si128 (&ixiy, intpart);
    /* Store to memory for use in C expressions
       xmm registers are stored to memory in reverse order */
    int ix = ((int *)&ixiy)[3];
    int iy = ((int *)&ixiy)[2];

    int xs = spline->x_stride;
    /* This macro is used to give the pointer to coefficient data.
       i and j should be in the range [0,3].  Coefficients are read four
       at a time, so no j value is needed. */
#define P(i) (spline->coefs+(ix+(i))*xs+iy)
    /* Prefetch the data from main memory into cache so it's available
       when we need to use it. */
    _mm_prefetch ((const char*)P(0), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(1), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(2), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(3), _MM_HINT_T0);

    /* Now compute the vectors:
       tpx = [t_x^3 t_x^2 t_x 1]
       tpy = [t_y^3 t_y^2 t_y 1]
       tpz = [t_z^3 t_z^2 t_z 1] */
    __m128 ipart  = _mm_cvtepi32_ps (intpart);
    __m128 txty   = _mm_sub_ps (uxuy, ipart);
    __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
    __m128 t2     = _mm_mul_ps (txty, txty);
    __m128 t3     = _mm_mul_ps (t2, txty);
    __m128 tpx    = t3;
    __m128 tpy    = t2;
    __m128 tpz    = txty;
    __m128 zero   = one;
    _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

    /* a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
       da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
       A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3] */
    __m128 a, b, da, db, bP, dbP, tmp0, tmp1, tmp2, tmp3;
    /* x-dependent vectors */
    _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpx,   a);
    _MM_MATVEC4_PS (A_s[4], A_s[5], A_s[6], A_s[7], tpx,  da);
    /* y-dependent vectors */
    _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpy,   b);
    _MM_MATVEC4_PS (A_s[4], A_s[5], A_s[6], A_s[7], tpy,  db);
    /* Compute cP, dcP, and d2cP products 1/4 at a time to maximize
       register reuse and avoid rerereading from memory or cache.
       1st quarter */
    tmp0 = _mm_loadu_ps (P(0));
    tmp1 = _mm_loadu_ps (P(1));
    tmp2 = _mm_loadu_ps (P(2));
    tmp3 = _mm_loadu_ps (P(3));
    _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   b,   bP);
    _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  db,  dbP);
    /* Compute value */
    _MM_DOT4_PS (a, bP, *val);
    /* Compute gradient */
    _MM_DOT4_PS (da, bP, grad[0]);
    _MM_DOT4_PS (a, dbP, grad[1]);
    /* Multiply gradients and hessians by appropriate grid inverses */
    float dxInv = spline->x_grid.delta_inv;
    float dyInv = spline->y_grid.delta_inv;
    grad[0] *= dxInv;
    grad[1] *= dyInv;
#undef P
}
#endif

/************************************************************/
/* 3D single-precision, real evaulation functions           */
/************************************************************/

/* Value only */
#ifndef HAVE_SSE
void eval_UBspline_3d_s (UBspline_3d_s * restrict spline, 
                         double x, double y, double z,
                         float* restrict val)
{
    intptr_t xs, ys;
    float ipartx, iparty, ipartz, tx, ty, tz, ux, uy, uz;
    intptr_t ix, iy, iz;
    float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4];
    float* restrict coefs = spline->coefs;

    x -= spline->x_grid.start;
    y -= spline->y_grid.start;
    z -= spline->z_grid.start;
    ux = x*spline->x_grid.delta_inv;
    uy = y*spline->y_grid.delta_inv;
    uz = z*spline->z_grid.delta_inv;

    tx = modff (ux, &ipartx);  ix = (intptr_t) ipartx;
    ty = modff (uy, &iparty);  iy = (intptr_t) iparty;
    tz = modff (uz, &ipartz);  iz = (intptr_t) ipartz;

    tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
    tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
    tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;

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
#else
void eval_UBspline_3d_s (UBspline_3d_s * restrict spline, 
                         double x, double y, double z,
                         float* restrict val)
{
    _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);  
    _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);

    /* SSE mesh point determination */
    __m128 xyz       = _mm_set_ps (x, y, z, 0.0);
    __m128 x0y0z0    = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 
				   spline->z_grid.start, 0.0);
    __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 
				   spline->z_grid.delta_inv, 0.0);
    xyz = _mm_sub_ps (xyz, x0y0z0);
    /* ux = (x - x0)/delta_x and same for y and z */
    __m128 uxuyuz    = _mm_mul_ps (xyz, delta_inv);
    /* intpart = trunc (ux, uy, uz) */
    __m128i intpart  = _mm_cvttps_epi32(uxuyuz);
    __m128i ixiyiz;
    _mm_storeu_si128 (&ixiyiz, intpart);
    /* Store to memory for use in C expressions
       xmm registers are stored to memory in reverse order */
    int ix = ((int *)&ixiyiz)[3];
    int iy = ((int *)&ixiyiz)[2];
    int iz = ((int *)&ixiyiz)[1];

    int xs = spline->x_stride;
    int ys = spline->y_stride;

    /* This macro is used to give the pointer to coefficient data.
       i and j should be in the range [0,3].  Coefficients are read four
       at a time, so no k value is needed. */
#define P(i,j) (spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz))
    /* Prefetch the data from main memory into cache so it's available
       when we need to use it. */
    _mm_prefetch ((const char*)P(0,0), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(0,1), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(0,2), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(0,3), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(1,0), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(1,1), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(1,2), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(1,3), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(2,0), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(2,1), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(2,2), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(2,3), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(3,0), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(3,1), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(3,2), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(3,3), _MM_HINT_T0);

    /* Now compute the vectors:
       tpx = [t_x^3 t_x^2 t_x 1]
       tpy = [t_y^3 t_y^2 t_y 1]
       tpz = [t_z^3 t_z^2 t_z 1] */
    __m128 ipart  = _mm_cvtepi32_ps (intpart);
    __m128 txtytz = _mm_sub_ps (uxuyuz, ipart);
    __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
    __m128 t2     = _mm_mul_ps (txtytz, txtytz);
    __m128 t3     = _mm_mul_ps (t2, txtytz);
    __m128 tpx    = t3;
    __m128 tpy    = t2;
    __m128 tpz    = txtytz;
    __m128 zero   = one;
    _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

    /* a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
       da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
       A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3] */
    __m128 a, b, c, cP[4],bcP,
	tmp0, tmp1, tmp2, tmp3;

    /* x-dependent vectors */
    _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpx,   a);
    /* y-dependent vectors */
    _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpy,   b);
    /* z-dependent vectors */
    _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpz,   c);

    /* Compute cP, dcP, and d2cP products 1/4 at a time to maximize
       register reuse and avoid rerereading from memory or cache.
       1st quarter */
    tmp0 = _mm_loadu_ps (P(0,0));  tmp1 = _mm_loadu_ps (P(0,1));
    tmp2 = _mm_loadu_ps (P(0,2));  tmp3 = _mm_loadu_ps (P(0,3));
    _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[0]);
    /* 2nd quarter */
    tmp0 = _mm_loadu_ps (P(1,0));  tmp1 = _mm_loadu_ps (P(1,1));
    tmp2 = _mm_loadu_ps (P(1,2));  tmp3 = _mm_loadu_ps (P(1,3));
    _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[1]);
    /* 3rd quarter */
    tmp0 = _mm_loadu_ps (P(2,0));  tmp1 = _mm_loadu_ps (P(2,1));
    tmp2 = _mm_loadu_ps (P(2,2));  tmp3 = _mm_loadu_ps (P(2,3));
    _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[2]);
    /* 4th quarter */
    tmp0 = _mm_loadu_ps (P(3,0));  tmp1 = _mm_loadu_ps (P(3,1));
    tmp2 = _mm_loadu_ps (P(3,2));  tmp3 = _mm_loadu_ps (P(3,3));
    _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[3]);
  
    /* Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products */
    _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3],   b,   bcP);

    /* Compute value */
    _MM_DOT4_PS (a, bcP, *val);

#undef P
}
#endif

/* Value and gradient */
#ifndef HAVE_SSE
void eval_UBspline_3d_s_vg (UBspline_3d_s * restrict spline, 
                            double x, double y, double z,
                            float* restrict val, float* restrict grad)
{
    intptr_t xs, ys;
    float ipartx, iparty, ipartz, tx, ty, tz, ux, uy, uz;
    intptr_t ix, iy, iz;
    float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], da[4], db[4], dc[4], 
	cP[16], bcP[4], dbcP[4];
    float* restrict coefs = spline->coefs;

    x -= spline->x_grid.start;
    y -= spline->y_grid.start;
    z -= spline->z_grid.start;
    ux = x*spline->x_grid.delta_inv;
    uy = y*spline->y_grid.delta_inv;
    uz = z*spline->z_grid.delta_inv;
    tx = modff (ux, &ipartx);  ix = (intptr_t) ipartx;  
    ty = modff (uy, &iparty);  iy = (intptr_t) iparty; 
    tz = modff (uz, &ipartz);  iz = (intptr_t) ipartz; 
  
    tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
    tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
    tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;

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
#else
/* Value and gradient */
void eval_UBspline_3d_s_vg (UBspline_3d_s * restrict spline, 
                            double x, double y, double z,
                            float* restrict val, float* restrict grad)
{
    _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);  
    _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
    _mm_prefetch ((const char*)  &A_s[ 4],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 5],_MM_HINT_T0);  
    _mm_prefetch ((const char*)  &A_s[ 6],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 7],_MM_HINT_T0);

    /* SSE mesh point determination */
    __m128 xyz       = _mm_set_ps (x, y, z, 0.0);
    __m128 x0y0z0    = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 
				   spline->z_grid.start, 0.0);
    __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 
				   spline->z_grid.delta_inv, 0.0);
    xyz = _mm_sub_ps (xyz, x0y0z0);
    /* ux = (x - x0)/delta_x and same for y and z */
    __m128 uxuyuz    = _mm_mul_ps (xyz, delta_inv);
    /* intpart = trunc (ux, uy, uz) */
    __m128i intpart  = _mm_cvttps_epi32(uxuyuz);
    __m128i ixiyiz;
    _mm_storeu_si128 (&ixiyiz, intpart);
    /* Store to memory for use in C expressions
       xmm registers are stored to memory in reverse order */
    int ix = ((int *)&ixiyiz)[3];
    int iy = ((int *)&ixiyiz)[2];
    int iz = ((int *)&ixiyiz)[1];

    int xs = spline->x_stride;
    int ys = spline->y_stride;

    /* This macro is used to give the pointer to coefficient data.
       i and j should be in the range [0,3].  Coefficients are read four
       at a time, so no k value is needed. */
#define P(i,j) (spline->coefs+(ix+(i))*xs+(iy+(j))*ys+(iz))
    /* Prefetch the data from main memory into cache so it's available
       when we need to use it. */
    _mm_prefetch ((const char*)P(0,0), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(0,1), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(0,2), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(0,3), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(1,0), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(1,1), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(1,2), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(1,3), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(2,0), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(2,1), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(2,2), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(2,3), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(3,0), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(3,1), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(3,2), _MM_HINT_T0);
    _mm_prefetch ((const char*)P(3,3), _MM_HINT_T0);

    /* Now compute the vectors:
       tpx = [t_x^3 t_x^2 t_x 1]
       tpy = [t_y^3 t_y^2 t_y 1]
       tpz = [t_z^3 t_z^2 t_z 1] */
    __m128 ipart  = _mm_cvtepi32_ps (intpart);
    __m128 txtytz = _mm_sub_ps (uxuyuz, ipart);
    __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
    __m128 t2     = _mm_mul_ps (txtytz, txtytz);
    __m128 t3     = _mm_mul_ps (t2, txtytz);
    __m128 tpx    = t3;
    __m128 tpy    = t2;
    __m128 tpz    = txtytz;
    __m128 zero   = one;
    _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

    /* a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
       da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
       A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3] */
    __m128 a, b, c, da, db, dc,
	cP[4], dcP[4], bcP, dbcP, bdcP,
	tmp0, tmp1, tmp2, tmp3;
  
    /* x-dependent vectors */
    _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpx,   a);
    _MM_MATVEC4_PS (A_s[4], A_s[5], A_s[6], A_s[7], tpx,  da);
    /* y-dependent vectors */
    _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpy,   b);
    _MM_MATVEC4_PS (A_s[4], A_s[5], A_s[6], A_s[7], tpy,  db);
    /* z-dependent vectors */
    _MM_MATVEC4_PS (A_s[0], A_s[1], A_s[2], A_s[3], tpz,   c);
    _MM_MATVEC4_PS (A_s[4], A_s[5], A_s[6], A_s[7], tpz,  dc);

    /* Compute cP, dcP, and d2cP products 1/4 at a time to maximize
       register reuse and avoid rerereading from memory or cache.
       1st quarter */
    tmp0 = _mm_loadu_ps (P(0,0));  tmp1 = _mm_loadu_ps (P(0,1));
    tmp2 = _mm_loadu_ps (P(0,2));  tmp3 = _mm_loadu_ps (P(0,3));
    _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[0]);
    _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[0]);
    /* 2nd quarter */
    tmp0 = _mm_loadu_ps (P(1,0));  tmp1 = _mm_loadu_ps (P(1,1));
    tmp2 = _mm_loadu_ps (P(1,2));  tmp3 = _mm_loadu_ps (P(1,3));
    _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[1]);
    _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[1]);
    /* 3rd quarter */
    tmp0 = _mm_loadu_ps (P(2,0));  tmp1 = _mm_loadu_ps (P(2,1));
    tmp2 = _mm_loadu_ps (P(2,2));  tmp3 = _mm_loadu_ps (P(2,3));
    _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[2]);
    _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[2]);
    /* 4th quarter */
    tmp0 = _mm_loadu_ps (P(3,0));  tmp1 = _mm_loadu_ps (P(3,1));
    tmp2 = _mm_loadu_ps (P(3,2));  tmp3 = _mm_loadu_ps (P(3,3));
    _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,   c,   cP[3]);
    _MM_MATVEC4_PS (tmp0, tmp1, tmp2, tmp3,  dc,  dcP[3]);
  
    /* Now compute bcP, dbcP, bdcP, d2bcP, bd2cP, and dbdc products */
    _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3],   b,   bcP);
    _MM_MATVEC4_PS (  cP[0],   cP[1],   cP[2],   cP[3],  db,  dbcP);
    _MM_MATVEC4_PS ( dcP[0],  dcP[1],  dcP[2],  dcP[3],   b,  bdcP);

    /* Compute value */
    _MM_DOT4_PS (a, bcP, *val);
    /* Compute gradient */
    _MM_DOT4_PS (da, bcP, grad[0]);
    _MM_DOT4_PS (a, dbcP, grad[1]);
    _MM_DOT4_PS (a, bdcP, grad[2]);
  
    /* Multiply gradients and hessians by appropriate grid inverses */
    float dxInv = spline->x_grid.delta_inv;
    float dyInv = spline->y_grid.delta_inv;
    float dzInv = spline->z_grid.delta_inv;
    grad[0] *= dxInv;
    grad[1] *= dyInv;
    grad[2] *= dzInv;
#undef P
}
#endif

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
    intptr_t Mx, Nx, N;
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
#ifdef HAVE_SSE
    if (N % 4) 
	N += 4 - (N % 4);
#endif 
    spline->x_stride = N;
    x_grid.delta_inv = 1.0/x_grid.delta;
    spline->x_grid   = x_grid;
#ifndef HAVE_SSE
    spline->coefs = malloc ((intptr_t)sizeof(float)*(intptr_t)Nx*(intptr_t)N);
#else
    posix_memalign ((void**)&spline->coefs, 64, (sizeof(float)*Nx*N));
    init_sse_data();    
#endif
    spline->nc = (intptr_t)sizeof(float)*(intptr_t)Nx*(intptr_t)N;

    return spline;
}

void set_multi_UBspline_1d_s (multi_UBspline_1d_s *spline, int num,
                              float *data)
{
    float *coefs = spline->coefs + num;
    intptr_t xs = spline->x_stride;

    find_coefs_1d_s (spline->x_grid, spline->xBC, data, 1, 
		     coefs, xs);
}

multi_UBspline_2d_s*
create_multi_UBspline_2d_s (Ugrid x_grid, Ugrid y_grid,
                            BCtype_s xBC, BCtype_s yBC, int num_splines)
{
    intptr_t Mx, My, Nx, Ny, N;
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
#ifdef HAVE_SSE
    if (N % 4) 
	N += 4 - (N % 4);
#endif
    spline->x_stride = Ny*N;
    spline->y_stride = N;
#ifndef HAVE_SSE
    spline->coefs = malloc ((intptr_t)sizeof(float)*(intptr_t)Nx*(intptr_t)Ny*(intptr_t)N);
#else
    posix_memalign ((void**)&spline->coefs, 64, 
		    sizeof(float)*Nx*Ny*N);
    init_sse_data();
#endif
    spline->nc = (intptr_t)sizeof(float)*(intptr_t)Nx*(intptr_t)Ny*(intptr_t)N;

    return spline;
}

void set_multi_UBspline_2d_s (multi_UBspline_2d_s* spline, int num, float *data)
{
    intptr_t Mx = spline->x_grid.num;
    intptr_t My = spline->y_grid.num;
    intptr_t Nx, Ny, ys, iy, ix;
    float *coefs = spline->coefs + num;
    intptr_t doffset;
    intptr_t coffset;

    if (spline->xBC.lCode == PERIODIC || spline->xBC.lCode == ANTIPERIODIC)     
	Nx = Mx+3;
    else                                   
	Nx = Mx+2;
    if (spline->yBC.lCode == PERIODIC || spline->yBC.lCode == ANTIPERIODIC)     
	Ny = My+3;
    else                                   
	Ny = My+2;

    ys = spline->y_stride;
    /* First, solve in the X-direction */
    for (iy=0; iy<My; iy++) {
	doffset = iy;
	coffset = iy*ys;
	find_coefs_1d_s (spline->x_grid, spline->xBC, data+doffset, (intptr_t)My,
			 coefs+coffset, (intptr_t)Ny*ys);
    }
  
    /* Now, solve in the Y-direction */
    for (ix=0; ix<Nx; ix++) {
	doffset = ix*Ny*ys;
	coffset = ix*Ny*ys;
	find_coefs_1d_s (spline->y_grid, spline->yBC, coefs+doffset, (intptr_t)ys, 
			 coefs+coffset, (intptr_t)ys);
    }
}

multi_UBspline_3d_s*
create_multi_UBspline_3d_s (Ugrid x_grid, Ugrid y_grid, Ugrid z_grid,
                            BCtype_s xBC, BCtype_s yBC, BCtype_s zBC,
                            int num_splines)
{
    intptr_t Mx, My, Mz, Nx, Ny, Nz, N;
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
#ifdef HAVE_SSE
    if (N % 4) 
	N += 4 - (N % 4);
#endif
    spline->x_stride      = Ny*Nz*N;
    spline->y_stride      = Nz*N;
    spline->z_stride      = N;
#ifndef HAVE_SSE
    spline->coefs = malloc ((intptr_t)sizeof(float)*(intptr_t)Nx*(intptr_t)Ny*(intptr_t)Nz*(intptr_t)N);
#else
    posix_memalign ((void**)&spline->coefs, 64, 
		    ((intptr_t)sizeof(float)*Nx*Ny*Nz*N));
    init_sse_data();
#endif
    spline->nc = (intptr_t)sizeof(float)*(intptr_t)Nx*(intptr_t)Ny*(intptr_t)Nz*(intptr_t)N;

    return spline;
}

void set_multi_UBspline_3d_s (multi_UBspline_3d_s* spline, int num, float *data)
{
    intptr_t Mx = spline->x_grid.num;
    intptr_t My = spline->y_grid.num;
    intptr_t Mz = spline->z_grid.num;
    intptr_t Nx, Ny, Nz, zs, ix, iy, iz;
    float *coefs = spline->coefs + num;
    intptr_t doffset;
    intptr_t coffset;

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

    zs = spline->z_stride;
    /* First, solve in the X-direction */
    for (iy=0; iy<My; iy++) 
	for (iz=0; iz<Mz; iz++) {
	    doffset = iy*Mz+iz;
	    coffset = (iy*Nz+iz)*zs;
	    find_coefs_1d_s (spline->x_grid, spline->xBC, data+doffset, (intptr_t)My*Mz,
			     coefs+coffset, (intptr_t)Ny*Nz*zs);
	}
  
    /* Now, solve in the Y-direction */
    for (ix=0; ix<Nx; ix++) 
	for (iz=0; iz<Nz; iz++) {
	    doffset = (ix*Ny*Nz + iz)*zs;
	    coffset = (ix*Ny*Nz + iz)*zs;
	    find_coefs_1d_s (spline->y_grid, spline->yBC, coefs+doffset, (intptr_t)Nz*zs, 
			     coefs+coffset, (intptr_t)Nz*zs);
	}

    /* Now, solve in the Z-direction */
    for (ix=0; ix<Nx; ix++) 
	for (iy=0; iy<Ny; iy++) {
	    doffset = ((ix*Ny+iy)*Nz)*zs;
	    coffset = ((ix*Ny+iy)*Nz)*zs;
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
    intptr_t n, i;
    float ux;
    float ipartx, tx;
 
    float tpx[4], a[4];
    float* restrict coefs = spline->coefs;
    intptr_t ix;
    intptr_t xs;

    x -= spline->x_grid.start;
    ux = x*spline->x_grid.delta_inv;
    tx = modff (ux, &ipartx);  ix = (intptr_t) ipartx;
 
    tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;

    a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
    a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
    a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
    a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);

    xs = spline->x_stride;

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
    intptr_t n, i;
    float ux;
    float ipartx, tx;
    float tpx[4], a[4], da[4];
    float* restrict coefs = spline->coefs;
    intptr_t ix, xs;
    float dxInv;

    x -= spline->x_grid.start;
    ux = x*spline->x_grid.delta_inv;

    tx = modff (ux, &ipartx);  ix = (intptr_t) ipartx;
  
    tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;

    a[0]  = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
    a[1]  = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
    a[2]  = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
    a[3]  = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);
    da[0] = (dAf[ 0]*tpx[0] + dAf[ 1]*tpx[1] + dAf[ 2]*tpx[2] + dAf[ 3]*tpx[3]);
    da[1] = (dAf[ 4]*tpx[0] + dAf[ 5]*tpx[1] + dAf[ 6]*tpx[2] + dAf[ 7]*tpx[3]);
    da[2] = (dAf[ 8]*tpx[0] + dAf[ 9]*tpx[1] + dAf[10]*tpx[2] + dAf[11]*tpx[3]);
    da[3] = (dAf[12]*tpx[0] + dAf[13]*tpx[1] + dAf[14]*tpx[2] + dAf[15]*tpx[3]);

    xs = spline->x_stride;

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
  
    dxInv = spline->x_grid.delta_inv;
    for (n=0; n<spline->num_splines; n++) 
	grads[n] *= dxInv;
}

/************************************************************/
/* 2D single precision, real evaulation functions           */
/************************************************************/
#ifndef HAVE_SSE
void eval_multi_UBspline_2d_s (multi_UBspline_2d_s *spline,
                               double x, double y,
                               float* restrict vals)
{
    intptr_t n, i, j;
    float ux;
    float uy;
    float ipartx, iparty, tx, ty;
    float tpx[4], tpy[4], a[4], b[4];
    float* restrict coefs = spline->coefs;
    intptr_t ix, iy;
    intptr_t xs, ys;
    float prefactor;

    x -= spline->x_grid.start;
    y -= spline->y_grid.start;
    ux = x*spline->x_grid.delta_inv;
    uy = y*spline->y_grid.delta_inv;

    tx = modff (ux, &ipartx);  ix = (intptr_t) ipartx;
    ty = modff (uy, &iparty);  iy = (intptr_t) iparty;
 
    tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
    tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
 

    a[0] = (Af[ 0]*tpx[0] + Af[ 1]*tpx[1] + Af[ 2]*tpx[2] + Af[ 3]*tpx[3]);
    a[1] = (Af[ 4]*tpx[0] + Af[ 5]*tpx[1] + Af[ 6]*tpx[2] + Af[ 7]*tpx[3]);
    a[2] = (Af[ 8]*tpx[0] + Af[ 9]*tpx[1] + Af[10]*tpx[2] + Af[11]*tpx[3]);
    a[3] = (Af[12]*tpx[0] + Af[13]*tpx[1] + Af[14]*tpx[2] + Af[15]*tpx[3]);

    b[0] = (Af[ 0]*tpy[0] + Af[ 1]*tpy[1] + Af[ 2]*tpy[2] + Af[ 3]*tpy[3]);
    b[1] = (Af[ 4]*tpy[0] + Af[ 5]*tpy[1] + Af[ 6]*tpy[2] + Af[ 7]*tpy[3]);
    b[2] = (Af[ 8]*tpy[0] + Af[ 9]*tpy[1] + Af[10]*tpy[2] + Af[11]*tpy[3]);
    b[3] = (Af[12]*tpy[0] + Af[13]*tpy[1] + Af[14]*tpy[2] + Af[15]*tpy[3]);

    xs = spline->x_stride;
    ys = spline->y_stride;

    for (n=0; n<spline->num_splines; n++)
	vals[n] = 0.0;

    for (i=0; i<4; i++)
	for (j=0; j<4; j++) {
	    prefactor = a[i]*b[j];
	    coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
	    for (n=0; n<spline->num_splines; n++) 
		vals[n] += prefactor*coefs[n];
	}
}
#else
void eval_multi_UBspline_2d_s (multi_UBspline_2d_s *spline,
                               double x, double y,
                               float* restrict vals)
{
    _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);  
    _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
    /* SSE mesh point determination */
    __m128 xy        = _mm_set_ps (x, y, 0.0, 0.0);
    __m128 x0y0      = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 0.0, 0.0);
    __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 0.0, 0.0);
    xy = _mm_sub_ps (xy, x0y0);
    /* ux = (x - x0)/delta_x and same for y */
    __m128 uxuy    = _mm_mul_ps (xy, delta_inv);
    /* intpart = trunc (ux, uy) */
    __m128i intpart  = _mm_cvttps_epi32(uxuy);
    __m128i ixiy;
    _mm_storeu_si128 (&ixiy, intpart);
    /* Store to memory for use in C expressions
       xmm registers are stored to memory in reverse order */
    int ix = ((int *)&ixiy)[3];
    int iy = ((int *)&ixiy)[2];

    intptr_t xs = spline->x_stride;
    intptr_t ys = spline->y_stride;

    /* Now compute the vectors:
       tpx = [t_x^3 t_x^2 t_x 1]
       tpy = [t_y^3 t_y^2 t_y 1]
       tpz = [t_z^3 t_z^2 t_z 1] */
    __m128 ipart  = _mm_cvtepi32_ps (intpart);
    __m128 txty   = _mm_sub_ps (uxuy, ipart);
    __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
    __m128 t2     = _mm_mul_ps (txty, txty);
    __m128 t3     = _mm_mul_ps (t2, txty);
    __m128 tpx    = t3;
    __m128 tpy    = t2;
    __m128 tpz    = txty;
    __m128 zero   = one;
    _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

    /* a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
       da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
       A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3] */
    __m128 a4, b4, da4, db4;
    /* x-dependent vectors */
    _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpx,   a4);
    _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpx,  da4);
    /* y-dependent vectors */
    _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpy,   b4);
    _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpy,  db4);

    __m128 a[4], b[4];
    __m128 tmp;

    /* Unpack a values */
    tmp=_mm_unpacklo_ps(  a4,   a4);   a[0]=_mm_unpacklo_ps(tmp, tmp);   a[1]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpackhi_ps(  a4,   a4);   a[2]=_mm_unpacklo_ps(tmp, tmp);   a[3]=_mm_unpackhi_ps(tmp, tmp);

    /* Unpack b values */
    tmp=_mm_unpacklo_ps(  b4,   b4);   b[0]=_mm_unpacklo_ps(tmp, tmp);   b[1]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpackhi_ps(  b4,   b4);   b[2]=_mm_unpacklo_ps(tmp, tmp);   b[3]=_mm_unpackhi_ps(tmp, tmp);

    int N = spline->num_splines;
    int Nm = (N+3)/4;

    __m128 mvals[Nm];

    /* Zero out values; */
    for (int n=0; n<Nm; n++)     mvals[n] = _mm_setzero_ps();

    /* Main compute loop */
    __m128 ab;
    for (int i=0; i<4; i++)
	for (int j=0; j<4; j++) {
	    ab      = _mm_mul_ps (  a[i],   b[j]);
      
	    __m128* restrict coefs = (__m128*)(spline->coefs + (ix+i)*xs + (iy+j)*ys);
	    for (int n=0; n<Nm; n++) 
		mvals[n]     = _mm_add_ps (mvals[n],     _mm_mul_ps(  ab   , coefs[n]));
	}
  
    /* Now, store results back */
    for (int n=0; n<N; n++) 
	vals[n] = ((float*)mvals)[n];
}
#endif

#ifndef HAVE_SSE
void eval_multi_UBspline_2d_s_vg (multi_UBspline_2d_s *spline,
                                  double x, double y,
                                  float* restrict vals,
                                  float* restrict grads)
{
    intptr_t n, i, j;
    float ux;
    float uy;
    float ipartx, iparty, tx, ty;
    intptr_t ix, iy;
    float tpx[4], tpy[4], a[4], b[4], da[4], db[4];
    float* restrict coefs = spline->coefs;
    intptr_t xs, ys;
    float ab, dab[2];
    float dxInv, dyInv;

    x -= spline->x_grid.start;
    y -= spline->y_grid.start;
    ux = x*spline->x_grid.delta_inv;
    uy = y*spline->y_grid.delta_inv;

    tx = modff (ux, &ipartx);  ix = (intptr_t) ipartx;
    ty = modff (uy, &iparty);  iy = (intptr_t) iparty;
  
    tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
    tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;

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

    xs = spline->x_stride;
    ys = spline->y_stride;

    for (n=0; n<spline->num_splines; n++) {
	vals[n] = 0.0;
	grads[2*n+0] = grads[2*n+1] = grads[2*n+2] = 0.0;
    }

    for (i=0; i<4; i++)
	for (j=0; j<4; j++) {
	    ab = a[i]*b[j];
	    dab[0] = da[i]* b[j];
	    dab[1] =  a[i]*db[j];
      
	    coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys);
	    for (n=0; n<spline->num_splines; n++) {
		vals [n]     +=   ab   *coefs[n];
		grads[2*n+0] +=  dab[0]*coefs[n];
		grads[2*n+1] +=  dab[1]*coefs[n];
	    }
	}

    dxInv = spline->x_grid.delta_inv;
    dyInv = spline->y_grid.delta_inv;
    for (n=0; n<spline->num_splines; n++) {
	grads[2*n+0] *= dxInv;
	grads[2*n+1] *= dyInv;
    }
}
#else
void eval_multi_UBspline_2d_s_vg (multi_UBspline_2d_s *spline,
                                  double x, double y,
                                  float* restrict vals,
                                  float* restrict grads)
{
    _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);  
    _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
    _mm_prefetch ((const char*)  &A_s[ 4],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 5],_MM_HINT_T0);  
    _mm_prefetch ((const char*)  &A_s[ 6],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 7],_MM_HINT_T0);
    /* SSE mesh point determination */
    __m128 xy        = _mm_set_ps (x, y, 0.0, 0.0);
    __m128 x0y0      = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 0.0, 0.0);
    __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 0.0, 0.0);
    xy = _mm_sub_ps (xy, x0y0);
    /* ux = (x - x0)/delta_x and same for y */
    __m128 uxuy    = _mm_mul_ps (xy, delta_inv);
    /* intpart = trunc (ux, uy) */
    __m128i intpart  = _mm_cvttps_epi32(uxuy);
    __m128i ixiy;
    _mm_storeu_si128 (&ixiy, intpart);
    /* Store to memory for use in C expressions
       xmm registers are stored to memory in reverse order */
    int ix = ((int *)&ixiy)[3];
    int iy = ((int *)&ixiy)[2];

    intptr_t xs = spline->x_stride;
    intptr_t ys = spline->y_stride;

    /* Now compute the vectors:
       tpx = [t_x^3 t_x^2 t_x 1]
       tpy = [t_y^3 t_y^2 t_y 1]
       tpz = [t_z^3 t_z^2 t_z 1] */
    __m128 ipart  = _mm_cvtepi32_ps (intpart);
    __m128 txty   = _mm_sub_ps (uxuy, ipart);
    __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
    __m128 t2     = _mm_mul_ps (txty, txty);
    __m128 t3     = _mm_mul_ps (t2, txty);
    __m128 tpx    = t3;
    __m128 tpy    = t2;
    __m128 tpz    = txty;
    __m128 zero   = one;
    _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

    /* a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
       da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
       A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3] */
    __m128 a4, b4, da4, db4;
    /* x-dependent vectors */
    _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpx,   a4);
    _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpx,  da4);
    /* y-dependent vectors */
    _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpy,   b4);
    _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpy,  db4);

    __m128 a[4], b[4], da[4], db[4];
    __m128 tmp;

    /* Unpack a values */
    tmp=_mm_unpacklo_ps(  a4,   a4);   a[0]=_mm_unpacklo_ps(tmp, tmp);   a[1]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpackhi_ps(  a4,   a4);   a[2]=_mm_unpacklo_ps(tmp, tmp);   a[3]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpacklo_ps( da4,  da4);  da[0]=_mm_unpacklo_ps(tmp, tmp);  da[1]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpackhi_ps( da4,  da4);  da[2]=_mm_unpacklo_ps(tmp, tmp);  da[3]=_mm_unpackhi_ps(tmp, tmp);

    /* Unpack b values */
    tmp=_mm_unpacklo_ps(  b4,   b4);   b[0]=_mm_unpacklo_ps(tmp, tmp);   b[1]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpackhi_ps(  b4,   b4);   b[2]=_mm_unpacklo_ps(tmp, tmp);   b[3]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpacklo_ps( db4,  db4);  db[0]=_mm_unpacklo_ps(tmp, tmp);  db[1]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpackhi_ps( db4,  db4);  db[2]=_mm_unpacklo_ps(tmp, tmp);  db[3]=_mm_unpackhi_ps(tmp, tmp);

    int N = spline->num_splines;
    int Nm = (N+3)/4;

    __m128 mvals[Nm], mgrad[2*Nm];

    /* Zero out values; */
    for (int n=0; n<Nm; n++)     mvals[n] = _mm_setzero_ps();
    for (int n=0; n<2*Nm; n++)   mgrad[n] = _mm_setzero_ps();

    /* Main compute loop */
    __m128 ab, dab[2];
    for (int i=0; i<4; i++)
	for (int j=0; j<4; j++) {
	    ab      = _mm_mul_ps (  a[i],   b[j]);
	    dab[0]  = _mm_mul_ps ( da[i],   b[j]);
	    dab[1]  = _mm_mul_ps (  a[i],  db[j]);
      
	    __m128* restrict coefs = (__m128*)(spline->coefs + (ix+i)*xs + (iy+j)*ys);
	    for (int n=0; n<Nm; n++) {
		mvals[n]     = _mm_add_ps (mvals[n],     _mm_mul_ps(  ab   , coefs[n]));
		mgrad[2*n+0] = _mm_add_ps (mgrad[2*n+0], _mm_mul_ps( dab[0], coefs[n]));
		mgrad[2*n+1] = _mm_add_ps (mgrad[2*n+1], _mm_mul_ps( dab[1], coefs[n]));
	    }
	}
  
    /* Now, store results back */
    float dxInv = spline->x_grid.delta_inv;
    float dyInv = spline->y_grid.delta_inv;
    for (int n=0; n<N; n++) {
	int nd4 = n>>2;
	int nm4 = n & 3;
	vals[n] = ((float*)mvals)[n];
	grads[2*n+0] = ((float*)mgrad)[nd4*8 + 4*0 + nm4] * dxInv;
	grads[2*n+1] = ((float*)mgrad)[nd4*8 + 4*1 + nm4] * dyInv;
    }
}
#endif

/************************************************************/
/* 3D single precision, real evaulation functions           */
/************************************************************/
#ifndef HAVE_SSE
void eval_multi_UBspline_3d_s (multi_UBspline_3d_s *spline,
                               double x, double y, double z,
                               float* restrict vals)
{
    intptr_t n, i, j, k;
    float ux;
    float uy;
    float uz;
    float ipartx, iparty, ipartz, tx, ty, tz;
    float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4];
    float* restrict coefs = spline->coefs;
    intptr_t ix, iy, iz;
    intptr_t xs, ys, zs;
    float abc;

    x -= spline->x_grid.start;
    y -= spline->y_grid.start;
    z -= spline->z_grid.start;
    ux = x*spline->x_grid.delta_inv;
    uy = y*spline->y_grid.delta_inv;
    uz = z*spline->z_grid.delta_inv;

    tx = modff (ux, &ipartx);  ix = (intptr_t) ipartx;
    ty = modff (uy, &iparty);  iy = (intptr_t) iparty;
    tz = modff (uz, &ipartz);  iz = (intptr_t) ipartz;
  
    tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
    tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
    tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;

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
    zs = spline->z_stride;

    for (n=0; n<spline->num_splines; n++)
	vals[n] = 0.0;

    for (i=0; i<4; i++)
	for (j=0; j<4; j++) 
	    for (k=0; k<4; k++) {
		abc = a[i]*b[j]*c[k];
		coefs = spline->coefs + ((ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
		for (n=0; n<spline->num_splines; n++) 
		    vals[n] += abc*coefs[n];
	    }
}
#else
void eval_multi_UBspline_3d_s (multi_UBspline_3d_s *spline,
                               double x, double y, double z,
                               float* restrict vals)
{
    _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);  
    _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);

    /* SSE mesh point determination */
    __m128 xyz       = _mm_set_ps (x, y, z, 0.0);
    __m128 x0y0z0    = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 
				   spline->z_grid.start, 0.0);
    __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 
				   spline->z_grid.delta_inv, 0.0);
    xyz = _mm_sub_ps (xyz, x0y0z0);
    /* ux = (x - x0)/delta_x and same for y and z */
    __m128 uxuyuz    = _mm_mul_ps (xyz, delta_inv);
    /* intpart = trunc (ux, uy, uz) */
    __m128i intpart  = _mm_cvttps_epi32(uxuyuz);
    __m128i ixiyiz;
    _mm_storeu_si128 (&ixiyiz, intpart);
    /* Store to memory for use in C expressions
       xmm registers are stored to memory in reverse order */
    int ix = ((int *)&ixiyiz)[3];
    int iy = ((int *)&ixiyiz)[2];
    int iz = ((int *)&ixiyiz)[1];

    intptr_t xs = spline->x_stride;
    intptr_t ys = spline->y_stride;
    intptr_t zs = spline->z_stride;

    /* Now compute the vectors:
       tpx = [t_x^3 t_x^2 t_x 1]
       tpy = [t_y^3 t_y^2 t_y 1]
       tpz = [t_z^3 t_z^2 t_z 1] */
    __m128 ipart  = _mm_cvtepi32_ps (intpart);
    __m128 txtytz = _mm_sub_ps (uxuyuz, ipart);
    __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
    __m128 t2     = _mm_mul_ps (txtytz, txtytz);
    __m128 t3     = _mm_mul_ps (t2, txtytz);
    __m128 tpx    = t3;
    __m128 tpy    = t2;
    __m128 tpz    = txtytz;
    __m128 zero   = one;
    _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

    /* a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
       A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3] */
    __m128 a4, b4, c4;
  
    /* x-dependent vectors */
    _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpx,   a4);
    /* y-dependent vectors */
    _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpy,   b4);
    /* z-dependent vectors */
    _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpz,   c4);

    __m128 a[4], b[4], c[4];
    __m128 tmp;

    /* Unpack a values */
    tmp=_mm_unpacklo_ps(  a4,   a4);   a[0]=_mm_unpacklo_ps(tmp, tmp);   a[1]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpackhi_ps(  a4,   a4);   a[2]=_mm_unpacklo_ps(tmp, tmp);   a[3]=_mm_unpackhi_ps(tmp, tmp);

    /* Unpack b values */
    tmp=_mm_unpacklo_ps(  b4,   b4);   b[0]=_mm_unpacklo_ps(tmp, tmp);   b[1]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpackhi_ps(  b4,   b4);   b[2]=_mm_unpacklo_ps(tmp, tmp);   b[3]=_mm_unpackhi_ps(tmp, tmp);

    /* Unpack c values */
    tmp=_mm_unpacklo_ps(  c4,   c4);   c[0]=_mm_unpacklo_ps(tmp, tmp);   c[1]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpackhi_ps(  c4,   c4);   c[2]=_mm_unpacklo_ps(tmp, tmp);   c[3]=_mm_unpackhi_ps(tmp, tmp);

    int N = spline->num_splines;
    int Nm = (N+3)/4;

    __m128 mvals[Nm];

    /* Zero out values; */
    for (int n=0; n<Nm; n++)     mvals[n] = _mm_setzero_ps();

    /* Main compute loop */
    __m128 abc;
    for (int i=0; i<4; i++)
	for (int j=0; j<4; j++)
	    for (int k=0; k<4; k++) {
		abc      = _mm_mul_ps (  a[i], _mm_mul_ps(  b[j],  c[k]));
		__m128* restrict coefs = (__m128*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
		for (int n=0; n<Nm; n++) 
		    mvals[n]     = _mm_add_ps (mvals[n], _mm_mul_ps(  abc   , coefs[n]));
	    }
  
    /* Now, store results back */
    for (int n=0; n<N; n++) 
	vals[n] = ((float*)mvals)[n];
  
}
#endif

#ifndef HAVE_SSE
void eval_multi_UBspline_3d_s_vg (multi_UBspline_3d_s *spline,
                                  double x, double y, double z,
                                  float* restrict vals,
                                  float* restrict grads)
{
    intptr_t n, i, j, k;
    float ux;
    float uy;
    float uz;
    float ipartx, iparty, ipartz, tx, ty, tz;
    intptr_t ix, iy, iz;
    float tpx[4], tpy[4], tpz[4], a[4], b[4], c[4], 
	da[4], db[4], dc[4];
    float* restrict coefs = spline->coefs;
    intptr_t xs, ys, zs;
    float abc, dabc[3];
    float dxInv, dyInv, dzInv;

    x -= spline->x_grid.start;
    y -= spline->y_grid.start;
    z -= spline->z_grid.start;
    ux = x*spline->x_grid.delta_inv;
    uy = y*spline->y_grid.delta_inv;
    uz = z*spline->z_grid.delta_inv;

    tx = modff (ux, &ipartx);  ix = (intptr_t) ipartx;
    ty = modff (uy, &iparty);  iy = (intptr_t) iparty;
    tz = modff (uz, &ipartz);  iz = (intptr_t) ipartz;
  
    tpx[0] = tx*tx*tx;  tpx[1] = tx*tx;  tpx[2] = tx;  tpx[3] = 1.0;
    tpy[0] = ty*ty*ty;  tpy[1] = ty*ty;  tpy[2] = ty;  tpy[3] = 1.0;
    tpz[0] = tz*tz*tz;  tpz[1] = tz*tz;  tpz[2] = tz;  tpz[3] = 1.0;

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

    xs = spline->x_stride;
    ys = spline->y_stride;
    zs = spline->z_stride;

    for (n=0; n<spline->num_splines; n++) {
	vals[n] = 0.0;
	grads[3*n+0] = grads[3*n+1] = grads[3*n+2] = 0.0;
    }

    for (i=0; i<4; i++)
	for (j=0; j<4; j++) 
	    for (k=0; k<4; k++) {
		abc = a[i]*b[j]*c[k];
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

    dxInv = spline->x_grid.delta_inv;
    dyInv = spline->y_grid.delta_inv;
    dzInv = spline->z_grid.delta_inv; 
    for (n=0; n<spline->num_splines; n++) {
	grads[3*n+0] *= dxInv;
	grads[3*n+1] *= dyInv;
	grads[3*n+2] *= dzInv;
    }
}
#else
void eval_multi_UBspline_3d_s_vg (multi_UBspline_3d_s *spline,
                                  double x, double y, double z,
                                  float* restrict vals,
                                  float* restrict grads)
{
    _mm_prefetch ((const char*)  &A_s[ 0],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 1],_MM_HINT_T0);  
    _mm_prefetch ((const char*)  &A_s[ 2],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 3],_MM_HINT_T0);
    _mm_prefetch ((const char*)  &A_s[ 4],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 5],_MM_HINT_T0);  
    _mm_prefetch ((const char*)  &A_s[ 6],_MM_HINT_T0);  _mm_prefetch ((const char*)  &A_s[ 7],_MM_HINT_T0);

    /* SSE mesh point determination */
    __m128 xyz       = _mm_set_ps (x, y, z, 0.0);
    __m128 x0y0z0    = _mm_set_ps (spline->x_grid.start,  spline->y_grid.start, 
				   spline->z_grid.start, 0.0);
    __m128 delta_inv = _mm_set_ps (spline->x_grid.delta_inv,spline->y_grid.delta_inv, 
				   spline->z_grid.delta_inv, 0.0);
    xyz = _mm_sub_ps (xyz, x0y0z0);
    /* ux = (x - x0)/delta_x and same for y and z */
    __m128 uxuyuz    = _mm_mul_ps (xyz, delta_inv);
    /* intpart = trunc (ux, uy, uz) */
    __m128i intpart  = _mm_cvttps_epi32(uxuyuz);
    __m128i ixiyiz;
    _mm_storeu_si128 (&ixiyiz, intpart);
    /* Store to memory for use in C expressions
       xmm registers are stored to memory in reverse order */
    int ix = ((int *)&ixiyiz)[3];
    int iy = ((int *)&ixiyiz)[2];
    int iz = ((int *)&ixiyiz)[1];

    intptr_t xs = spline->x_stride;
    intptr_t ys = spline->y_stride;
    intptr_t zs = spline->z_stride;

    /* Now compute the vectors:
       tpx = [t_x^3 t_x^2 t_x 1]
       tpy = [t_y^3 t_y^2 t_y 1]
       tpz = [t_z^3 t_z^2 t_z 1] */
    __m128 ipart  = _mm_cvtepi32_ps (intpart);
    __m128 txtytz = _mm_sub_ps (uxuyuz, ipart);
    __m128 one    = _mm_set_ps (1.0, 1.0, 1.0, 1.0);
    __m128 t2     = _mm_mul_ps (txtytz, txtytz);
    __m128 t3     = _mm_mul_ps (t2, txtytz);
    __m128 tpx    = t3;
    __m128 tpy    = t2;
    __m128 tpz    = txtytz;
    __m128 zero   = one;
    _MM_TRANSPOSE4_PS(zero, tpz, tpy, tpx);

    /* a  =  A * tpx,   b =  A * tpy,   c =  A * tpz
       da = dA * tpx,  db = dA * tpy,  dc = dA * tpz, etc.
       A is 4x4 matrix given by the rows A_s[0], A_s[1], A_s[ 2], A_s[ 3] */
    __m128 a4, b4, c4, da4, db4, dc4;
  
    /* x-dependent vectors */
    _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpx,   a4);
    _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpx,  da4);
    /* y-dependent vectors */
    _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpy,   b4);
    _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpy,  db4);
    /* z-dependent vectors */
    _MM_MATVEC4_PS (A_s[ 0], A_s[ 1], A_s[ 2], A_s[ 3], tpz,   c4);
    _MM_MATVEC4_PS (A_s[ 4], A_s[ 5], A_s[ 6], A_s[ 7], tpz,  dc4);

    __m128 a[4], b[4], c[4], da[4], db[4], dc[4];
    __m128 tmp;

    /* Unpack a values */
    tmp=_mm_unpacklo_ps(  a4,   a4);   a[0]=_mm_unpacklo_ps(tmp, tmp);   a[1]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpackhi_ps(  a4,   a4);   a[2]=_mm_unpacklo_ps(tmp, tmp);   a[3]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpacklo_ps( da4,  da4);  da[0]=_mm_unpacklo_ps(tmp, tmp);  da[1]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpackhi_ps( da4,  da4);  da[2]=_mm_unpacklo_ps(tmp, tmp);  da[3]=_mm_unpackhi_ps(tmp, tmp);

    /* Unpack b values */
    tmp=_mm_unpacklo_ps(  b4,   b4);   b[0]=_mm_unpacklo_ps(tmp, tmp);   b[1]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpackhi_ps(  b4,   b4);   b[2]=_mm_unpacklo_ps(tmp, tmp);   b[3]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpacklo_ps( db4,  db4);  db[0]=_mm_unpacklo_ps(tmp, tmp);  db[1]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpackhi_ps( db4,  db4);  db[2]=_mm_unpacklo_ps(tmp, tmp);  db[3]=_mm_unpackhi_ps(tmp, tmp);

    /* Unpack c values */
    tmp=_mm_unpacklo_ps(  c4,   c4);   c[0]=_mm_unpacklo_ps(tmp, tmp);   c[1]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpackhi_ps(  c4,   c4);   c[2]=_mm_unpacklo_ps(tmp, tmp);   c[3]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpacklo_ps( dc4,  dc4);  dc[0]=_mm_unpacklo_ps(tmp, tmp);  dc[1]=_mm_unpackhi_ps(tmp, tmp);
    tmp=_mm_unpackhi_ps( dc4,  dc4);  dc[2]=_mm_unpacklo_ps(tmp, tmp);  dc[3]=_mm_unpackhi_ps(tmp, tmp);

    int N = spline->num_splines;
    int Nm = (N+3)/4;

    __m128 mvals[Nm], mgrad[3*Nm];

    /* Zero out values; */
    for (int n=0; n<Nm; n++)     mvals[n] = _mm_setzero_ps();
    for (int n=0; n<3*Nm; n++)   mgrad[n] = _mm_setzero_ps();

    /* Main compute loop */
    __m128 abc, dabc[3];
    for (int i=0; i<4; i++)
	for (int j=0; j<4; j++)
	    for (int k=0; k<4; k++) {
		abc      = _mm_mul_ps (  a[i], _mm_mul_ps(  b[j],  c[k]));
		dabc[0]  = _mm_mul_ps ( da[i], _mm_mul_ps(  b[j],  c[k]));
		dabc[1]  = _mm_mul_ps (  a[i], _mm_mul_ps( db[j],  c[k]));
		dabc[2]  = _mm_mul_ps (  a[i], _mm_mul_ps(  b[j],  dc[k]));

		__m128* restrict coefs = (__m128*)(spline->coefs + (ix+i)*xs + (iy+j)*ys + (iz+k)*zs);
		for (int n=0; n<Nm; n++) {
		    mvals[n]     = _mm_add_ps (mvals[n], _mm_mul_ps(  abc   , coefs[n]));
		    mgrad[3*n+0] = _mm_add_ps (mgrad[3*n+0], _mm_mul_ps( dabc[0], coefs[n]));
		    mgrad[3*n+1] = _mm_add_ps (mgrad[3*n+1], _mm_mul_ps( dabc[1], coefs[n]));
		    mgrad[3*n+2] = _mm_add_ps (mgrad[3*n+2], _mm_mul_ps( dabc[2], coefs[n]));
		}
	    }
  
    /* Now, store results back */
    float dxInv = spline->x_grid.delta_inv;
    float dyInv = spline->y_grid.delta_inv;
    float dzInv = spline->z_grid.delta_inv;
    for (int n=0; n<N; n++) {
	int nd4 = n>>2;
	int nm4 = n & 3;
	vals[n] = ((float*)mvals)[n];
	grads[3*n+0] = ((float*)mgrad)[nd4*12 + 4*0 + nm4] * dxInv;
	grads[3*n+1] = ((float*)mgrad)[nd4*12 + 4*1 + nm4] * dyInv;
	grads[3*n+2] = ((float*)mgrad)[nd4*12 + 4*2 + nm4] * dzInv;
    }
}
#endif

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

