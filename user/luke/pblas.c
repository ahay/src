/* Supplying BLAS interface */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif


/*^*/

/* using dumb implementations */

void pblas_saxpy(int n, float a, const float *x, int sx, float *y, int sy)
/*< y += a*x >*/
{
    int i, ix, iy;
#ifdef _OPENMP
#pragma omp parallel for private(ix,iy)
#endif
    for (i=0; i < n; i++) {
	ix = i*sx;
	iy = i*sy;
	y[iy] += a * x[ix];
    }
}

void pblas_saxpyA(int n, float a, const float *x, int sx, float *y, int sy)
/*< y += a*x >*/
{
    int i, ix, iy;
#ifdef _OPENMP
#pragma omp for private(ix,iy)
#endif
    for (i=0; i < n; i++) {
	ix = i*sx;
	iy = i*sy;
	y[iy] += a * x[ix];
    }
}

void pblas_sswap(int n, float *x, int sx, float* y, int sy) 
/*< swap x and y >*/
{
    int i, ix, iy;
    float t;
#ifdef _OPENMP
#pragma omp parallel for private(ix,iy,t)
#endif
    for (i=0; i < n; i++) {
	ix = i*sx;
	iy = i*sy;
	t = x[ix];
	x[ix] = y[iy];
	y[iy] = t;
    }
}

void pblas_sswapA(int n, float *x, int sx, float* y, int sy) 
/*< swap x and y >*/
{
    int i, ix, iy;
    float t;
#ifdef _OPENMP
#pragma omp for private(ix,iy,t)
#endif
    for (i=0; i < n; i++) {
	ix = i*sx;
	iy = i*sy;
	t = x[ix];
	x[ix] = y[iy];
	y[iy] = t;
    }
}

float pblas_sdot(int n, const float *x, int sx, const float *y, int sy)
/*< x'y float -> complex >*/
{
    int i, ix, iy;
    float dot=0;

#ifdef _OPENMP
#pragma omp parallel for  private (ix,iy) reduction(+:dot)
#endif
    for (i=0; i < n; i++) {
	ix = i*sx;
	iy = i*sy;
	dot += x[ix] * y[iy];
    }

    return dot;
}


float  pblas_dsdot(int n, const float *x, int sx, const float *y, int sy)
/*< x'y float -> complex >*/
{
    int i, ix, iy;
    double dot=0;
#ifdef _OPENMP
#pragma omp parallel for private(ix,iy) reduction(+:dot) 
#endif
   for (i=0; i < n; i++) {
	ix = i*sx;
	iy = i*sy;
	dot += (double) x[ix] * y[iy];
    }
    return dot;

}

float pblas_snrm2 (int n, const float* x, int sx) 
/*< sum x_i^2 >*/
{
    int i, ix;
    float xn;

    xn = 0.0;
#ifdef _OPENMP
#pragma omp parallel for private(ix) reduction(+:xn)
#endif
    for (i=0; i < n; i++) {
	ix = i*sx;
	xn += x[ix]*x[ix];
    }
    return xn;
}

float pblas_scnrm2 (int n, const void* x, int sx) 
/*< sum |x_i|^2 >*/
{
    int i, ix;
    float xn;
    const sf_complex* c;

    c = (const sf_complex*) x;

    xn = 0.0;
#ifdef _OPENMP
#pragma omp parallel for private(ix) reduction(+:xn)
#endif
    for (i=0; i < n; i++) {
	ix = i*sx;
#ifdef SF_HAS_COMPLEX_H
	xn += crealf(c[ix]*conjf(c[ix]));
#else
	xn += crealf(sf_cmul(c[ix],conjf(c[ix])));
#endif
    }
    return xn;
}

void pblas_sscal(int n, float alpha, float *x, int sx)
/*< x = alpha*x >*/
{
    int i, ix;
#ifdef _OPENMP
#pragma omp parallel for private(ix)
#endif
    for (i=0; i < n; i++) {
        ix = i*sx;
	x[ix] *= alpha;
    }
}

void pblas_csscal(int n, float alpha, void *x, int sx)
/*< x = alpha*x >*/
{
    int i, ix;
    sf_complex* c;

    c = (sf_complex*) x;
#ifdef _OPENMP
#pragma omp parallel for private(ix)
#endif
    for (i=0; i < n; i++) {
        ix = i*sx;
#ifdef SF_HAS_COMPLEX_H
	c[ix] *= alpha;
#else
	c[ix] = sf_crmul(c[ix],alpha);
#endif
    }
}


void pblas_cdotc_sub(int n, 
		     const void *x, int sx,
		     const void *y, int sy, void *dot)
/*< complex hermitian dot product >*/
{
    sf_complex *cx, *cy, xi, yi;
    sf_double_complex prod, pi;
    int i, ix, iy;
    
    cx = (sf_complex*) x;
    cy = (sf_complex*) y;

    prod = sf_dcmplx(0.,0.);
#ifdef _OPENMP
#pragma omp parallel for private(ix,iy,xi,yi,pi) reduction(+:prod)
#endif
    for (i=0; i < n; i++) {
	ix = i*sx;
	iy = i*sy;
	xi = cx[ix];
	yi = cy[iy];
	pi = sf_dcmplx((double) crealf(xi)*crealf(yi) + 
		       (double) cimagf(xi)*cimagf(yi), 
		       (double) crealf(xi)*cimagf(yi) - 
		       (double) cimagf(xi)*crealf(yi));
#ifdef SF_HAS_COMPLEX_H
	prod += pi;
#else
	prod = sf_dcadd(prod,pi);
#endif
    }

    *((sf_complex *) dot) = sf_cmplx(creal(prod),cimag(prod));
}

