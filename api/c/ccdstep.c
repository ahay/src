/* Conjugate-direction iteration for complex numbers. */
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

#include <float.h>

#include "cdstep.h"
#include "clist.h"
#include "alloc.h"
#include "blas.h"

#include "_bool.h"
#include "c99.h"
#include "komplex.h"
/*^*/

static sf_clist steps;

static void mysaxpy(int n, sf_double_complex a, 
		    const sf_complex *x, 
		    sf_complex *y);

static sf_double_complex mydsdot(int n, 
				 const sf_complex *cx, 
				 const sf_complex *cy);

void sf_ccdstep_init(void) 
/*< initialize internal storage >*/
{
    steps = sf_clist_init();
}

void sf_ccdstep_close(void) 
/*< free internal storage >*/
{
    sf_clist_close(steps);
}

void sf_ccdstep(bool forget          /* restart flag */, 
	       int nx               /* model size */, 
	       int ny               /* data size */, 
	       sf_complex* x        /* current model [nx] */, 
	       const sf_complex* g  /* gradient [nx] */, 
	       sf_complex* rr       /* data residual [ny] */, 
	       const sf_complex* gg /* conjugate gradient [ny] */) 
/*< Step of conjugate-direction iteration. 
  The data residual is rr = A x - dat
>*/
{
    sf_complex *s, *si, *ss;
    sf_double_complex alpha;
    double beta;
    int i, n, ix, iy;

    s = sf_complexalloc(nx+ny);
    ss = s+nx;

    for (ix=0; ix < nx; ix++) {  s[ix] = g[ix]; }
    for (iy=0; iy < ny; iy++) { ss[iy] = gg[iy]; }

    sf_clist_rewind (steps);
    n = sf_clist_depth (steps);

    for (i=0; i < n; i++) {
	sf_clist_down (steps, &si, &beta);
#ifdef SF_HAS_COMPLEX_H
	alpha = - mydsdot(ny, gg, si+nx) / beta;
#else
	alpha = sf_dcrmul(mydsdot(ny, gg, si+nx),-1./beta);
#endif
	mysaxpy(nx+ny,alpha,si,s);
    }
    
    beta = creal(mydsdot(ny, s+nx, s+nx));
    if (beta < DBL_EPSILON) return;

    sf_clist_add (steps, s, beta);
    if (forget) sf_clist_chop (steps);
#ifdef SF_HAS_COMPLEX_H
    alpha = -  mydsdot(ny, rr, ss) / beta;
#else
    alpha = sf_dcrmul(mydsdot(ny, rr, ss),-1./beta);
#endif    

    mysaxpy(nx,alpha,s,x);
    mysaxpy(ny,alpha,ss,rr);
}

static void mysaxpy(int n, sf_double_complex a, 
		    const sf_complex *x, 
		    sf_complex *y)
/* y += a*x */
{
    int i;

#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
    for (i=0; i < n; i++) {
#ifdef SF_HAS_COMPLEX_H
	y[i] += a * x[i];
#else
	y[i] = sf_cadd(y[i],sf_cmul(x[i],sf_cmplx(sf_creal(a),sf_cimag(a))));
#endif
    }
}

static sf_double_complex mydsdot(int n, 
				 const sf_complex *cx, 
				 const sf_complex *cy)
/* Hermitian dot product */
{
    sf_complex xi, yi;
    sf_double_complex prod, pi;
    int i;
    
    prod = sf_dcmplx(0.,0.);
#if defined(_OPENMP) && defined(SF_HAS_COMPLEX_H)
#pragma omp parallel for private(i,xi,yi,pi) reduction(+:prod)
#endif
    for (i=0; i < n; i++) {
	xi = cx[i];
	yi = cy[i];
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

    return prod;
}


