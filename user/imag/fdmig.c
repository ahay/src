/* Constant-velocity implicit finite-difference migration */
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
/*^*/

#include "fdmig.h"
#include "ctridiagonal.h"

static float alpha, dw, beta;
static sf_complex *cd, *cu;
static ctris slv;
static bool hi;
static int nx, nz, nw;

void fdmig_init(bool hi1    /* if use 45-degree */, 
		int nx1     /* lateral samples */, 
		int nz1     /* vertical samples */, 
		int nw1     /* frequency samples */, 
		float dx    /* lateral sampling */, 
		float dz    /* vertical sampling */, 
		float dw1   /* frequency sampling */, 
		float vel   /* velocity */, 
		float beta1 /* "one-sixth" term */)
/*< Initialize >*/
{
    hi = hi1;
    nx = nx1;
    nz = nz1;
    nw = nw1;
    beta = beta1;
    dw = 2.*SF_PI*dw1*dz; /* dimensionless */
    
    alpha = vel*dz/(2*dx); /* dimensionless */
    alpha *= alpha;
    slv = ctridiagonal_init (nx);    
    cd = sf_complexalloc(nx);
    cu = sf_complexalloc(nx);
}

void fdmig_close(void) 
/*< Free allocated storage >*/
{
    ctridiagonal_close(slv);
    free (cd);
    free (cu);
}

void fdmig (sf_complex **dat /* input data */, 
	    float **img         /* output image */, 
	    sf_file movie       /* save movie (if not NULL) */)
/*< Migrate >*/
{
    int iz, ix, iw;
    float omega;
    sf_complex aa, bb, shift; 

    for (iz=0; iz < nz; iz++) { 
	for (ix=0; ix < nx; ix++) {
	    img[iz][ix] = 0.;
	    cu[ix] = sf_cmplx(0.,0.);
	}
	if (NULL != movie) sf_complexwrite(cu,nx,movie);
    }

    
    for (iw=1; iw < nw; iw++) { 
	omega = dw*iw;
	aa = sf_cmplx(beta*omega,alpha);
#ifdef SF_HAS_COMPLEX_H
	if (hi) aa += alpha/omega; /* add the third derivative term */	
	bb = omega - 2.*aa;
#else
	if (hi) aa.r += alpha/omega; /* add the third derivative term */
	bb.r = omega - 2.*aa.r;
	bb.i = - 2.*aa.i;
#endif

	ctridiagonal_const_define (slv,conjf(bb),conjf(aa));
	    
	for (ix=0; ix < nx; ix++) {
	    cu[ix] = dat[iw][ix];
	}

	for (iz=0; iz < nz; iz++) {
	    cd[0] = sf_cmplx(0.,0.);      
	    for (ix=1; ix < nx-1; ix++) {
#ifdef SF_HAS_COMPLEX_H
		cd[ix] = aa*(cu[ix+1] + cu[ix-1]) + bb*cu[ix];
#else
		cd[ix] = sf_cadd(sf_cmul(aa,sf_cadd(cu[ix+1],cu[ix-1])),
				 sf_cmul(bb,cu[ix]));
#endif
	    }
	    cd[nx-1] = sf_cmplx(0.,0.);
	    
	    ctridiagonal_solve (slv, cd);
	    
	    for (ix=0; ix < nx; ix++) {
		shift = sf_cmplx(cosf(omega),sinf(omega));
#ifdef SF_HAS_COMPLEX_H
		cu[ix] = cd[ix]*shift;
#else
		cu[ix] = sf_cmul(cd[ix],shift);
#endif
		img[iz][ix] += crealf(cu[ix]);
	    }

	    if (NULL != movie) sf_complexwrite(cu,nx,movie);
	}
    }
}

/* 	$Id$	 */

