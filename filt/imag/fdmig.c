#include <rsf.h>

#include "fdmig.h"
#include "ctridiagonal.h"

static float alpha, dw, beta;
static float complex *cd, *cu;
static ctris slv;
static bool hi;
static int nx, nz, nw;

void fdmig_init(bool hi1, int nx1, int nz1, int nw1, 
		float dx, float dz, float dw1, float vel, float beta1)
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
{
    ctridiagonal_close(slv);
    free (cd);
    free (cu);
}

void fdmig (float complex **dat, float **img, sf_file movie)
{
    int iz, ix, iw;
    float omega;
    float complex aa, bb; 

    for (iz=0; iz < nz; iz++) { 
	for (ix=0; ix < nx; ix++) {
	    img[iz][ix] = 0.;
	    cu[ix] = 0.;
	}
	if (NULL != movie) sf_complexwrite(cu,nx,movie);
    }

    
    for (iw=1; iw < nw; iw++) { 
	omega = dw*iw;
	aa = beta*omega + I*alpha;
	if (hi) aa += alpha/omega; /* add the third derivative term */	
	bb = omega - 2.*aa;

	ctridiagonal_const_define (slv,conjf(bb),conjf(aa));
	    
	for (ix=0; ix < nx; ix++) {
	    cu[ix] = dat[iw][ix];
	}

	for (iz=0; iz < nz; iz++) {
	    cd[0] = 0.;      
	    for (ix=1; ix < nx-1; ix++) {
		cd[ix] = aa*(cu[ix+1] + cu[ix-1]) + bb*cu[ix];
	    }
	    cd[nx-1] = 0.;
	    
	    ctridiagonal_solve (slv, cd);
	    
	    for (ix=0; ix < nx; ix++) {
		cu[ix] = cd[ix]*cexpf(I*omega);
		img[iz][ix] += crealf(cu[ix]);
	    }

	    if (NULL != movie) sf_complexwrite(cu,nx,movie);
	}
    }
}
