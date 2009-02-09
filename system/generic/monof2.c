/* Fit a Gaussian to 2-D spectrum by nonlinear weighted least squares */
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
#include <math.h>

#include <rsf.h>

#include "monof2.h"

/* formulas for separable least squares:

a = f.e/(e.e + eps)
ap = (f.ep-2*a*e.ep)/(e.e + eps)
num = a*(a*e.ep + ap*(e.e + 2.*eps))
den = ap*ap*e.e + 2*a*ap*e.ep + a*a*ep.ep
da = num/den

need: 

fe -> f.e
ee -> e.e + eps
fep -> f.ep
eep -> e.ep
epep -> ep.ep

*/

void monof2(float **data               /* input [ny][nx] */,
	    float **pred               /* prediction [ny][nx] */,
	    int nliter                 /* number of reweighting iterations */,
	    int niter                  /* number of iterations */, 
	    float* a                   /* estimated parameters [3] */, 
	    int nx, float dx, float x0 /* first axis */, 
	    int ny, float dy, float y0 /* second axis */, 
	    bool verb                  /* verbosity flag */)
/*< Estimate shape. (Caution: data gets corrupted) >*/ 
{
    int ix, iy, iter, i, j;
    float f2, fe, ee, fep[3], eep[3], epep[6], x, y, da[3], aa, f, e2;
    float ep[3], eps, num[3], den[6], r2, e, x2, y2, xy, det, w;
    
    eps = 10.*FLT_EPSILON;
    eps *= eps;
    
    if (verb) sf_warning("got a0=%g b0=%g c0=%g\n"
			 "nx=%d dx=%g x0=%g ny=%d dy=%g y0=%g",
			 a[0],a[1],a[2],nx,dx,x0,ny,dy,y0);

    if (nliter > 1) sf_irls_init(nx*ny);
	
    for (i = 0; i < nliter; i++) {
	f2 = 0.;
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		f = data[iy][ix];
		w = (0 == i)? 1.:pred[iy][ix];

		f2 += w*f*f;		    
	    }
	}
	aa = 1.;

	for (iter = 0; iter < niter; iter++) {
	    ee = eps;
	    fe = 0.;
	    for (j=0; j < 3; j++) {	    
		fep[j] = 0.;
		eep[j] = 0.;
		epep[j] = 0.;
		epep[3+j] = 0.;
	    }
	    for (iy=0; iy < ny; iy++) {
		y = y0+iy*dy;
		y2 = y*y;
		for (ix=0; ix < nx; ix++) {
		    x = x0+ix*dx;
		    x2 = x*x;
		    xy = x*y;
		    e = expf(-a[0]*x2-a[1]*xy-a[2]*y2);
		    e2 = e*e;
		    ep[0] = -x2*e;
		    ep[1] = -xy*e;
		    ep[2] = -y2*e;
	    
		    f = data[iy][ix];
		    w = (0 == i)? 1.:pred[iy][ix];
	    
		    ee += w*e2;
		    fe += w*f*e;
		    for (j=0; j < 3; j++) {
			fep[j] += w*f*ep[j];
			eep[j] += w*e*ep[j];
		    }
		    epep[0] += w*ep[0]*ep[0];
		    epep[1] += w*ep[0]*ep[1];
		    epep[2] += w*ep[0]*ep[2];
		    epep[3] += w*ep[1]*ep[1];
		    epep[4] += w*ep[1]*ep[2];
		    epep[5] += w*ep[2]*ep[2];
		}
	    }
              
	    aa = fe/ee;
	    for (j=0; j < 3; j++) {
		da[j] = (fep[j] - 2.*aa*eep[j])/ee;
		num[j] = aa*(aa*eep[j] + da[j]*(ee+eps));
	    }
	    den[0] = aa*aa*epep[0] + da[0]*(2.*aa*eep[0] + da[0]*(ee-eps));
	    den[1] = aa*(aa*epep[1] + da[1]*eep[0]) + 
		da[0]*(aa*eep[1] + da[1]*(ee-eps));
	    den[2] = aa*(aa*epep[2] + da[2]*eep[0]) + 
		da[0]*(aa*eep[2] + da[2]*(ee-eps));
	    den[3] = aa*aa*epep[3] + da[1]*(2.*aa*eep[1] + da[1]*(ee-eps));
	    den[4] = aa*(aa*epep[4] + da[1]*eep[2]) + 
		da[2]*(aa*eep[1] + da[1]*(ee-eps));
	    den[5] = aa*aa*epep[5] + da[2]*(2.*aa*eep[2] + da[2]*(ee-eps));
	    
	    det = 
		den[2]*(den[2]*den[3] - den[1]*den[4]) + 
		den[1]*(den[1]*den[5] - den[2]*den[4]) + 
		den[0]*(den[4]*den[4] - den[3]*den[5]);
	    
	    if (det > 0. && eps > det) {
		det = eps;
	    } else if (det < 0. && -eps < det) {
		det = -eps;
	    }
	
	    da[0] = (
		num[2]*(den[2]*den[3] - den[1]*den[4]) + 
		num[1]*(den[1]*den[5] - den[2]*den[4]) + 
		num[0]*(den[4]*den[4] - den[3]*den[5])
		)/det;
	    da[1] = ( 
		num[1]*(den[2]*den[2] - den[0]*den[5]) + 
		num[0]*(den[1]*den[5] - den[2]*den[4]) + 
		num[2]*(den[4]*den[0] - den[1]*den[2])
		)/det;
	    da[2] = (
		num[0]*(den[2]*den[3] - den[1]*den[4]) + 
		num[2]*(den[1]*den[1] - den[3]*den[0]) + 
		num[1]*(den[4]*den[0] - den[1]*den[2])
		)/det;
	    
	    r2 = f2 - aa*aa*(ee+eps); /* residual squared */
	    if (verb) sf_warning("monof: iter=%d r2=%g da=(%g,%g,%g) "
				 "aa=%g a=(%g,%g,%g)",
				 iter,r2,da[0],da[1],da[2],aa,a[0],a[1],a[2]);
	    
	    for (j=0; j < 3; j++) {
		a[j] += da[j];
	    }
	    if (a[0] < 0.) a[0] = 0.;
	    if (a[2] < 0.) a[2] = 0.;
	    
	    if (r2 < eps || da[0]*da[0] + da[1]*da[1] + da[2]*da[2] < eps) 
		break;
	} /* iter */
    
	for (iy=0; iy < ny; iy++) {
	    y = y0+iy*dy;
	    y2 = y*y;
	    for (ix=0; ix < nx; ix++) {
		x = x0+ix*dx;
		x2 = x*x;
		xy = x*y;
		e = expf(-a[0]*x2-a[1]*xy-a[2]*y2);

		pred[iy][ix] = (i < nliter-1)? aa*e-data[iy][ix]: aa*e;
	    }
	}

	if (i < nliter-1) sf_cauchy (nx*ny,pred[0],pred[0]);
	
	if (verb) sf_warning ("%d iterations", iter);
    } /* reweighting iterations */

    if (nliter > 1) sf_irls_close();
}

/* 	$Id$	 */
