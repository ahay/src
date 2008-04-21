/* Multiple arrivals by marching down. */
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

#include <math.h>

#include <rsf.h>

#include "ztrace2.h"

#ifndef _ztrace_h

#define NS 4 /* number of outputs */
/*^*/

#endif

static int nx, nz, na, nax;
static float dx,dz,da, x0,z0,a0, x1,z1;
static sf_eno2 cvel, fslc[NS];
static float **slice;

static int snap(float *f, int n);

void ztrace2_init (int order        /* interpolation order for velocity */, 
		   int iorder       /* interpolation order for values */,
		   int nx1          /* lateral samples */, 
		   int nz1          /* depth samples */, 
		   int na1          /* angle samples */, 
		   float dx1        /* lateral sampling */, 
		   float dz1        /* depth sampling */, 
		   float da1        /* angle sampling */,
		   float x01        /* lateral origin */, 
		   float z01        /* depth origin */, 
		   float a01        /* angle origin */,
		   float** vel      /* slowness [nx][nz] */, 
		   float** slice_in /* depth slice [nx][na][NS] */)
/*< Initialize >*/
{
    int is;

    slice = slice_in;
    nx = nx1; nz = nz1; na = na1;
    dx = dx1; dz = dz1; da = da1;
    x0 = x01; z0 = z01; a0 = a01;
    nax = na*nx; 
    x1 = x0 + (nx-1)*dx;
    z1 = z0 + (nz-1)*dz;

    cvel = sf_eno2_init (order, nz, nx);
    sf_eno2_set (cvel, vel);

    for (is=0; is < NS; is++) {
	fslc[is] = sf_eno2_init (iorder, na, nx);
    }
}

void ztrace_close (void)
/*< Free allocated storage >*/
{
    int is;
    
    sf_eno2_close (cvel);
    
    for (is = 0; is < NS; is++) {
	sf_eno2_close (fslc[is]);
    }
}

void ztrace2_step (int kz) 
/*< Step in depth >*/
{
    float v, v0, g[2], g0[2], t, z, x, a, p[2], s, sx, sz, sx1, sx2, sz1, sz2, fx, fz, fa;
    int is, kx, kp, k, ix, iz, ia;

    /* assign the previous slice for interpolation */
    for (is=0; is < NS; is++) {
	sf_eno2_set1 (fslc[is], slice[is]);
    }
    
    for (kx=0; kx < nx; kx++) { /* loop in x */

	/* get slowness and gradient */
	sf_eno2_apply(cvel,kz,kx,0.,0.,&v0,g0,BOTH);
	g0[1] /= dx;
	g0[0] /= dz;

	x = x0+kx*dx; ix=kx;
	z = z0+kz*dz; iz=kz;
	
	for (kp=0; kp < na; kp++) { /* loop in angle */
	    k = kp + kx*na; /* index in a slice */

	    a = a0+kp*da;
	    p[1] = sinf(a);
	    p[0] = -cosf(a);

            t = sf_cell_update2 (2, 0.,v0, p, g0);
            /* p is normal vector now ||p|| = 1 */

	    /* decide if we are out already */
	    if ((iz == 0    && p[0] < 0.) ||
		(iz == nz-1 && p[0] > 0.) ||
		(ix == 0    && p[1] < 0.) ||
		(ix == nx-1 && p[1] > 0.)) {
		slice[0][k] = 0.;
		slice[1][k] = x;
		slice[2][k] = z;
		slice[3][k] = a*180./SF_PI;
		continue;
	    } 

	    /* find the nearest intersection of ray and box */
	    sx1 = sf_quadratic_solve (g0[1],p[1],2*(x-x0)/v0);
	    sx2 = sf_quadratic_solve (g0[1],p[1],2*(x-x1)/v0);
	    sz1 = sf_quadratic_solve (g0[0],p[0],-2*dz/v0);
	    sz2 = sf_quadratic_solve (g0[0],p[0],2*(z-z1)/v0);

	    
	    sx = SF_MIN(sx1,sx2);
	    sz = SF_MIN(sz1,sz2);
	    s = SF_MIN(sx,sz);

	    t =	sf_cell_update1 (2, s, v0, p, g0);
	    /* p is slowness vector now ||p||=v */

	    if (s == sx) { /* exited from the side */
		if (s == sx1) {
		    ix = 0;
		    x = x0;
		} else {
		    ix = nx-1;
		    x = x1;
		}
		fx = 0.;
		z += p[0]*s;
		fz = (z-z0)/dz;
		iz = snap(&fz,nz);
	    } else {
		if (s == sz1) {
		    iz = kz-1;
		    z -= dz;
		} else {
		    iz = nz-1;
		    z = z1;
		}
		fz = 0.;
		x += p[1]*s;
		fx = (x-x0)/dx;
		ix = snap(&fx,nx);
	    }

	    /* slowness and gradient at new location */
	    sf_eno2_apply(cvel,iz,ix,fz,fx,&v,g,BOTH);
	    g[1] /= dx;
	    g[0] /= dz;

	    t += sf_cell_update2 (2, s, v, p, g);
	    /* p is normal vector now ||p||=1 */

	    if (s == sz1) { /* to previous level */
		/* interpolate */
		fa = (sf_cell_p2a(p)-a0)/da;
		ia = snap(&fa,na);
		for (is=0; is < NS; is++) {
		    sf_eno2_apply(fslc[is],ia,ix,fa,fx,
				  &slice[is][k],NULL,FUNC);
		}
		slice[0][k] += t;
	    } else {
		slice[0][k] = t;
		slice[1][k] = x;
		slice[2][k] = z;
		slice[3][k] = sf_cell_p2a(p)*180./SF_PI;
	    }
	} /* kp */
    } /* kx */
}

static int snap(float *f, int n)
{
    int i;

    i = floorf(*f);
    if (i < 0) {
	i=0;
	*f=0.;
    } else if (i >= n) {
	i=n-1;
	*f=0.;
    } else {
	*f -= i;
    }

    return i;
}
