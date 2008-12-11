/* Multiple arrivals by marching down using Hamiltonian dynamics. */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

#include "hdtrace.h"

#include "enogrid.h"
#include "grid1.h"
#include "analytical.h"

#ifndef _hdtrace_h

#define NS 4 /* number of outputs */
/*^*/

#endif

static int nx, nz, na, nax;
static float dx,dz,da, x0,z0,a0, xm,zm,am, r2a;
static sf_eno2 cvel;
static enogrid slice;
static grid1 *prev, *curr;

static int snap(float *f, int n);

void hdtrace_init (int order        /* interpolation order for velocity */, 
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
		   float** vel      /* slowness [nx][nz] */)
/*< Initialize >*/
{
    int ix, kx, ka;
    float a, f[NS];

    nx = nx1; nz = nz1; na = na1;
    dx = dx1; dz = dz1; da = da1;
    x0 = x01; z0 = z01; a0 = a01;
    nax = na*nx; 
    xm = x0 + (nx-1)*dx;
    zm = z0 + (nz-1)*dz;
    am = a0 + (na-1)*da;

    cvel = sf_eno2_init (order, nz, nx);
    sf_eno2_set (cvel, vel);

    prev = (grid1*) sf_alloc(nx,sizeof(grid1));
    curr = (grid1*) sf_alloc(nx,sizeof(grid1));
    for (ix=0; ix < nx; ix++) {
	prev[ix] = grid1_init();
	curr[ix] = grid1_init();
    }

    r2a = 180./SF_PI; /* radian to degree */

    for (kx=0; kx < nx; kx++) { /* loop in x */
	for (ka=0; ka < na; ka++) { /* loop in angle */
	    a = a0+ka*da;

            /* f is vector (time, x, z, angle) for one way dynamic rays */
	    f[0] = 0.;
	    f[1] = x0+kx*dx; 
	    f[2] = z0;
	    f[3] = a*r2a;
	    
	    grid1_insert(curr[kx],a,4,f);
	}
    } 

    slice = enogrid_init (iorder, nx, NS, prev);
}

void hdtrace_close (void)
/*< Free allocated storage >*/
{
    sf_eno2_close (cvel);
    enogrid_close(slice);

    /* close grid */
    free(prev);
    free(curr);
}

void hdtrace_write (sf_file out)
/*< Write it out >*/
{
    int ix;
    for (ix=0; ix < nx; ix++) {
	grid1_write(curr[ix],NS,out);
    }
}

void hdtrace_step (int kz) 
/*< Step in depth >*/
{
    float v2, v1, g1[2], g2[2], t, z1, x1, z2, x2, a1, a2, p2[2], p1[2], f[4];
    float s, sx, sz, sx1, sx2, sz1, sz2, fx, fz, fa, stepz;
    int kx, ka, k, ix, iz, ia;

    /* assign the previous slice for interpolation */
    for (ix=0; ix < nx; ix++) {
	grid1_close(prev[ix]);
	prev[ix] = curr[ix];
	curr[ix] = grid1_init();
    }
    
    for (kx=0; kx < nx; kx++) { /* loop in x */

	/* get slowness and gradient */
	sf_eno2_apply(cvel,kz,kx,0.,0.,&v1,g1,BOTH);
	g1[0] /= dz;
	g1[1] /= dx;
	
	for (ka=0; ka < na; ka++) { /* loop in angle */
	    k = ka + kx*na; /* index in a slice */

            /* initial angle */
	    a1 = a0+ka*da;

	    /* p1 is dimensionless */
	    p1[0] = -cosf(a1); /* pz */
	    p1[1] =  sinf(a1); /* px */

            /* initial position */
	    x1 = x0+kx*dx; 
	    z1 = z0+kz*dz; 

	    /* decide if we are out already */
	    if ((kz == 0    && p1[0] < 0.) ||
		(kz == nz-1 && p1[0] > 0.) ||
		(kx == 0    && p1[1] < 0.) ||
		(kx == nx-1 && p1[1] > 0.)) {
		f[0] = 0.;
		f[1] = x1;
		f[2] = z1;
		f[3] = a1*r2a;
		grid1_insert(curr[kx],a1,4,f);
		continue;
	    } 

	    /* find the nearest intersection of ray and box */
	    /*
	    sx1 = sf_quadratic_solve (g1[1],p1[1],2*(x1-x0)/v1);
	    sx2 = sf_quadratic_solve (g1[1],p1[1],2*(x1-xm)/v1);
	    sz1 = sf_quadratic_solve (g1[0],p1[0],2*dz/v1);
	    sz2 = sf_quadratic_solve (g1[0],p1[0],2*(z1-zm)/v1);
	    */



	    sx1 = (x0-x1)/p1[1]; if (sx1 <= 0.) sx1 = SF_HUGE;
	    sx2 = (xm-x1)/p1[1]; if (sx2 <= 0.) sx2 = SF_HUGE;
	    sz1 = -dz/p1[0];     if (sz1 <= 0.) sz1 = SF_HUGE;
	    sz2 = (zm-z1)/p1[0]; if (sz2 <= 0.) sz2 = SF_HUGE;

	    sx = SF_MIN(sx1,sx2);
	    sz = SF_MIN(sz1,sz2);
	    s = SF_MIN(sx,sz);

	    if (s == sx) { /* exited from the side */
		if (s == sx1) {
		    ix = 0;
		    x2 = x0;
		} else {
		    ix = nx-1;
		    x2 = xm;
		}
		fx = 0.;
		z2 = z1 + p1[0]*s;
		fz = (z2-z0)/dz;
		iz = snap(&fz,nz);
	    } else {
		if (s == sz1) {
		    iz = kz-1;
		    z2 = z1 - dz;
		} else {
		    iz = nz-1;
		    z2 = zm;
		}
		fz = 0.;
		stepz = dz/p1[0];
                /* integrate dx/dz */
		x2 = x1 + p1[1]*stepz;
		fx = (x2-x0)/dx;
		ix = snap(&fx,nx);
                
                /* integrate dpx/dz */
		a2 = a1 + g1[1]*stepz;
		fa = (a2-a0)/da;
		ia = snap(&fa,na);
                
                /* integrate dT/dz */
		t = -v1*stepz;

                /* Neri and Candy 4th order symplectic algorithm */
		/* traveltime integration */
                /* t = analytical */
	    }

            /* convert ray parameter to angle (filt/lib/cell.c) */
	    a1 = sf_cell_p2a(p1);
	    a2 = sf_cell_p2a(p2);

	    if (s == sz1) { /* to previous level */
		/* interpolate */
		enogrid_apply(slice,ix,fx,a2,f);
		if (f[0] < 0.) f[0]=0.;
		f[0] += t;
		if (f[1] < x0) f[1]=x0;
		if (f[1] > xm) f[1]=xm;
		if (f[2] < z0) f[2]=z0;
		if (f[2] > zm) f[2]=zm;
		if (f[3] < a0*r2a) f[3]=a0*r2a;
		if (f[3] > am*r2a) f[3]=am*r2a;
	    } else {
		f[0] = t;
		f[1] = x2;
		f[2] = z2;
		f[3] = a2*r2a;
	    }
	    
	    grid1_insert(curr[kx],a1,4,f);
	} /* ka */
    } /* kx */
}

static int snap(float *f, int n)
{
    int i;

    i = floorf(*f);
    if (i < 0) {
	i=0;
	*f=0.;
    } else if (i >= n-1) {
	i=n-1;
	*f=0.;
    } else {
	*f -= i;
    }

    return i;
}
