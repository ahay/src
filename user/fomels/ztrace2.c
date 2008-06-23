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

#include "enogrid.h"
#include "grid1.h"
#include "analytical.h"

#ifndef _ztrace_h

#define NS 4 /* number of outputs */
/*^*/

#endif

static int nx, nz, na, nax;
static float dx,dz,da, x0,z0,a0, xm,zm,am;
static sf_eno2 cvel;
static enogrid slice;
static grid1 *prev, *curr;

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

    for (kx=0; kx < nx; kx++) { /* loop in x */
	for (ka=0; ka < na; ka++) { /* loop in angle */
	    a = a0+ka*da;
	    
	    f[0] = 0.;
	    f[1] = x0+kx*dx; 
	    f[2] = z0;
	    f[3] = a*180./SF_PI;
	    
	    grid1_insert(curr[kx],a,4,f);
	}
    } 

    slice = enogrid_init (iorder, nx, NS, prev);
}

void ztrace2_close (void)
/*< Free allocated storage >*/
{
    sf_eno2_close (cvel);
    enogrid_close(slice);

    /* close grid */
    free(prev);
    free(curr);
}

void ztrace2_write (sf_file out)
/*< Write it out >*/
{
    int ix;
    for (ix=0; ix < nx; ix++) {
	grid1_write(curr[ix],NS,out);
    }
}

void ztrace2_step (int kz) 
/*< Step in depth >*/
{
    float v2, v1, g1[2], g2[2], t, z1, x1, z2, x2, a1, a2, p2[2], p1[2], f[4];
    float s, sx, sz, sx1, sx2, sz1, sz2, fx, fz;
    int kx, ka, k, ix, iz;

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

	    a1 = a0+ka*da;
	    p1[1] = sinf(a1);
	    p1[0] = -cosf(a1);
	    /* p0 is dimensionless */

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
		f[3] = a1*180./SF_PI;
		grid1_insert(curr[kx],a1,4,f);
		continue;
	    } 

	    /* find the nearest intersection of ray and box */
	    sx1 = sf_quadratic_solve (g1[1],p1[1],2*(x1-x0)/v1);
	    sx2 = sf_quadratic_solve (g1[1],p1[1],2*(x1-xm)/v1);
	    sz1 = sf_quadratic_solve (g1[0],p1[0],2*dz/v1);
	    sz2 = sf_quadratic_solve (g1[0],p1[0],2*(z1-zm)/v1);
	    
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
		z2 = z1 + (p1[0]+v1*g1[0]*s*0.5)*s;
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
		x2 = x1 + (p1[1]+v1*g1[0]*s*0.5)*s;
		fx = (x2-x0)/dx;
		ix = snap(&fx,nx);
	    }

	    /* slowness and gradient at new location */
	    sf_eno2_apply(cvel,iz,ix,fz,fx,&v2,g2,BOTH);
	    g2[1] /= dx;
	    g2[0] /= dz;

	    t = analytical(x1,z1,v1,g1,p1,
			   x2,z2,v2,g2,p2);

	    if (s == sz1) { /* to previous level */
		/* interpolate */
		a1 = (sf_cell_p2a(p1)-a0)/da;
		a2 = (sf_cell_p2a(p2)-a0)/da;
		enogrid_apply(slice,ix,fx,a2,f);
		if (f[0] < 0.) f[0]=0.;
		f[0] += t;
		if (f[1] < x0) f[1]=x0;
		if (f[1] > xm) f[1]=xm;
		if (f[2] < z0) f[2]=z0;
		if (f[2] > zm) f[2]=zm;
		if (f[3] < a0*180./SF_PI) f[3]=a0*180./SF_PI;
		if (f[3] > am*180./SF_PI) f[3]=am*180./SF_PI;
	    } else {
		f[0] = t;
		f[1] = x2;
		f[2] = z2;
		f[3] = a2*180./SF_PI;
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
