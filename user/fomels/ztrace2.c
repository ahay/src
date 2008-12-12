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

#ifndef _ztrace_h

#define NS 4 /* number of outputs */
/*^*/

#endif

static int nx, nz, na, nax;
static float dx,dz,da, x0,z0,a0, xm,zm,am, r2a;
static sf_eno2 cvel;
static enogrid slice;
static grid1 *prev, *curr;

static int snap(float *f, int n);

void ztrace2_init (int order        /* interpolation order for velocity */, 
		   int iorder       /* interpolation order for values */,
		   float** vel      /* slowness [nx][nz] */,
		   int nx1          /* lateral samples */, 
		   int nz1          /* depth samples */, 
		   int na1          /* angle samples */, 
		   float dx1        /* lateral sampling */, 
		   float dz1        /* depth sampling */, 
		   float da1        /* angle sampling */,
		   float x01        /* lateral origin */, 
		   float z01        /* depth origin */, 
		   float a01        /* angle origin */)
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
	    
	    f[0] = 0.;
	    f[1] = x0+kx*dx; 
	    f[2] = z0;
	    f[3] = a*r2a;
	    
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
    float v, g[2], t, z1, x1, x2, a1, a2, delz, tn, cs, f[4], s, fx;
    int kx, ka, k, ix;

    /* assign the previous slice for interpolation */
    for (ix=0; ix < nx; ix++) {
	grid1_close(prev[ix]);
	prev[ix] = curr[ix];
	curr[ix] = grid1_init();
    }
    
    for (kx=0; kx < nx; kx++) { /* loop in x */

	/* get slowness squared and gradient */
	sf_eno2_apply(cvel,kz,kx,0.,0.,&v,g,BOTH);
	g[0] /= dz;
	g[1] /= dx;
	
	for (ka=0; ka < na; ka++) { /* loop in angle */
	    k = ka + kx*na; /* index in a slice */
	    a1 = a0+ka*da;

	    tn = tanf(a1);
	    cs = cosf(a1);

	    x1 = x0+kx*dx; 
	    z1 = z0+kz*dz; 

	    /* decide if we are out already */
	    if ((kz == 0    && cs > 0.) ||
		(kz == nz-1 && cs < 0.) ||
		(kx == 0    && tn < 0.) ||
		(kx == nx-1 && tn > 0.)) {
		f[0] = 0.;
		f[1] = x1;
		f[2] = z1;
		f[3] = a1*r2a;
	    } else {
		x2 = x1+tn*dz;
		if (x2 < x0) {
		    /* exit from the left side */
		    delz = (x0-x1)/tn;
		    t = s*delz/cs;
		    a2 = a1 + (g[1]-g[0]*tn)*delz/v;
		    
		    enogrid_apply(slice,0,0.,a2,f);
		    f[1]=x0;
		    f[2]=z1-delz;
		} else if (x2 > xm) {
		    /* exit from the right side */
		    delz = (xm-x1)/tn;
		    t = s*delz/cs;
		    a2 = a1 + (g[1]-g[0]*tn)*delz/v;
		    
		    enogrid_apply(slice,nx-1,0.,a2,f);
		    f[1]=xm;
		    f[2]=z1-delz;
		} else {
		    /* exit on previous level */
		    t = s*dz/cs;
		    a2 = a1 + (g[1]-g[0]*tn)*dz/v;
		    
		    fx = (x2-x0)/dx;
		    ix = snap(&x2,nx);
		    enogrid_apply(slice,ix,fx,a2,f);
		    
		    if (f[1] < x0) f[1]=x0;
		    if (f[1] > xm) f[1]=xm;
		    if (f[2] < z0) f[2]=z0;
		    if (f[2] > zm) f[2]=zm;
		} 
		
		if (f[0] < 0.) f[0]=0.;
		f[0] += t;
		if (f[3] < a0*r2a) f[3]=a0*r2a;
		if (f[3] > am*r2a) f[3]=am*r2a;
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
