/* Multiple arrivals by marching down/up using Hamiltonian dynamics. */
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

#include "enogrid.h"
#include "grid1.h"
#include "symplectrace.h"

#ifndef _hdtrace_h

#define NS 4 /* number of outputs */
/*^*/

#endif


static int nx, nz, np, npx;
static float dx,dz,dp,x0,z0,p0,xm,zm,pm;
static pqv *hvec;
static sf_eno2 cvel;
static enogrid slice;
static grid1 *prev, *curr;
static float** slow /* slowness [nx][nz] */


static int snap(float *f, int n);


void hdtrace_init (int order        /* interpolation order for velocity */, 
		   int iorder       /* interpolation order for values */,
		   int nx1          /* lateral samples */, 
		   int nz1          /* depth samples */, 
		   int np1          /* horizontal slowness samples */, 
		   float dx1        /* lateral sampling */, 
		   float dz1        /* depth sampling */, 
		   float dp1        /* horizontal slowness sampling */,
		   float x01        /* lateral origin */, 
		   float z01        /* depth origin */, 
		   float p01        /* horizontal slowness origin */)
/*< Initialize >*/
{
    int ix, kx, kp;
    float p, f[NS];

    nx = nx1; nz = nz1; np = np1;
    dx = dx1; dz = dz1; dp = dp1;
    x0 = x01; z0 = z01; p0 = p01;
    npx = np*nx; 
    xm = x0 + (nx-1)*dx;
    zm = z0 + (nz-1)*dz;
    pm = p0 + (np-1)*dp;

    cvel = sf_eno2_init (order, nz, nx);
    sf_eno2_set (cvel, slow);

    prev = (grid1*) sf_alloc(nx,sizeof(grid1));
    curr = (grid1*) sf_alloc(nx,sizeof(grid1));
    for (ix=0; ix < nx; ix++) {
	prev[ix] = grid1_init();
	curr[ix] = grid1_init();
    }

    for (kx=0; kx < nx; kx++) { /* loop in x */
	for (kp=0; kp < np; kp++) { /* loop in horizontal slowness */
	    p = p0+kp*dp;

            /* f is vector (traveltime, x, z, p) for one way dynamic rays */
	    f[0] = 0.;
	    f[1] = x0+kx*dx; 
	    f[2] = z0;
	    f[3] = p;
	    
	    grid1_insert(curr[kx],p,4,f);
	}
    } 

    slice = enogrid_init (iorder, nx, NS, prev);
    nc4_init();
    hvec = (pqv*) sf_alloc(1,sizeof(pqv));
}

void hdtrace_close (void)
/*< Free allocated storage >*/
{
    free(hvec);

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
    float x1, z1, p1[2], f[4], ss, ds;
    int kx, kp, k, ix;






    float t;
    float fx, fz, fa, stepz;
    int iz, ia;








    /* grid dimension restrictions */
    /* dx.dp > 1/(4.pi.f) */
    /* dz.dp > 1/(4.pi.f) */

    /* assign the previous slice for interpolation */
    for (ix=0; ix < nx; ix++) {
	grid1_close(prev[ix]);
	prev[ix] = curr[ix];
	curr[ix] = grid1_init();
    }

    for (kx=0; kx < nx; kx++) { /* loop in x */
	
	for (kp=0; kp < np; kp++) { /* loop in horizontal slowness */
	    k = kp + kx*np; /* index in a slice */

            /* initial dimensionless horizontal slowness */
	    p1[1] = p0+kp*dp;

            /* initial dimensionless vertical one-way slowness */
	    p1[0] = -fsqrt(1-p1[1]*p1[1]);

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
		f[3] = p1[1];
		grid1_insert(curr[kx],p1[1],4,f);
		continue;
	    }

            /* Hamiltonian vector */
	    ss = slow_bilininterp(x1,z1,slow,nx,nz,dx,dz,x0,z0);
	    hvec = hvec_init(0.,0.,x1,z1,ss*p1[1],ss*p1[0]);



            /* Loop until previous depth level or exiting from phase space grid */
	    while (pqvec->p[0] < xm or ){

     	        /* predict sigma step size to the next cell */
		ds = nc4_cellstep(hvec,slow,nx,nz,dx,dz,x0,z0,dp,dp);

                /* symplectic cell step and traveltime integration */
		nc4_sigmastep(hvec,ds,slow,nx,nz,dx,dz,x0,z0);






		if ((hvec->step) == 1) continue;


            /*
	    if (s == sx) { */
                /* exited from the side */
	    /*	if (s == sx1) {
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
	    }
            */

            /* space grid position */
	    fx = (hvec->q[1]-x0)/dx;
	    ix = snap(&fx,nx);

            /* angle grid position */
            /* convert ray parameter to angle (filt/lib/cell.c) */
	    a2 = sf_cell_p2a(hvec->p);
	    fa = (a2-a0)/da;
	    ia = snap(&fa,na);

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

	    }
	    
	    grid1_insert(curr[kx],a1,4,f);





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
    } else if (i >= n-1) {
	i=n-1;
	*f=0.;
    } else {
	*f -= i;
    }

    return i;
}
