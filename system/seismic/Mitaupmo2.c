/* Inverse normal moveout in tau-p-x domain. */
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

int main (int argc, char* argv[])
{
    int it,iy,ip,ix, nt,nx,ny, np, nw, ntx;
    float dt, t0, p, p0, f0, f1, ft, dp, x0, dx, p2, v2;
    float *plane=NULL, *vel=NULL, **coord=NULL, *plane2=NULL, *trace=NULL;
    sf_file cmp=NULL, nmod=NULL, velocity=NULL;

    sf_init (argc,argv);
    cmp = sf_input("in");
    velocity = sf_input("velocity");
    nmod = sf_output("out");

    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(cmp,"n2",&np)) sf_error("No n2= in input");
    if (!sf_histfloat(cmp,"d2",&dp)) sf_error("No d2= in input");
    if (!sf_histfloat(cmp,"o2",&p0)) sf_error("No o2= in input");

    if (!sf_getint("nx",&nx)) sf_error("Need nx=");
    /* number of offsets */
    if (!sf_getfloat("dx",&dx)) sf_error("Need dx=");
    /* offset sampling */
    if (!sf_getfloat("x0",&x0)) x0=0.0;
    /* first offset */

    if (!sf_getint("nw",&nw)) nw=4;
    /* interpolator size (2,3,4,6,8) */

    ntx = nt*nx;
    ny = sf_leftsize(cmp,2);

    sf_putint(nmod,"n2",nx);
    sf_putfloat(nmod,"d2",dx);
    sf_putfloat(nmod,"o2",x0);
    sf_shiftdim(cmp, nmod, 2);

    plane = sf_floatalloc(ntx);
    plane2 = sf_floatalloc(ntx);

    vel = sf_floatalloc(nt);
    coord = sf_floatalloc2(2,nt);
    trace = sf_floatalloc(nt);

    for (iy = 0; iy < ny; iy++) {
	sf_floatread (vel,nt,velocity);	

	/* square velocity */
	for (it=0; it < nt; it++) {
	    vel[it] *= vel[it];
	}

	for (ip = 0; ip < np; ip++) {
	    p = p0 + ip*dp;
	    p2 = p*p;

	    f0 = f1 = 0.;
	    for (it=0; it < nt; it++) {
		coord[it][0] = t0+f0*dt;
		coord[it][1] = t0+f1*dt;

		v2 = vel[it];
		ft = 1.0-p2*v2;
		
		if (ft <= SF_EPS) {
		    for (; it < nt; it++) {
			coord[it][0]=t0-10.*dt;
			coord[it][1]=x0-10.*dx;
		    }
		    break;
		}
		
		ft = sqrtf(ft);
		f0 += ft;
		f1 += p*v2/ft;
	    }

	    sf_floatread (trace,nt,cmp);

	    sf_int2_init (coord, t0,x0, dt,dx, nt,nx, sf_spline_int, nw, nt);
	    sf_int2_lop (true,false,ntx,nt,plane,trace);	 
	
	    /* from spline coefficients to model */
	    if (nw > 2) { 
		for (ix=0; ix < nx; ix++) {
		    sf_spline_post (nw, ix*nt, 1, nt, plane, plane2);
		}
		for (it=0; it < nt; it++) {
		    sf_spline_post (nw, it, nt, nx, plane2, plane);
		}
	    }

	    sf_floatwrite (plane,ntx,nmod);
	}
    }

    exit (0);
}
