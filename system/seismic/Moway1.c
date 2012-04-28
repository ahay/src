/* Oriented one-way wave equation. 

Axis order: t, p, x
*/
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
#include "warp3.h"

int main(int argc, char* argv[])
{
    bool lagrange;
    int it, nt, ip, np, ix, nx, ntpx, iz, nz, ny;
    float eps, t, t0, dt, p, p0, dp, x, x0, dx, sx, st, sp, z0, dz, sq, v, g;
    float ***slice, ***slice0, ***tstr, ***pstr, ***xstr, *vv, *vg;
    sf_file in, out, vel, vgrad;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    vel = sf_input("velocity");
    vgrad = sf_input("vgrad");
    
    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&np)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&nx)) sf_error("No n3= in input");

    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&dp)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"d3",&dx)) sf_error("No d3= in input");

    if (!sf_histfloat(in,"o1",&t0)) t0=0.;
    if (!sf_histfloat(in,"o2",&p0)) p0=0.;
    if (!sf_histfloat(in,"o3",&x0)) x0=0.;

    ntpx = nt*np*nx;

    if (!sf_histint(vel,"n1",&ny) || ny != nx) sf_error("Need n1=%d in velocity",nx);
    if (!sf_histint(vel,"n2",&nz)) sf_error("No n2= in velocity");
    if (!sf_histfloat(vel,"d2",&dz))  sf_error("No d2= in velocity");
    if (!sf_histfloat(vel,"o2",&z0)) z0=0.;

    if (!sf_getfloat("eps",&eps)) eps=0.1;
    /* stretch regularization */

    if (!sf_getbool("lagrange",&lagrange)) lagrange=false;
    /* Use Lagrangian method */    

    sf_putint(out,"n4",nz);
    sf_putfloat(out,"o4",z0);
    sf_putfloat(out,"d4",dz);
    
    tstr = sf_floatalloc3(nt,np,nx);
    pstr = sf_floatalloc3(nt,np,nx);
    xstr = sf_floatalloc3(nt,np,nx); 

    vv = sf_floatalloc(nx);
    vg = sf_floatalloc(nx);

    warp3_init(nt, t0, dt,
	       np, p0, dp,
               nx, x0, dx,
	       nt, np, nx, eps); 

    slice  = sf_floatalloc3(nt,np,nx);

    if (lagrange) {
	slice0 = sf_floatalloc3(nt,np,nx);
	sx = sp = st = 0.;
    } else {
	slice0 = slice;
    } 
    
    sf_floatread(slice0[0][0],ntpx,in);

    for (iz=0; iz < nz; iz++) {
	sf_warning("depth %d of %d;",iz+1,nz);

	sf_floatwrite (slice[0][0],ntpx,out);
	if (iz==nz-1) break;

	sf_floatread(vv,nx,vel);
	sf_floatread(vg,nx,vgrad);
 	
	for (ix=0; ix < nx; ix++) {
	    x = x0+ix*dx;
	    v = vv[ix];
	    g = vg[ix];

	    for (ip=0; ip < np; ip++) {
		p = p0+ip*dp;

		sq = 1./sqrt(1.-v*v*p*p);

		if (lagrange) {
		    sx -= sq*p*v;
		    sp += sq*g/(v*v);
		    st -= sq/v;
		} else {
		    sx = -sq*p*v;
		    sp = sq*g/(v*v);
		    st = -sq/v;
		}
		
		for (it=0; it < nt; it++) {
		    t = t0+it*dt;

		    tstr[ix][ip][it] = t+st*dz;
		    pstr[ix][ip][it] = p+sp*dz;
		    xstr[ix][ip][it] = x+sx*dz;
		}
	    }
	}

	warp3(slice0,tstr,pstr,xstr,slice);
    }
    sf_warning(".");
    
    exit(0);
}
