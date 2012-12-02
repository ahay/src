/* Oriented velocity continuation. 

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
    int it, nt, ip, np, ix, nx, ntpx, iv, nv;
    float eps, t, t0, dt, p, p0, dp, x, x0, dx, v0, v, dv, dv2, sq;
    float ***slice, ***slice0, ***tstr, ***pstr, ***xstr;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    
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

    if (!sf_getfloat("eps",&eps)) eps=0.1;
    /* stretch regularization */

    if (!sf_getint("nv",&nv)) nv=1; /* number of velocity steps */

    if (!sf_getbool("lagrange",&lagrange)) lagrange=false;
    /* Use Lagrangian method */

    if (!sf_getfloat("v0",&v0)) v0=0.; /* starting velocity */
    if (!sf_getfloat("vmax",&v)) sf_error("Need vmax="); /* end velocity */

    dv = (v-v0)/nv;

    sf_putint(out,"n4",nv);
    sf_putfloat(out,"o4",v0+dv);
    sf_putfloat(out,"d4",dv);
    
    tstr   = sf_floatalloc3(nt,np,nx);
    pstr   = sf_floatalloc3(nt,np,nx);
    xstr   = sf_floatalloc3(nt,np,nx); 

    warp3_init(nt, t0, dt,
	       np, p0, dp,
               nx, x0, dx,
	       nt, np, nx, eps); 

    if (lagrange) {
	slice0 = sf_floatalloc3(nt,np,nx);
	slice  = sf_floatalloc3(nt,np,nx);
    } else {
	slice  = sf_floatalloc3(nt,np,nx);	
	slice0 = slice;
    } 
		
    sf_floatread(slice0[0][0],ntpx,in);

    for (iv=1; iv <= nv; iv++) {
	sf_warning("step %d of %d;",iv,nv);

	v = v0+iv*dv;

	if (lagrange) {
	    dv2 = 0.25*(v*v-v0*v0);
	} else {
	    dv2 = 0.25*(2*v-dv)*dv;
	}
	    
	for (ix=0; ix < nx; ix++) {
	    x = x0+ix*dx;
	    for (ip=0; ip < np; ip++) {
		p = p0+ip*dp;
		sq = sqrtf(1-p*p*dv2);		
		for (it=0; it < nt; it++) {
		    t = t0+it*dt;
		    
		    tstr[ix][ip][it] = t*sq;
		    pstr[ix][ip][it] = p/sq;
		    xstr[ix][ip][it] = x - t*p*dv2;
		}
	    }
	}

	warp3(slice0,tstr,pstr,xstr,slice);	
	sf_floatwrite (slice[0][0],ntpx,out);
    }
    sf_warning(".");

    exit(0);
}
