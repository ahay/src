/* Velocity analysis using T-X-P domain. */
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

#include "fint1.h"

static float v,p,x;

static float vscan(float t, int it)
{
    return t+p*x*x/(x+t*p*v*v);
}

int main(int argc, char* argv[])
{
    fint1 nmo;
    int it, nt, ix, nx, ip, np, iv, nv, nw, mute;
    float t0, x0, p0, v0, dt, dx, dp, dv, str;
    float *trace, **stack;
    sf_file cmp, scan;

    sf_init (argc,argv);
    cmp = sf_input("in");
    scan = sf_output("out");

    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(cmp,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(cmp,"n3",&np)) sf_error("No n3= in input");
    
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o2",&x0)) sf_error("No o2= in input");
    if (!sf_histfloat(cmp,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(cmp,"o3",&p0)) sf_error("No o3= in input");
    if (!sf_histfloat(cmp,"d3",&dp)) sf_error("No d3= in input");

    if (!sf_getfloat("v0",&v0)) sf_error("Need v0=");
    /* first scanned velocity */
    if (!sf_getfloat("dv",&dv)) sf_error("Need dv=");
    /* step in velocity */
    if (!sf_getint("nv",&nv)) sf_error("Need nv=");
    /* number of scanned velocities */

    sf_putfloat(scan,"o2",v0);
    sf_putfloat(scan,"d2",dv);
    sf_putint(scan,"n2",nv);

    sf_putint(scan,"n3",1);

   if (!sf_getint("extend",&nw)) nw=4;
    /* trace extension */

    if (!sf_getint("mute",&mute)) mute=12;
    /* mute zone */

    if (!sf_getfloat("str",&str)) str=0.5;
    /* maximum stretch allowed */
    
    trace = sf_floatalloc(nt);
    nmo = fint1_init(nw,nt,mute);

    stack =  sf_floatalloc2(nt,nv);
    
    for (it=0; it < nt*nv; it++) {
	stack[0][it] = 0.;
    }

    for (ip=0; ip < np; ip++) {
	p = p0+ip*dp;

	for (ix=0; ix < nx; ix++) {
	    x = x0+ix*dx;

	    sf_floatread(trace,nt,cmp); 

	    /* normalize */
	    for (it=0; it < nt; it++) {
		trace[it] /= nt*np*nx;
	    }

	    for (iv=0; iv < nv; iv++) {
		v = v0 + iv * dv;

		stretch(nmo,vscan,nt,dt,t0,nt,dt,t0,trace,str);

		for (it=0; it < nt; it++) {
		    stack[iv][it] += trace[it];
		}
	    } /* v */
	} /* x */
    } /* p */

    sf_floatwrite (stack[0],nt*nv,scan);
    exit(0);
}
