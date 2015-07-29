/* Velocity transform.

Inverse of sfvscan.
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

#include <math.h>
#include <rsf.h>
#include "fint1.h"

static float h;

static float hyperb(float t, int it) 
{
    if (t > h) {
	return sqrtf(t*t-h*h);
    } else {
	return 0.;
    }
}

int main(int argc, char* argv[])
{
    fint1 nmo;
    bool half, slow;
    int it,ih,ix,iv, nt,nh,nx,nv, nw, CDPtype;
    float dt, dh, t0, h0, v0, dv, v, dy;
    float *trace=NULL, **stack=NULL;
    sf_file cmp=NULL, scan=NULL;

    sf_init (argc,argv);
    scan = sf_input("in");
    cmp = sf_output("out");

    if (!sf_histint(scan,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(scan,"n2",&nv)) sf_error("No n2= in input");
    nx = sf_leftsize(scan,2);

    if (!sf_histfloat(scan,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histfloat(scan,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(scan,"o2",&v0)) sf_error("No o2= in input");
    if (!sf_histfloat(scan,"d2",&dv)) sf_error("No d2= in input");

    if (!sf_getfloat("h0",&h0) && !sf_histfloat(scan,"h0",&h0)) 
	sf_error("Need h0=");
    if (!sf_getfloat("dh",&dh) && !sf_histfloat(scan,"dh",&dh)) 
	sf_error("Need dh=");
    if (!sf_getint("nh",&nh) && !sf_histint(scan,"nh",&nh)) 
	sf_error("Need nh=");

    sf_putfloat(cmp,"o2",h0);
    sf_putfloat(cmp,"d2",dh);
    sf_putint(cmp,"n2",nh);

    sf_putfloat(cmp,"v0",v0);
    sf_putfloat(cmp,"dv",dv);
    sf_putint(cmp,"nv",nv);

    if (!sf_getbool("half",&half)) half=true;
    /* if y, the second axis is half-offset instead of full offset */

    sf_putstring(cmp,"label2",half? "half-offset": "offset");

    if (half) {
	dh *= 2.;
	h0 *= 2.;
    }

    CDPtype=1;
    if (sf_histfloat(scan,"d3",&dy)) {
	CDPtype=0.5+0.5*dh/dy;
	if (0 == CDPtype) CDPtype=1;
	if (1 != CDPtype) sf_histint(scan,"CDPtype",&CDPtype);
    }
    sf_warning("CDPtype=%d",CDPtype);

    trace = sf_floatalloc(nt);
    stack =  sf_floatalloc2(nt,nh);

    if (!sf_getint("extend",&nw)) nw=4;
    /* trace extension */

    trace = sf_floatalloc(nt);
    nmo = fint1_init(nw,nt,0);

    if (!sf_getbool("slowness",&slow)) slow=false;
    /* if y, use slowness instead of velocity */

    for (ix=0; ix < nx; ix++) {
	sf_warning("scan %d of %d;",ix+1,nx);

	for (it=0; it < nt*nh; it++) {
	    stack[0][it] = 0.;
	}

	for (iv=0; iv < nv; iv++) {
	    v = v0 + iv * dv;
	    sf_floatread(trace,nt,scan); 

	    for (it=0; it < nt; it++) {
		trace[it] /= nt*nh;
	    }
	    fint1_set(nmo,trace);

	    for (ih=0; ih < nh; ih++) {
		h = h0 + ih * dh + (dh/CDPtype)*(ix%CDPtype);
		h = slow? h*v: h/v;

		stretch(nmo,hyperb,nt,dt,t0,nt,dt,t0,trace,0.);

		for (it=0; it < nt; it++) {
		    stack[ih][it] += trace[it];
		}
	    } /* h */
	} /* v */
	
	sf_floatwrite (stack[0],nt*nh,cmp);
    } /* x */
    sf_warning(".");

    exit(0);
}
