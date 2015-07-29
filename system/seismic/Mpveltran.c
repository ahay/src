/* Slope-based velocity transform. */
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
#include <float.h>
#include <rsf.h>

#include "warp2.h"

int main (int argc, char* argv[])
{
    bool half, inter;
    int it,ix,ih, nt,nx, nh, CDPtype, nv, ntv;
    float dt, t0, h, h0, t, f, g, dh, dy, v0, dv, eps;
    float *p, *pt, **tstr, **vstr, **slice, **vtr;
    sf_file cmp, vel, dip, dipt;

    sf_init (argc,argv);
    cmp = sf_input("in");
    dip = sf_input("dip");
    vel = sf_output("out");

    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(cmp,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histfloat(cmp,"d2",&dh)) sf_error("No d2= in input");
    if (!sf_histfloat(cmp,"o2",&h0)) sf_error("No o2= in input");

    if (!sf_getbool("half",&half)) half=true;
    /* if y, the second axis is half-offset instead of full offset */

    if (half) {
	dh *= 2.;
	h0 *= 2.;
    }

    CDPtype=1;
    if (sf_histfloat(cmp,"d3",&dy)) {
	CDPtype=0.5+0.5*dh/dy;
	if (1 != CDPtype) sf_histint(cmp,"CDPtype",&CDPtype);
    }
    if (CDPtype < 1) CDPtype=1;
    sf_warning("CDPtype=%d",CDPtype);

    if (!sf_getint("nv",&nv)) sf_error("Need nv=");
    /* number of velocities */
    if (!sf_getfloat("v0",&v0)) sf_error("Need v0=");
    /* velocity origin */
    if (!sf_getfloat("dv",&dv)) sf_error("Need dv=");
    /* velocity sampling */

    sf_putint(vel,"n2",nv);
    sf_putfloat(vel,"o2",v0);
    sf_putfloat(vel,"d2",dv);

    nx = sf_leftsize(cmp,2);

    ntv = nt*nv;

    p = sf_floatalloc(nt);
    tstr = sf_floatalloc2(nt,nh);
    vstr = sf_floatalloc2(nt,nh);
    slice = sf_floatalloc2(nt,nh);
    vtr = sf_floatalloc2(nt,nv);

    if (!sf_getbool("interval",&inter)) inter=false;
    /* if y, compute interval velocity */

    if (inter) {
	dipt = sf_input("dipt");
	pt = sf_floatalloc(nt);
    } else {
	dipt = NULL;
	pt = NULL;
    }

    if (!sf_getfloat("eps",&eps)) eps=0.1;
    /* stretch regularization */

    warp2_init (nt, t0, dt,
		nv, v0, dv,
		nt, nh, eps);

    for (ix = 0; ix < nx; ix++) {
	sf_floatread (slice[0], nt*nh, cmp);

	for (ih = 0; ih < nh; ih++) {
	    h = h0 + (ih+0.5)*dh + (dh/CDPtype)*(ix%CDPtype);

	    sf_floatread (p, nt, dip);
	    if (inter) sf_floatread (pt, nt, dipt);

	    for (it=0; it < nt; it++) {
		t = t0 + it*dt;
		f = t - p[it]*h*dt/dh;

		if (f < 0. || f > t) {
		    tstr[ih][it]=t0-10.*dt;
		    vstr[ih][it]=0.;
		} else {
		    tstr[ih][it] = sqrtf(t*f);
		    if (inter) {
			g = dt*p[it]+t*pt[it];
			vstr[ih][it] = 
			    sqrtf(fabsf(
				      h*dh*(p[it]*h*g*dt - 2*pt[it]*t*t*dh)/
				      (p[it]*p[it]*t*(2*t*dh-h*g)+FLT_EPSILON)
				      ))/dt;
		    } else {
			vstr[ih][it] = h/sqrtf(t*(t-f)+FLT_EPSILON);
		    }
		}
	    }
	} /* ih */

	warp2(slice,tstr,vstr,vtr);
	sf_floatwrite (vtr[0],ntv,vel);
    }

    exit (0);
}
