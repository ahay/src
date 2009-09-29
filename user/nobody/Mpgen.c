/* Generate stereopicks from time-migration velocities and slopes. */
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

int main (int argc, char* argv[])
{
    bool half;
    int it,ix,ih, nt,nx, nh, CDPtype;
    float z, dt, t0, h, h0, dh, x, dx, x0, r, q, vv, pp, t2;
    float *t, *y, *py, *ph, *p, *v;
    sf_file vel, dip, tim, cmp, ydip, hdip;

    sf_init (argc,argv);
    vel = sf_input("in");
    dip = sf_input("dip");
 
    tim = sf_output("out");
    cmp = sf_output("cmp");
    ydip = sf_output("ydip");
    hdip = sf_output("hdip");

    if (SF_FLOAT != sf_gettype(vel)) sf_error("Need float input");
    if (!sf_histint(vel,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(vel,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(vel,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(vel,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(vel,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(vel,"o2",&x0)) sf_error("No o2= in input");

    if (!sf_getint("nh",&nh)) nh=1;
    if (!sf_getfloat("dh",&dh)) dh=dx;
    if (!sf_getfloat("h0",&h0)) h0=0.;

    sf_putint(tim,"n3",nh);
    sf_putfloat(tim,"d3",dh);
    sf_putfloat(tim,"o3",h0);

    sf_putint(cmp,"n3",nh);
    sf_putfloat(cmp,"d3",dh);
    sf_putfloat(cmp,"o3",h0);

    sf_putint(ydip,"n3",nh);
    sf_putfloat(ydip,"d3",dh);
    sf_putfloat(ydip,"o3",h0);

    sf_putint(hdip,"n3",nh);
    sf_putfloat(hdip,"d3",dh);
    sf_putfloat(hdip,"o3",h0);

    if (!sf_getbool("half",&half)) half=true;
    /* if y, the second axis is half-offset instead of full offset */
	
    if (!half) {
	dh *= 0.5;
	h0 *= 0.5;
    }

    CDPtype=0.5+dh/dx;
    if (1 != CDPtype) sf_histint(cmp,"CDPtype",&CDPtype);
    sf_warning("CDPtype=%d",CDPtype);

    t = sf_floatalloc(nt);
    y = sf_floatalloc(nt);
    py = sf_floatalloc(nt);
    ph = sf_floatalloc(nt);
    p = sf_floatalloc(nt);
    v = sf_floatalloc(nt);


    for (ih = 0; ih < nh; ih++) {
	for (ix = 0; ix < nx; ix++) {
	    x = x0 + ix*dx;
	    h = h0 + (ih+0.5)*dh + (dh/CDPtype)*(ix%CDPtype); 
	    
	    sf_floatread (v, nt, vel);
	    sf_floatread (p, nt, dip);
	    
	    for (it=0; it < nt; it++) {
		z = t0 + it*dt;
		
		vv = 0.25*v[it]*v[it];
		pp = p[it];

		r = 1.+pp*pp*vv;
		t2 = h*h/vv + 0.5*z*(z*r+hypotf(z*r,2.*h*pp));
		q = sqrtf(t2*r-h*h/vv);

		t[it] = sqrtf(t2);
		y[it] = x + pp*vv*t2/q;
		py[it] = pp*q/(t[it]*r+FLT_EPSILON);
		ph[it] = h/(t[it]*r*vv+FLT_EPSILON);
	    }

	    sf_floatwrite(t, nt, tim);
	    sf_floatwrite(y, nt, cmp);
	    sf_floatwrite(py, nt, ydip);
	    sf_floatwrite(ph, nt, hdip);
	}
    }

    exit (0);
}

/* 	$Id: Minmo.c 729 2004-07-29 18:22:16Z fomels $	 */
