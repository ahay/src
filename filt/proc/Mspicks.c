/* Generate stereotomography picks from time migration. */
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
    int nt, nx, nh, it, ix, ih;
    const int na=5; /* number of data attributes in a pick */
    float dt, t0, dh, h0, dx, x0, h, h2, r, v, p, t, t2, x, y, q, py, ph;
    float *imgt, *velt, *dipt, ***pikt;
    sf_file img, vel, dip, pik;

    sf_init (argc,argv);
    img = sf_input("in");
    vel = sf_input("vel");
    dip = sf_input("dip");
    pik = sf_output("out");

    if (SF_FLOAT != sf_gettype(img)) sf_error("Need float input");
    if (!sf_histint(img,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(img,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(img,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histint(img,"n2",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(img,"d2",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(img,"o2",&x0)) sf_error("No o1= in input");

    if (!sf_getint("nh",&nh)) sf_error("Need nh=");
    /* number of offsets */
    if (!sf_getfloat("dh",&dh)) sf_error("Need dh=");
    /* offset sampling */
    if (!sf_getfloat("h0",&h0)) h0=0.;
    /* offset origin */
	
    sf_putint(pik,"n1",na);
    sf_putint(pik,"n2",nh);
    sf_putfloat(pik,"d2",dh);
    sf_putfloat(pik,"o2",h0);
    sf_putint(pik,"n3",nt);
    sf_putfloat(pik,"d3",dt);
    sf_putfloat(pik,"o3",t0);
    sf_putint(pik,"n4",nx);
    sf_putfloat(pik,"d4",dx);
    sf_putfloat(pik,"o4",x0);

    imgt = sf_floatalloc(nt);
    velt = sf_floatalloc(nt);
    dipt = sf_floatalloc(nt);
    pikt = sf_floatalloc3(na,nh,nt);

    for (ix=0; ix < nx; ix++) {
	x = x0 + ix*dx;

	sf_floatread(imgt,nt,img);
	sf_floatread(velt,nt,vel);
	sf_floatread(dipt,nt,dip);
	for (it=0; it < nt; it++) {
	    t = t0 + it*dt;

	    v = 0.5*velt[it];
	    v *= v;

	    p = dipt[it]*dt/dx;
	    r = 1.+p*p*v;

	    for (ih=0; ih < nh; ih++) {
		h = h0 + ih*dh;
		h2 = h*h/v;

		/* time  squared */
		t2 = h2+0.5*t*(t*r+hypotf(t*r,2*h*p)); 

		q = t2*r - h2;

		if (q > 0.01*dt*dt) {
		    q = sqrtf(q);		    
		} else {
		    q = 0.1*dt;
		}

		/* midpoint */
		y = x + p*v*t2/q;

		t2 = sqrtf(t2);

		/* time */
		pikt[it][ih][0] = t2;

		/* source */
		pikt[it][ih][1] = y - h;

		/* receiver */
		pikt[it][ih][2] = y + h;

		/* offset slope */
		ph = dt*h/(dt*t2*r*v+0.0001*dx*dx);

		/* midpoint slope */
		py = p*q/(t2*r+0.01*dt);

		/* source slope */
		pikt[it][ih][3] = 0.5*(py-ph);
		
		/* receiver slope */
		pikt[it][ih][4] = 0.5*(py+ph);
	    }
	}
	sf_floatwrite(pikt[0][0],na*nh*nt,pik);
    }

    exit (0);
}

/* 	$Id: Mspicks.c 729 2004-07-29 18:22:16Z fomels $	 */
