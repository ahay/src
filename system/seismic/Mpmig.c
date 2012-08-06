/* Slope-based prestack time migration. */
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

#include "warp2.h"

int main (int argc, char* argv[])
{
    bool half, mzo;
    int it,ix,ih, nt,nx, nh, CDPtype, ntx;
    float dt, t0, h, h0, t, tm, tp, tx, tq, dh, dx, x0, x;
    float *px, *ph, **xstr, **tstr, **slice, **img, eps;
    sf_file cmp, mig, xdip, hdip;

    sf_init (argc,argv);
    cmp = sf_input("in");
    xdip = sf_input("xdip");
    hdip = sf_input("hdip");
    mig = sf_output("out");

    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(cmp,"n3",&nh)) sf_error("No n3= in input");
    if (!sf_histfloat(cmp,"d3",&dh)) sf_error("No d3= in input");
    if (!sf_histfloat(cmp,"o3",&h0)) sf_error("No o3= in input");

    if (!sf_histint(cmp,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(cmp,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(cmp,"o2",&x0)) sf_error("No o2= in input");

    if (!sf_getbool("half",&half)) half=true;
    /* if y, the second axis is half-offset instead of full offset */
	
    if (!half) {
	dh *= 0.5;
	h0 *= 0.5;
    }

    CDPtype=0.5+dh/dx;
    if (1 != CDPtype) sf_histint(cmp,"CDPtype",&CDPtype);
    sf_warning("CDPtype=%d",CDPtype);

    ntx = nt*nx;

    px = sf_floatalloc(nt);
    ph = sf_floatalloc(nt);
    tstr = sf_floatalloc2(nt,nx);
    xstr = sf_floatalloc2(nt,nx);
    slice = sf_floatalloc2(nt,nx);
    img = sf_floatalloc2(nt,nx);

    if (!sf_getbool("mzo",&mzo)) mzo=false;
    /* do migration to zero offset */

    if (!sf_getfloat("eps",&eps)) eps=1.0;
    /* stretch regularization */    

    warp2_init (nt, t0, dt,
		nx, x0, dx,
		nt, nx, eps);

    for (ih = 0; ih < nh; ih++) {
	sf_floatread (slice[0], ntx, cmp);

	for (ix = 0; ix < nx; ix++) {
	    x = x0 + ix*dx;
	    h = h0 + (ih+0.5)*dh + (dh/CDPtype)*(ix%CDPtype); 

	    sf_floatread (px, nt, xdip);
	    sf_floatread (ph, nt, hdip);
	    
	    for (it=0; it < nt; it++) {
		t = t0 + it*dt;
		tm = t - ph[it]*h*dt/dh;
		tx = h*px[it]*dt/dx;
		tp = tm*ph[it] + px[it]*tx*dh/dx;
		tq = tm*tm-tx*tx;

		if (mzo) {
		    tstr[ix][it] = 
			sqrtf(fabsf(t*tq*tq)/
			      (fabsf(tm*tm*tm)+SF_EPS));
		    xstr[ix][it] = x - h*tx/(tm+SF_EPS);
		} else {
		    tstr[ix][it] = 
			sqrtf(fabsf(t*ph[it]*tq*tq)/
			      (fabsf(tm*tm*tp)+SF_EPS));
		    xstr[ix][it] = x - t*h*px[it]/(tp+SF_EPS);
		}
	    }
	}
	
	warp2(slice,tstr,xstr,img);
	sf_floatwrite (img[0],ntx,mig);
    }


    exit (0);
}

