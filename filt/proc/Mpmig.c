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
#include <float.h>

#include <rsf.h>

#include "int2.h"
#include "interp_spline.h"
#include "spline.h"

int main (int argc, char* argv[])
{
    bool half, mzo;
    int it,ix,ih, nt,nx, nh, CDPtype, ntx, nw;
    float dt, t0, h, h0, t, tm, tp, tx, tq, dh, dx, x0, x;
    float *px, *ph, **coord, *ord, *img, *img2;
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

    if (!sf_histint(cmp,"n3",&nh)) sf_error("No n2= in input");
    if (!sf_histfloat(cmp,"d3",&dh)) sf_error("No d2= in input");
    if (!sf_histfloat(cmp,"o3",&h0)) sf_error("No o2= in input");

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

    if (!sf_getint("nw",&nw)) nw=4;
    /* interpolator size (2,3,4,6,8) */

    ntx = nt*nx;

    px = sf_floatalloc(nt);
    ph = sf_floatalloc(nt);
    coord = sf_floatalloc2(2,nt);
    ord = sf_floatalloc(nt);
    img = sf_floatalloc(ntx);
    img2 = sf_floatalloc(ntx);

    if (!sf_getbool("mzo",&mzo)) mzo=false;
    /* do migration to zero offset */

    for (ih = 0; ih < nh; ih++) {
	h = h0 + (ih+0.5)*dh + (dh/CDPtype)*(ix%CDPtype); 

	for (it=0; it < ntx; it++) {
	    img[it]=0.;
	}

	for (ix = 0; ix < nx; ix++) {
	    x = x0 + ix*dx;
	    	
	    sf_floatread (ord, nt, cmp);
	    sf_floatread (px, nt, xdip);
	    sf_floatread (ph, nt, hdip);
	    
	    for (it=0; it < nt; it++) {
		t = t0 + it*dt;
		tm = t - ph[it]*h*dt/dh;
		tx = h*px[it]*dt/dh;
		tp = tm*ph[it] + px[it]*tx;
		tq = tm*tm-tx*tx;

		if (mzo) {
		    coord[it][0] = 
			sqrtf(fabsf(t*tq*tq)/
			      (fabsf(tm*tm*tm)+FLT_EPSILON));
		    coord[it][1] = x - h*tx/(tm+FLT_EPSILON);
		} else {
		    coord[it][0] = 
			sqrtf(fabsf(t*ph[it]*tq*tq)/
			      (fabsf(tm*tm*tp)+FLT_EPSILON));
		    coord[it][1] = x - t*h*px[it]/(tp+FLT_EPSILON);
		}
	    }

	    int2_init (coord, t0, x0, dt, dx, nt, nx, spline_int, nw, nt);
	    int2_lop (true,true,ntx,nt,img,ord);	    
	}
	
	/* from spline coefficients to model */
	if (nw > 2) { 
	    for (ix=0; ix < nx; ix++) {
		spline_post (nw, ix*nt, 1, nt, img, img2);
	    }
	    for (it=0; it < nt; it++) {
		spline_post (nw, it, nt, nx, img2, img);
	    }
	}
	
	sf_floatwrite (img,ntx,mig);
    }
	
    exit (0);
}

/* 	$Id: Minmo.c 729 2004-07-29 18:22:16Z fomels $	 */
