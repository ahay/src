/* Slope-based prestack (t,xs,h) time migration . */
/*
  Copyright (C) 2004 University of Texas at Austin
  Copyright (C) 2010 Politecnico di milano	
  
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
    bool half, mzo;
    int it,ix,ih,is, nt, ns, nh, nx, CDPtype, ntx, nw;
    float dt, t0, h, h0, t, tm, tp, tx, tq, dh, ds, s0, x, dx, x0, xe, px,s, ph0;
    float *ps, *ph, **coord, *ord, *img, *img2;
    sf_file csg, mig, sdip, hdip;

    sf_init (argc,argv);
    csg = sf_input("in");
    sdip = sf_input("sdip");
    hdip = sf_input("hdip");
    mig = sf_output("out");

    if (SF_FLOAT != sf_gettype(csg)) sf_error("Need float input");
    if (!sf_histint(csg,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(csg,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(csg,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(csg,"n3",&nh)) sf_error("No n3= in input");
    if (!sf_histfloat(csg,"d3",&dh)) sf_error("No d3= in input");
    if (!sf_histfloat(csg,"o3",&h0)) sf_error("No o3= in input");

    if (!sf_histint(csg,"n2",&ns)) sf_error("No n2= in input");
    if (!sf_histfloat(csg,"d2",&ds)) sf_error("No d2= in input");
    if (!sf_histfloat(csg,"o2",&s0)) sf_error("No o2= in input");

    if (!sf_getbool("half",&half)) half=true;
    /* if y, the second axis is half-offset instead of full offset */
	
    if (!half) {
	dh *= 0.5;
	h0 *= 0.5;
   
    }
    x0 = s0+h0;
    xe = x0 + (ns-1)*ds + (nh-1)*dh;
    dx = dh;
    nx = (int) (round ( (xe-x0)/dx )) + 1;

    sf_warning("x0=%f xe=%f  dx=%f nCMP=%d",x0,xe,dx,nx);
    
    sf_putint(mig,"n2",nx);
    sf_putfloat(mig,"o2",x0);
    sf_putfloat(mig,"d2",dx);

    CDPtype=0.5+dh/dx;
    if (1 != CDPtype) sf_histint(csg,"CDPtype",&CDPtype);
    sf_warning("CDPtype=%d",CDPtype);

    if (!sf_getint("nw",&nw)) nw=4;
    /* interpolator size (2,3,4,6,8) */

    ntx = nt*nx;

    ps = sf_floatalloc(nt);
    ph = sf_floatalloc(nt);
    coord = sf_floatalloc2(2,nt);
    ord = sf_floatalloc(nt);
    img = sf_floatalloc(ntx);
    img2 = sf_floatalloc(ntx);

    if (!sf_getbool("mzo",&mzo)) mzo=false;
    /* do migration to zero offset */


    for (ih = 0; ih < nh; ih++) {
	sf_warning("offset=%d/%d",ih,nh);
	for (it=0; it < ntx; it++) {
	    img[it]=0.;
	}

	for (is = 0; is < ns; is++) {
	    s = s0 + is*ds;
	    //h = h0 + (ih+0.5)*dh + (dh/CDPtype)*(ix%CDPtype); 
	    h = h0+(ih+0.5)*dh;

	    x = s + h ;
	    
	    sf_floatread (ord, nt, csg);
	    sf_floatread (ps, nt, sdip);
	    sf_floatread (ph, nt, hdip);
	    
	    
	    for (it=0; it < nt; it++) {
		t = t0 + it*dt;
		
	
		px = ps[it]; 
	    ph0 =  ph[it]-ps[it];


		tm = t - ph0*h*dt/dh;
		tx = h*px*dt/dh;
		tp = tm*ph0 + px*tx;
		tq = tm*tm-tx*tx;
		
		if (mzo) {
		    coord[it][0] = 
			sqrtf(fabsf(t*tq*tq)/
			      (fabsf(tm*tm*tm)+SF_EPS));
		    coord[it][1] = x - h*tx/(tm+SF_EPS);
		} else {
		    coord[it][0] = 
			sqrtf(fabsf(t*ph0*tq*tq)/
			      (fabsf(tm*tm*tp)+SF_EPS));
		    coord[it][1] = x - t*h*px/(tp+SF_EPS);
		}
	    }
	    
	    sf_int2_init (coord, t0,x0, dt,dx, nt,nx, sf_spline_int, nw, nt);
	    sf_int2_lop (true,true,ntx,nt,img,ord);	    
	}
	
	/* from spline coefficients to model */
	if (nw > 2) { 
	    for (ix=0; ix < nx; ix++) {
		sf_spline_post (nw, ix*nt, 1, nt, img, img2);
	    }
	    for (it=0; it < nt; it++) {
		sf_spline_post (nw, it, nt, nx, img2, img);
	    }
	}
	
	sf_floatwrite (img,ntx,mig);
    }
    
    
    exit (0);
}

/* 	$Id: Mpmig.c 4796 2009-09-29 19:39:07Z ivlad $	 */
