/* 3-D Inverse normal moveout.*/
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

#include "stretch4.h"

int main (int argc, char* argv[])
{
    map4 nmo; /* using cubic spline interpolation */
    bool half, slow, ellipse;
    int it,ix,iy, ih, nt, nx, ny, nh, nhx, nhy, nw, CDPtype;
    float dt, t0, hx, hy, h0x, h0y, h, f, dhx, dhy, eps, dc, vx, vy, vxy;
    float *trace, *vel, *off,*offx,*offy, *str, *out;
    sf_file cmp, nmod, velocity, offset;

    sf_init (argc,argv);
    cmp = sf_input("in");
    nmod = sf_output("out");

    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(cmp,"n2",&nhx)) sf_error("No n2= in input");
    if (!sf_histint(cmp,"n3",&nhy)) sf_error("No n3= in input");

    offx = sf_floatalloc(nhx);
    offy = sf_floatalloc(nhy);

    if (NULL != sf_getstring("velocity")) {
	velocity = sf_input("velocity");
    }else{
      sf_warning("No velocity file specified. Using velocity ellipse parameters.");
      ellipse=true;
      if (NULL != sf_getstring("vx")) sf_error("No vx= in input. Must specify velocity file, or vx= vy= vxy= parameters.");
      if (NULL != sf_getstring("vy")) sf_error("No vy= in input. Must specify velocity file, or vx= vy= vxy= parameters.");
      if (NULL != sf_getstring("vxy")) sf_error("No vxy= in input. Must specify velocity file, or vx= vy= vxy= parameters.");
    }

    CDPtype=1;
    if (NULL != sf_getstring("offset")) {
	offset = sf_input("offset");
	sf_floatread (off,nh,offset);
	sf_fileclose(offset);
    } else {
	if (!sf_histfloat(cmp,"d2",&dhx)) sf_error("No d2= in input");
	if (!sf_histfloat(cmp,"o2",&h0x)) sf_error("No o2= in input");
	if (!sf_histfloat(cmp,"d3",&dhy)) sf_error("No d3= in input");
	if (!sf_histfloat(cmp,"o3",&h0y)) sf_error("No o3= in input");
    

	if (!sf_getbool("half",&half)) half=true;
	/* if y, the second and third axes are half-offset instead of full offset */
	
	if (half) {
	    dhx *= 2.;
	    dhy *= 2.;
	    h0x *= 2.;
	    h0y *= 2.;
	}
	
	if (sf_histfloat(cmp,"d4",&dc)) {
	    CDPtype=0.5+0.5*dhx/dc;
	    if (1 != CDPtype) sf_histint(cmp,"CDPtype",&CDPtype);
	} 	   
	if (1 > CDPtype) CDPtype=1;
	sf_warning("CDPtype=%d",CDPtype);

	for (ih = 0; ih < nh; ih++) {
	  offx[ih] = h0x + ih*dhx;
	}
	for (ih = 0; ih < nh; ih++) {
	  offy[ih] = h0y + ih*dhy;
	}

    }

    if (!sf_getbool("slowness",&slow)) slow=false;
    /* if y, use slowness instead of velocity */

    if (!sf_getbool("ellipse",&ellipse)) ellipse=false;
    /* if y, use ellipse instead of velocity */
    

    nx = sf_leftsize(cmp,2);
    ny = sf_leftsize(cmp,3);

    if (!sf_getfloat ("h0x",&h0x)) h0x=0.;
    if (!sf_getfloat ("h0y",&h0y)) h0y=0.;
    /* reference offset */
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    trace = sf_floatalloc(nt);
    vel = sf_floatalloc(nt);
    str = sf_floatalloc(nt);
    out = sf_floatalloc(nt);

    if (!sf_getint("extend",&nw)) nw=8;
    /* trace extension */

    nmo = stretch4_init (nt, t0, dt, nt, eps);

    if (ellipse == true) {

      for (ix = 0; ix < nx; ix++) {
	for (iy = 0; iy < ny; iy++) {
	  for (ih = 0; ih < nhx; ih++) {
	    /*Should probably window data into 2-D gathers and loop.*/
	    /*As long as x and y coords are there, should be 3-D ok.*/
	    sf_floatread (trace,nt,cmp);
	    
	    hx = offx[ih] + (dhx/CDPtype)*(ix%CDPtype) - h0x;
	    hy = offy[ih] + (dhy/CDPtype)*(iy%CDPtype) - h0y;
	    /*h = h*h - h0*h0;*/
	    
	    for (it=0; it < nt; it++) {
	      f = t0 + it*dt;
	      f = f*f + hx*hx/(vx*vx)+hy*hy/(vy*vy)+hx*hy/(vxy*vxy);
		
	      if (f < 0.) {
		str[it]=t0-10.*dt;
	      } else {
		str[it] = sqrtf(f);
	      }
	    }

	    stretch4_define (nmo,str);
	    stretch4_apply (nmo,trace,out);
	    
	    sf_floatwrite (out,nt,nmod);
	  }
	}
      }

    } else {

      for (ix = 0; ix < nx; ix++) {
	for (iy = 0; iy < ny; iy++) {
	  sf_floatread (vel,nt,velocity);	

	  for (ih = 0; ih < nh; ih++) {
	    sf_floatread (trace,nt,cmp);
	    
	    hx = offx[ih] + (dhx/CDPtype)*(ix%CDPtype); 
	    hy = offx[ih] + (dhy/CDPtype)*(iy%CDPtype);
	    h = (hx-h0x)*(hx-h0x)+(hy-h0y)*(hy-h0y);
	    
	    for (it=0; it < nt; it++) {
		f = t0 + it*dt;
		if (slow) {
		    f = f*f + h*vel[it]*vel[it];
		} else {
		    f = f*f + h/(vel[it]*vel[it]);
		}
		if (f < 0.) {
		    str[it]=t0-10.*dt;
		} else {
		    str[it] = sqrtf(f);
		}
	    }

	    stretch4_define (nmo,str);
	    stretch4_apply (nmo,trace,out);
	    
	    sf_floatwrite (out,nt,nmod);
	  }
	}
      }
    }
	
    exit (0);
}

/* 	$Id: Minmo.c 1806 2006-04-20 13:21:47Z fomels $	 */
