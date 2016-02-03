/* Slope-based normal moveout for 3-D CMP geometry. */
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
#include "stretch4.h"

int main (int argc, char* argv[])
{
    map4 nmo;
    bool half;
    int it,ix,ihx,ihy,ih, nt,nx, nhx,nhy, nw;
    float dt, t0, hx,hy, h0x, h0y, t, f, dhx, dhy, eps;
    float *trace=NULL, *px=NULL, *py=NULL, *offx=NULL, *offy=NULL,  *str=NULL, *out=NULL, *vtr=NULL;

    sf_file cmp, nmod, dipx, dipy, vel;

    sf_init (argc,argv);
    cmp = sf_input("in");
    /*(Axis:Label--> 1:t, 2:hx, 3:hy, 4:cmp)*/
    dipx = sf_input("dipx");
    dipy = sf_input("dipy");
    nmod = sf_output("out");
    vel = sf_output("vel");

    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(cmp,"n2",&nhx)) sf_error("No n2= in input");
    if (!sf_histint(cmp,"n3",&nhy)) sf_error("No n2= in input");

    nx = sf_leftsize(cmp,3);
	
    if (!sf_histfloat(cmp,"d2",&dhx)) sf_error("No d2= in input");
    if (!sf_histfloat(cmp,"o2",&h0x)) sf_error("No o2= in input");
    if (!sf_histfloat(cmp,"d3",&dhy)) sf_error("No d3= in input");
    if (!sf_histfloat(cmp,"o3",&h0y)) sf_error("No o3= in input");
    
    
    if (!sf_getbool("half",&half)) half=true;
    /* if y, the second axis is half-offset instead of full offset */
    
    if (half) {
	dhx *= 2.;
	h0x *= 2.;
	dhy *= 2.;
	h0y *= 2.;
    }
    
    offx = sf_floatalloc(nhx*nhy);
    offy = sf_floatalloc(nhy*nhx);
    ih=0;
    for (ihx = 0; ihx < nhx; ihx++) {
	for (ihy = 0; ihy < nhy; ihy++) {
	    offx[ih] = h0x + ihx*dhx; 
	    offy[ih] = h0y + ihy*dhy;
	    ih++;
	}
    }

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    trace = sf_floatalloc(nt);
    px = sf_floatalloc(nt);
    py = sf_floatalloc(nt);
    str = sf_floatalloc(nt);
    out = sf_floatalloc(nt);
    vtr = sf_floatalloc(nt);

    if (!sf_getint("extend",&nw)) nw=8;
    /* trace extension */

    nmo = stretch4_init (nt, t0, dt, nt, eps);

    eps = 100.*FLT_EPSILON;

    /* BUG: figure out dh for irregular off[ih] */

    for (ix = 0; ix < nx; ix++) { /* midpoints */
/*       if(offpar){ */
/* 	sf_floatread (offx,nt,offset); */
/* 	sf_floatread (offy,nt,offset); */
/* 	nh = (sizeof(offx) / sizeof(float))*(sizeof(offy) / sizeof(float)); */
/*       } */
      sf_warning("CMP %d of %d",ix+1,nx);
      ih = 0;
      for (ihx = 0; ihx < nhx; ihx++) {
	for (ihy = 0; ihy < nhy; ihy++) {/* offsets */
	  hx=offx[ih];
	  hy=offy[ih];
	  ih++;
	  sf_floatread (trace,nt,cmp); /* data */
	  sf_floatread (px, nt, dipx);   /* slopex */
	  sf_floatread (py, nt, dipy);   /* slopey */

	  for (it=0; it < nt; it++) { /* time */
	    t = t0 + it*dt;
	    f = t - px[it]*hx*dt/dhx - py[it]*hy*dt/dhy; 

	    if (f < 0. || f > t) {
	      str[it] = t0-10.*dt;
	      vtr[it] = 0.;
	    } else {
	      str[it] = sqrtf(t*f); /* t -> tau */
	      vtr[it] = sqrtf((hx*hx+hy*hy)/(t*(t-f)+eps));
	    }
	  }
	
	
	  stretch4_define (nmo,str);

	  stretch4_apply (false,nmo,trace,out);
	  sf_floatwrite (out,nt,nmod);

	  stretch4_apply (false,nmo,vtr,vtr);
	  sf_floatwrite (vtr,nt,vel);
	}
      }
    }

    exit(0);
}
