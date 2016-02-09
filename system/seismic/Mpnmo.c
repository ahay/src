/* Slope-based normal moveout. */
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
    sf_map4 nmo;
    bool half;
    int it,ix,ih, nt,nx, nh, CDPtype;
    float dt, t0, h, h0, t, f, g, dh, eps, dy;
    float *trace, *p, *q, *off, *str, *out, *vtr, *etr;
    sf_file cmp, nmod, dip, offset, crv, vel, eta;

    sf_init (argc,argv);
    cmp = sf_input("in");
    dip = sf_input("dip");
    nmod = sf_output("out");
    vel = sf_output("vel");

    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(cmp,"n2",&nh)) sf_error("No n2= in input");

    off = sf_floatalloc(nh);

    CDPtype=1;
    if (NULL != sf_getstring("offset")) {
	offset = sf_input("offset");
	sf_floatread (off,nh,offset);
	sf_fileclose(offset);
    } else {
	if (!sf_histfloat(cmp,"d2",&dh)) sf_error("No d2= in input");
	if (!sf_histfloat(cmp,"o2",&h0)) sf_error("No o2= in input");

	if (!sf_getbool("half",&half)) half=true;
	/* if y, the second axis is half-offset instead of full offset */
	
	if (half) {
	    dh *= 2.;
	    h0 *= 2.;
	}
	
	if (sf_histfloat(cmp,"d3",&dy)) {
	    CDPtype=0.5+0.5*dh/dy;
	    if (1 != CDPtype) sf_histint(cmp,"CDPtype",&CDPtype);
	} 	 
	if (CDPtype < 1) CDPtype=1;
	sf_warning("CDPtype=%d",CDPtype);

	for (ih = 0; ih < nh; ih++) {
	    off[ih] = h0 + ih*dh; 
	}
    }

    nx = sf_leftsize(cmp,2);

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    trace = sf_floatalloc(nt);
    p = sf_floatalloc(nt);
    str = sf_floatalloc(nt);
    out = sf_floatalloc(nt);
    vtr = sf_floatalloc(nt);
    
    if (NULL != sf_getstring("crv")) {
	crv = sf_input("crv");
	eta = sf_output("eta");
	q = sf_floatalloc(nt);
	etr = sf_floatalloc(nt);
    } else {
	crv = NULL;
	eta = NULL;
	q = NULL;
	etr = NULL;
    }

    nmo = sf_stretch4_init (nt, t0, dt, nt, eps);

    eps = 100.*FLT_EPSILON;
    
    /* BUG: figure out dh for irregular off[ih] */

    for (ix = 0; ix < nx; ix++) { /* midpoints */
	for (ih = 0; ih < nh; ih++) { /* offset */
	    h = off[ih] + 0.5*dh + (dh/CDPtype)*(ix%CDPtype); 

	    sf_floatread (trace,nt,cmp); /* data */
	    sf_floatread (p, nt, dip);   /* slope */

	    if (NULL != crv)  {
		sf_floatread (q, nt, crv); /* curvature */

		for (it=0; it < nt; it++) {
		    t = t0 + it*dt;
		    f = fabsf(p[it]);
		    g = q[it]*h*dh/(f+eps);
		    if (g < 0.) {
			str[it]=t0-10.*dt;
			vtr[it]=0.;
			etr[it]=0.;
		    } else {
			str[it] = t - f*h*dt/(sqrtf(g)+dh);
			vtr[it] = 
			    sqrtf(fabsf(sqrtf(g)*h/
					(dt*(f*str[it]+eps))));
			etr[it] = 
			    - ((g-dh*dh)*t/(dt*(f*h+eps*dh))+dh)/
			    (8.*sqrtf(g)+eps*dh);
		    }
		}
	    } else {
		for (it=0; it < nt; it++) { /* time */
		    t = t0 + it*dt;
		    f = t - p[it]*h*dt/dh; 
		    
		    if (f < 0. || f > t) {
			str[it] = t0-10.*dt;
			vtr[it] = 0.;
		    } else {
			str[it] = sqrtf(t*f); /* t -> tau */
			vtr[it] = h/sqrtf(t*(t-f)+eps);
		    }
		}
	    }		

	    sf_stretch4_define (nmo,str);

	    sf_stretch4_apply (false,nmo,trace,out);	    
	    sf_floatwrite (out,nt,nmod);

	    sf_stretch4_apply (false,nmo,vtr,vtr);	    
	    sf_floatwrite (vtr,nt,vel);

	    if (NULL != crv)  {
		sf_stretch4_apply (false,nmo,etr,etr);	    
		sf_floatwrite (etr,nt,eta);
	    }
	}
    }


    exit (0);
}

/* 	$Id$	 */
