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

#include "stretch4.h"

int main (int argc, char* argv[])
{
    map4 nmo;
    bool half;
    int it,ix,ih, nt,nx, nh, nw, CDPtype;
    float dt, t0, h, h0, t, f, g, dh, eps, dy;
    float *trace, *p, *q, *off, *str, *out;
    sf_file cmp, nmod, dip, offset, crv;

    sf_init (argc,argv);
    cmp = sf_input("in");
    dip = sf_input("dip");
    nmod = sf_output("out");

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

    if (NULL != sf_getstring("crv")) {
	crv = sf_input("crv");
	q = sf_floatalloc(nt);
    } else {
	crv = NULL;
	q = NULL;
    }

    if (!sf_getint("extend",&nw)) nw=8;
    /* trace extension */

    nmo = stretch4_init (nt, t0, dt, nt, eps);
    
    /* BUG: figure out dh for irregular off[ih] */

    for (ix = 0; ix < nx; ix++) {
	for (ih = 0; ih < nh; ih++) {
	    h = off[ih] + 0.5*dh + (dh/CDPtype)*(ix%CDPtype); 

	    sf_floatread (trace,nt,cmp);
	    sf_floatread (p, nt, dip);

	    if (NULL != crv)  {
		sf_floatread (q, nt, crv);

		for (it=0; it < nt; it++) {
		    t = t0 + it*dt;
		    f = p[it] - 0.5*q[it];
		    g = q[it]*h*dh/(f+FLT_EPSILON);
		    if (g < 0.) {
			str[it]=t0-10.*dt;
		    } else {
			str[it] = t - f*h*dt/(dh+sqrtf(g));
		    }
		}
	    } else {
		for (it=0; it < nt; it++) {
		    t = t0 + it*dt;
		    f = t - p[it]*h*dt/dh;
		    
		    if (f < 0.) {
			str[it]=t0-10.*dt;
		    } else {
			str[it] = sqrtf(t*f);
		    }
		}
	    }		

	    stretch4_define (nmo,str);
	    stretch4_apply (nmo,trace,out);
	    
	    sf_floatwrite (out,nt,nmod);
	}
    }
	
    exit (0);
}

/* 	$Id: Minmo.c 729 2004-07-29 18:22:16Z fomels $	 */
