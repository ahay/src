/* Inverse normal moveout. */
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
    bool half, slow;
    int it,ix,ih, nt,nx, nh, CDPtype;
    float dt, t0, h, h0, f, dh, eps, dy;
    float *trace, *vel, *off, *str, *out;
    sf_file cmp, nmod, velocity, offset;

    sf_init (argc,argv);
    cmp = sf_input("in");
    velocity = sf_input("velocity");
    nmod = sf_output("out");

    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(cmp,"n2",&nh)) sf_error("No n2= in input");

    off = sf_floatalloc(nh);

    if (!sf_getbool("half",&half)) half=true;
    /* if y, the second axis is half-offset instead of full offset */


    CDPtype=1;
    if (NULL != sf_getstring("offset")) {
	offset = sf_input("offset");
	sf_floatread (off,nh,offset);
	sf_fileclose(offset);
    } else {
	if (!sf_histfloat(cmp,"d2",&dh)) sf_error("No d2= in input");
	if (!sf_histfloat(cmp,"o2",&h0)) sf_error("No o2= in input");
	
	if (sf_histfloat(cmp,"d3",&dy)) {
	    CDPtype=half? 0.5+dh/dy : 0.5+0.5*dh/dy;
	    if (CDPtype < 1) {
		CDPtype=1;
	    } else if (1 != CDPtype) {
		sf_histint(cmp,"CDPtype",&CDPtype);
	    	sf_warning("CDPtype=%d",CDPtype);
	    }
	} 

	for (ih = 0; ih < nh; ih++) {
	    off[ih] = h0 + ih*dh; 
	}
    }

    if (!sf_getbool("slowness",&slow)) slow=false;
    /* if y, use slowness instead of velocity */

    nx = sf_leftsize(cmp,2);

    if (!sf_getfloat ("h0",&h0)) h0=0.;
    /* reference offset */
    if (half) h0 *= 2.;
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    trace = sf_floatalloc(nt);
    vel = sf_floatalloc(nt);
    str = sf_floatalloc(nt);
    out = sf_floatalloc(nt);

    nmo = stretch4_init (nt, t0, dt, nt, eps);
    
    for (ix = 0; ix < nx; ix++) {
	sf_floatread (vel,nt,velocity);	

	for (ih = 0; ih < nh; ih++) {
	    sf_floatread (trace,nt,cmp);
	    
	    h = off[ih] + (dh/CDPtype)*(ix%CDPtype); 
	    if (half) h *= 2;
	    h = h*h - h0*h0;
	    
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
	    stretch4_apply (false,nmo,trace,out);
	    
	    sf_floatwrite (out,nt,nmod);
	}
    }


    exit (0);
}

/* 	$Id$	 */
