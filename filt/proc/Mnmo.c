/* Normal moveout.

Compatible with sfvscan.
*/
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

#include "fint1.h"

int main (int argc, char* argv[])
{
    fint1 nmo;
    bool half, slow;
    int it,ix,iz,ih, nt,nx,nw, nh, nh2, CDPtype;
    float dt, t0, h, h0, f, dh, dy, v;
    float *trace, *vel, *off;
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
    nx = sf_leftsize(cmp,2);

    CDPtype=1;
    if (NULL != sf_getstring("offset")) {
	offset = sf_input("offset");
	if (!sf_histint(offset,"n1",&nh2)) sf_error("No n1= in offset");
	if (nh2 != nh && (nh2 != nx || (nh2 == nx && nh != 1)))
	    sf_error("Wrong dimensions in offset");

	off = sf_floatalloc(nh2);
	sf_floatread (off,nh2,offset);
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

	nh2 = nh;
	off = sf_floatalloc(nh2);
	for (ih = 0; ih < nh; ih++) {
	    off[ih] = h0 + ih*dh; 
	}	
    }

    if (!sf_getbool("slowness",&slow)) slow=false;
    /* if y, use slowness instead of velocity */

    if (!sf_getfloat ("h0",&h0)) h0=0.;
    /* reference offset */
    if (!sf_getint("extend",&nw)) nw=4;
    /* trace extension */

    trace = sf_floatalloc(nt);
    vel = sf_floatalloc(nt);

    nmo = fint1_init (nw, nt);
    
    for (ix = 0; ix < nx; ix++) {
	sf_floatread (vel,nt,velocity);	

	for (ih = 0; ih < nh; ih++) {
	    sf_floatread (trace,nt,cmp);
	    fint1_set(nmo,trace);

	    h = (nh2 == nh)? off[ih] + (dh/CDPtype)*(ix%CDPtype) : 
		off[ix]; 
	    h = h*h - h0*h0;
	    
	    for (it=0; it < nt; it++) {
		f = t0 + it*dt;
		v = vel[it];
		v = slow? h*v*v: h/(v*v);
		f = f*f + v;
		if (f < 0.) {
		    trace[it]=0.;
		} else {
		    f = (sqrtf(f) - t0)/dt;
		    iz = f;
		    if (iz >= 0 && iz < nt) {
			trace[it] = fint1_apply(nmo,iz,f-iz,false);
		    } else {
			trace[it] = 0.;
		    }
		}
	    }
	    
	    sf_floatwrite (trace,nt,nmod);
	}
    }
	
    exit (0);
}

/* 	$Id$	 */
