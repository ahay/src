/* Slopes for constant-velocity normal moveout. */
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
#include <rsf.h>

int main (int argc, char* argv[])
{
    bool half, squared, slow;
    int ix,ih,it,iv, nt,nx,nw, nh, nh2, CDPtype, nv;
    float dt, t0, h0,dh,h, dy, v0,dv,v, v2, delh, hp=0.0f, t;
    float *trace, *off;
    sf_file cmp, dip, offset;

    sf_init (argc,argv);
    cmp = sf_input("in");
    dip = sf_output("out");
    
    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    if (!sf_histint  (cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histint  (cmp,"n2",&nh)) sf_error("No n2= in input");
    nx = sf_leftsize(cmp,2);

    if (!sf_getbool("half",&half)) half=true;
    /* if y, the second axis is half-offset instead of full offset */

    if (!sf_getint("nv",&nv)) sf_error("Need nv=");
    /* number of velocities */
    
    if (!sf_getfloat("v0",&v0)) sf_error("Need v0=");
    /* first velocity */
    
    if (!sf_getfloat("dv",&dv)) sf_error("Need dv=");
    /* step in velocity */

    sf_shiftdim(cmp, dip, 2);

    sf_putint(dip,"n2",nv);
    sf_putfloat(dip,"o2",v0);
    sf_putfloat(dip,"d2",dv);

    CDPtype=1;
    if (NULL != sf_getstring("offset")) {
	offset = sf_input("offset");
	nh2 = sf_filesize(offset);
	if (nh2 != nh && nh2 != nh*nx) sf_error("Wrong dimensions in offset");

	off = sf_floatalloc(nh2);	
	sf_floatread (off,nh2,offset);
	sf_fileclose(offset);
    } else {
	if (!sf_histfloat(cmp,"d2",&dh)) sf_error("No d2= in input");
	if (!sf_histfloat(cmp,"o2",&h0)) sf_error("No o2= in input");

	if (sf_histfloat(cmp,"d3",&dy) && !sf_getint("CDPtype",&CDPtype)) {
	    CDPtype=half? 0.5+dh/dy : 0.5+0.5*dh/dy;
	    if (CDPtype < 1) {
		CDPtype=1;
	    } else if (1 != CDPtype) {
		sf_histint(cmp,"CDPtype",&CDPtype);
	    	sf_warning("CDPtype=%d",CDPtype);
	    }
	} 	    

	nh2 = nh;
	off = sf_floatalloc(nh2);
	for (ih = 0; ih < nh; ih++) {
	    off[ih] = h0 + ih*dh; 
	}

	offset = NULL;
    }

    if (!sf_getbool("slowness",&slow)) slow=false;
    /* if y, use slowness instead of velocity */

    if (!sf_getbool("squared",&squared)) squared=false;
    /* if y, the slowness or velocity is squared */

    if (!sf_getfloat ("h0",&h0)) h0=0.;
    /* reference offset */
    if (half) h0 *= 2.;
    if (!sf_getint("extend",&nw)) nw=4;
    /* trace extension */

    trace  = sf_floatalloc(nt);
    
    for (ix = 0; ix < nx; ix++) {
	for (iv=0; iv < nv; iv++) {
	    v = v0+iv*dv;
	    if (!squared) v *=v;

	    for (ih = 0; ih < nh; ih++) {
		h = (nh2 == nh)? off[ih] + (dh/CDPtype)*(ix%CDPtype) : 
		    off[ix*nh+ih];
		if (half) h *= 2;		
		if (ih==0) {
		    delh=h0;
		} else {
		    delh=h-hp;
		}
		hp=h;
		v2 = slow ? h*v: h/v;

		for (it=0; it < nt; it++) {
		    t = t0+it*dt;
		    t = SF_MAX(sqrtf(t*t + h*v2),dt);
		    
		    trace[it] = v2*delh/(t*dt);
		}
		sf_floatwrite(trace,nt,dip);
	    }
	}
    }

    exit (0);
}

