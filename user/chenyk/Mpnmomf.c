/* Normal moveout.

Compatible with sfvscan.

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

static float *vel, *par, h;
static bool slow, squared;

static float nmo_map(float t, int it) {
    float v;
    
    v = vel[it];
    if (!squared) v *= v;
    v = slow ? h*v : h/v;
    t = t*t + v;
    if (t > 0.) t=sqrtf(t);
    return t;
}

static float shifted_nmo_map(float t, int it) {
    float v, q, s;
    
    v = vel[it];
    s = par[it];

    if (!squared) v *= v;
    v = slow ? h*v : h/v;
    q = t*t + s*v;
    if (q > 0.) {
	t=t*(1.-1./s)+sqrtf(q)/s;
	return t;
    } else {
	return q;
    }
}

static float taner_nmo_map(float t, int it) {
    float v, a;
    
    v = vel[it];
    a = par[it];

    if (!squared) v *= v;
    v = slow ? h*v : h/v;
    t = t*t + v/(1+a*v);
    if (t > 0.) t=sqrtf(t);
    return t;
}

int main (int argc, char* argv[])
{
    fint1 nmo;
    bool half;
    int ix,ih, nt,nx,nw, nh, nh2, CDPtype, mute, *mask;
    float dt, t0, h0, dh, dy, str;
    float *trace, *off;
    mapfunc nmofunc;
    sf_file cmp, nmod, velocity, offset, msk, het;

    sf_init (argc,argv);
    cmp      = sf_input("in");
    velocity = sf_input("velocity");
    nmod     = sf_output("out");

    if (NULL != sf_getstring("s")) {
      	het = sf_input("s");
	nmofunc = shifted_nmo_map;	
    } else if (NULL != sf_getstring("a")) {
        het = sf_input("a");
        nmofunc = taner_nmo_map;
    } else {
	het = NULL;
	nmofunc =  nmo_map;
    } 

    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    if (!sf_histint  (cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histint  (cmp,"n2",&nh)) sf_error("No n2= in input");
    nx = sf_leftsize(cmp,2);

    if (!sf_getbool("half",&half)) half=true;
    /* if y, the second axis is half-offset instead of full offset */
    if (!sf_getfloat("str",&str)) str=0.5;
    /* maximum stretch allowed */

    if (!sf_getint("mute",&mute)) mute=12;
    /* mute zone */

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
    
    if (NULL != sf_getstring("mask")) {
	msk = sf_input("mask");
	mask = sf_intalloc(nh);
    } else {
	msk = NULL;
	mask = NULL;
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

    trace = sf_floatalloc(nt);
    vel   = sf_floatalloc(nt);
    if (NULL != het) {
	par = sf_floatalloc(nt);
    } else {
	par = NULL;
    }

    nmo = fint1_init (nw, nt, mute);
    
    for (ix = 0; ix < nx; ix++) {
	sf_warning("CMP %d of %d;",ix+1,nx);

	sf_floatread (vel,nt,velocity);	
	if (NULL != het) sf_floatread(par,nt,het);
	if (NULL != msk) sf_intread(mask,nh,msk);

	for (ih = 0; ih < nh; ih++) {
	    sf_floatread (trace,nt,cmp);
	    
	    /* skip dead traces */
	    if (NULL != msk && 0==mask[ih]) {
		sf_floatwrite (trace,nt,nmod);
		continue;
	    }
	    
	    fint1_set(nmo,trace);

	    h = (nh2 == nh)? off[ih] + (dh/CDPtype)*(ix%CDPtype) : 
		off[ix*nh+ih];
	    if (half) h *= 2;
	    h = h*h - h0*h0;
	    
	    stretch(nmo,nmofunc,nt,dt,t0,nt,dt,t0,trace,str);
	    sf_floatwrite (trace,nt,nmod);
	}
    }
    sf_warning(".");

    exit (0);
}

