/* Constant-velocity stacks. */
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

static float v2;

static float nmo_map(float t, int it) {
    t = t*t + v2;
    if (t > 0.) t=sqrtf(t);
    return t;
}

int main (int argc, char* argv[])
{
    fint1 nmo;
    bool half, squared, slow;
    int ix,ih,it,iv, nt,nx,nw, nh, nh2, m, CDPtype, mute, *mask, nv, *fold;
    float dt, t0, h0,dh,h, dy, str, v0,dv,v;
    float **traces, *trace, *off, *stack;
    double *dstack;
    mapfunc nmofunc;
    sf_file cmp, stk, offset, msk;

    sf_init (argc,argv);
    cmp = sf_input("in");
    stk = sf_output("out");
    nmofunc =  nmo_map;
    
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

    if (!sf_getint("nv",&nv)) sf_error("Need nv=");
    /* number of velocities */
    
    if (!sf_getfloat("v0",&v0)) sf_error("Need v0=");
    /* first velocity */
    
    if (!sf_getfloat("dv",&dv)) sf_error("Need dv=");
    /* step in velocity */

    sf_putint(stk,"n2",nv);
    sf_putfloat(stk,"o2",v0);
    sf_putfloat(stk,"d2",dv);

    CDPtype=1;
    if (NULL != sf_getstring("offset")) {
	offset = sf_input("offset");
	if (SF_FLOAT != sf_gettype(offset)) sf_error("Need float offset");
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
	if (SF_INT != sf_gettype(msk)) sf_error("Need integer mask");
	nh2 = sf_filesize(msk);
	if (nh2 != nh && nh2 != nh*nx) sf_error("Wrong dimensions in mask");
	mask = sf_intalloc(nh2);
	sf_intread (mask,nh2,msk);
	sf_fileclose(msk);
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

    traces = sf_floatalloc2(nt,nh);
    trace  = sf_floatalloc(nt);
    stack  = sf_floatalloc(nt);
    fold   = sf_intalloc(nt);
    dstack = (double*) sf_alloc(nt,sizeof(double));

    nmo = fint1_init (nw, nt, mute);
    
    for (ix = 0; ix < nx; ix++) {
	sf_warning("CMP %d of %d;",ix+1,nx);
	sf_floatread (traces[0],nt*nh,cmp);

	for (iv=0; iv < nv; iv++) {
	    v = v0+iv*dv;
	    if (!squared) v *=v;

	    for (it=0; it < nt; it++) {
		dstack[it] = 0.0;
		fold[it] = 0;
	    }

	    for (ih = 0; ih < nh; ih++) {
		/* skip dead traces */
		if (NULL != msk) {
		    m = (nh2 == nh)? mask[ih] + (dh/CDPtype)*(ix%CDPtype) : 
			mask[ix*nh+ih];	
		    if (0==m) continue;
		}
		
		for (it=0; it < nt; it++) {
		    trace[it] = traces[ih][it];
		}
	    
		fint1_set(nmo,trace);

		h = (nh2 == nh)? off[ih] + (dh/CDPtype)*(ix%CDPtype) : 
		    off[ix*nh+ih];
		if (half) h *= 2;
		h = h*h - h0*h0;
		v2 = slow ? h*v : h/v;

		stretch(nmo,nmofunc,nt,dt,t0,nt,dt,t0,trace,str);

		for (it=0; it < nt; it++) {
		    if (trace[it] != 0.0f) {
			fold[it]++;
			dstack[it] += trace[it];
		    }
		}
	    }

	    for (it=0; it < nt; it++) {
		if (fold[it] > 0) {
		    stack[it] = dstack[it]/fold[it];
		} else {
		    stack[it] = 0.0f;
		}
	    }
			
	    sf_floatwrite (stack,nt,stk);
	}
    }
    sf_warning(".");

    exit (0);
}

