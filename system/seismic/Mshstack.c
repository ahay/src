/* Shaping stack. */
/*
  Copyright (C) 2015 University of Texas at Austin
  
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

#include "inmo.h"

int main (int argc, char* argv[])
{
    bool half, slow;
    int ix,ih, nd, it,nt,nx, nh, CDPtype, jump, niter, restart;
    float dt, t0, h0, dh, eps, dy, tol;
    float *trace, *trace2, *vel, *off, **gather, **dense;
    sf_file cmp, stack, velocity, offset;

    sf_init (argc,argv);
    cmp = sf_input("in");
    velocity = sf_input("velocity");
    stack = sf_output("out");

    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histint(cmp,"n2",&nh)) sf_error("No n2= in input");
    
    sf_putint(stack,"n2",1);

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

    if (!sf_getint("jump",&jump)) jump=1;
    /* subsampling */

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of iterations */

    if (!sf_getint("restart",&restart)) restart=niter;
    /* GMRES memory */

    if (!sf_getfloat("tol",&tol)) tol=1e-5;
    /* GMRES tolerance */

    nd = (nt-1)*jump+1;
    
    sf_putint(stack,"n1",nd);
    sf_putfloat(stack,"d1",dt/jump);

    trace = sf_floatalloc(nd);
    trace2 = sf_floatalloc(nd);
    gather = sf_floatalloc2(nt,nh);
    dense = sf_floatalloc2(nd,nh);

    vel = sf_floatalloc(nd);

    for (ix = 0; ix < nx; ix++) { /* loop over midpoint */
	sf_floatread (vel,nd,velocity);	
	
	inmo_init(vel, off, nh, 
		  h0, dh, CDPtype, ix,
		  nt, slow,
		  t0, dt, eps, half,
		  jump);

	sf_floatread (gather[0],nt*nh,cmp);

	/* apply backward operator */
	interpolate(gather, dense);
	nmostack(dense,trace);
	
	sf_gmres_init(nd,restart);

	for (it=0; it < nd; it++) {
	    trace2[it] = 0.0f;
	}

	/* run GMRES */
	sf_gmres(trace,trace2,inmo_oper,NULL,niter,tol,true);

	sf_floatwrite (trace2,nd,stack);

	sf_gmres_close();

	inmo_close();
    }


    exit (0);
}


