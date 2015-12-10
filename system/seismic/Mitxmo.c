/* Forward and inverse normal moveout with interval velocity. */
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

#include "warp2.h"

int main (int argc, char* argv[])
{
    bool inv;
    int it,ip, nt,nx, np, iy, ny;
    float dt, t0, p, p0, p2, v, sx, st, ft, dp, eps, x0, dx;
    float **tp, *vel, **tx, **t, **x;
    sf_file inp, out, velocity;

    sf_init (argc,argv);

    inp = sf_input("in");
    out = sf_output("out");

    velocity = sf_input("velocity");


    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(inp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_getbool("inv",&inv)) inv=false;

    if (inv) {
	if (!sf_histint(inp,"n2",&np)) sf_error("No n2= in input");
	if (!sf_histfloat(inp,"d2",&dp)) sf_error("No d2= in input");
	if (!sf_histfloat(inp,"o2",&p0)) sf_error("No o2= in input");
	
	if (!sf_getint("nx",&nx)) sf_error("Need nx=");
	/* offset samples */
	if (!sf_getfloat("dx",&dx)) sf_error("Need dx=");
	/* offset sampling */
	if (!sf_getfloat("x0",&x0)) x0=0.;
	/* first offset */
	
	sf_putint(out,"n2",nx);
	sf_putfloat(out,"d2",dx);
	sf_putfloat(out,"o2",x0);
    } else {
	if (!sf_histint(inp,"n2",&nx)) sf_error("No n2= in input");
	if (!sf_histfloat(inp,"d2",&dx)) sf_error("No d2= in input");
	if (!sf_histfloat(inp,"o2",&x0)) sf_error("No o2= in input");
	
	if (!sf_getint("np",&np)) sf_error("Need np=");
	/* slope samples */
	if (!sf_getfloat("dp",&dp)) sf_error("Need dp=");
	/* slope sampling */
	if (!sf_getfloat("p0",&p0)) p0=0.;
	/* first slope */
	
	sf_putint(out,"n2",np);
	sf_putfloat(out,"d2",dp);
	sf_putfloat(out,"o2",p0);
    }

    ny = sf_leftsize(inp,2);

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    warp2_init(nt, t0, dt,
	       nx, x0, dx /* output grid */,
	       nt, np     /* input  grid */,
	       eps     /* regularization */);

    tp = sf_floatalloc2(nt,np);
    t = sf_floatalloc2(nt,np);
    x = sf_floatalloc2(nt,np);

    tx   = sf_floatalloc2(nt,nx);
    vel = sf_floatalloc(nt);

    for (iy = 0; iy < ny; iy++) {
	sf_floatread (vel,nt,velocity);	

	for (ip = 0; ip < np; ip++) {
	    p = p0 + ip*dp;
	    p2 = p*p;

	    st = 0.;
	    sx = 0.;
	    for (it=0; it < nt; it++) {
		v = vel[it];
		v *= v;

		ft = 1.-p2*v;

		if (ft < 0.) {
		    for (it=0; it < nt; it++) {
			t[ip][it]=t0-10.*dt;
			x[ip][it]=x0-10.*dx;
		    }
		    break;
		}

		ft = sqrtf(1.0/ft);

		if (it==0) {
		    st = t0*ft/dt;
		    sx = p*v*st;
		}

		t[ip][it] = st*dt;
		x[ip][it] = sx*dt;
		
		st += ft;
		sx += p*v*ft;
	    }
	}

	if (inv) {
	    sf_floatread (tp[0],nt*np,inp);	    
	    warp2(tp,t,x,tx);	    
	    sf_floatwrite (tx[0],nt*nx,out);
	} else {
	    sf_floatread (tx[0],nt*nx,inp);
	    fwarp2(tx,t,x,tp);
	    sf_floatwrite (tp[0],nt*np,out);
	}
    }

    exit (0);
}
