/* Phase-space Radon transform */
/*
  Copyright (C) 2007 University of Texas at Austin
  
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

#include "triangle2.h"

int main(int argc, char* argv[])
{
    bool verb;
    int nt, nx, np, n3, i, it, ix, ip, i3, ntx, ntp;
    int interp, rect1, rect2, niter;
    float t0, dt, t, x0, dx, x, p0, dp, p, eps;
    float *dtx, *dtp, **tp, **tx, *pp;
    sf_file in, dip, out, dip2;

    sf_init(argc,argv);
    in = sf_input("in");
    dip = sf_input("dip");
    out = sf_output("out");

    /* allocations */
    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    ntx = nt*nx;
 
    if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o2",&x0)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");

    if (!sf_getint("np",&np)) sf_error("Need np=");
    /* number of slopes */
    if (!sf_getfloat("p0",&p0)) sf_error("Need p0=");
    /* first slope */
    if (!sf_getfloat("dp",&dp)) sf_error("Need dp=");
    /* slope increment */
    ntp = nt*np;

    if (!sf_getint("interp",&interp)) interp=2;
    /* interpolation length */

    if (!sf_getint("rect1",&rect1)) rect1=3;
    if (!sf_getint("rect2",&rect2)) rect2=3;
    /* smoothing length for shaping */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */

    if (!sf_getfloat("eps",&eps)) eps=1./ntx;
    /* regularization parameter */

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    triangle2_init(rect1, rect2, nt, np, 1);
    sf_conjgrad_init(ntp, ntp, ntx, ntx, eps, FLT_EPSILON, verb, false);

    sf_putint(out,"n2",np);
    sf_putfloat(out,"o2",p0);
    sf_putfloat(out,"d2",dp);

    if (NULL != sf_getstring("dip2")) {
	dip2 = sf_output("dip2");
	sf_putint(dip2,"n2",np);
	sf_putfloat(dip2,"o2",p0);
	sf_putfloat(dip2,"d2",dp);
    } else {
	dip2 = NULL;
    }

    dtx = sf_floatalloc(ntx);
    dtp = sf_floatalloc(ntp);
    pp  = sf_floatalloc(ntp);
    tp = sf_floatalloc2(2,ntx);
    tx = sf_floatalloc2(2,ntp);

    for (i3=0; i3 < n3; i3++) {
	/* push slope */
	sf_floatread(dtx,ntx,dip);
	
	for (i=ix=0; ix < nx; ix++) {
	    x = x0+ix*dx;
	    for (it=0; it < nt; it++, i++) {
		t = t0+it*dt;

		p = dtx[i]*dt/dx;
		dtx[i] = -x*dp/dt;

		tp[i][0] = t-p*x;
		tp[i][1] = p;
	    }
	}

	sf_int2_init(tp, t0,p0, dt,dp, nt,np, sf_spline_int, interp, ntx);
	sf_conjgrad(NULL, sf_int2_lop, triangle2_lop, 
		    pp, dtp, dtx, niter);
	if (NULL != dip2) sf_floatwrite(dtp,ntp,dip2);
	sf_int2_close();

	/* pull data */
	sf_floatread(dtx,ntx,in);
	for (i=ip=0; ip < np; ip++) {
	    p = p0+ip*dp;
	    for (it=0; it < nt; it++, i++) {
		t = t0+it*dt;
		x = -dtp[i]*dt/dp;

		tx[i][0] = t+p*x;
		tx[i][1] = x;
	    }
	}
	sf_int2_init(tx, t0,x0, dt,dx, nt,nx, sf_spline_int, interp, ntp);
	sf_int2_lop(false,false,ntx,ntp,dtx,dtp);
	sf_floatwrite(dtp,ntp,out);
	sf_int2_close();
    }

    exit(0);
}
