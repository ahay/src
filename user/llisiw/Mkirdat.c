/* 2-D Pre-stack Kirchhoff redatuming. */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "kirdat.h"

int main(int argc, char* argv[])
{
    bool verb;
    int it, nt, ih, nh, is, ns, ng, left, right, ic, aper, shift, c, cc, hh;
    int ir, nr, sleft, sright;
    float datum, length, t0, dt, h0, dh, s0, ds, g0, dg, dist, tau, delta;
    float dr, s, h;
    float ***tr_in, ***tr_out, **table;
    sf_file in, out, green;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (!sf_getfloat("datum",&datum)) sf_error("Need datum=");
    /* datum depth */

    if (!sf_getint("aperture",&aper)) aper=50;
    /* aperture (number of near traces) */

    if (!sf_getfloat("length",&length)) length=0.025;
    /* filter length (in seconds) */

    /* read input */
    if (!sf_histint(in,"n1",&nt)) sf_error("No nt=");
    if (!sf_histint(in,"n2",&nh)) sf_error("No nh=");
    if (!sf_histint(in,"n3",&ns)) sf_error("No ns=");

    if (!sf_histfloat(in,"o1",&t0)) sf_error("No t0=");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No dt=");

    if (!sf_histfloat(in,"o2",&h0)) sf_error("No h0=");
    if (!sf_histfloat(in,"d2",&dh)) sf_error("No dh=");

    if (!sf_histfloat(in,"o3",&s0)) sf_error("No s0=");
    if (!sf_histfloat(in,"d3",&ds)) sf_error("No ds=");

    tr_in = sf_floatalloc3(nt,nh,ns);
    sf_floatread(tr_in[0][0],nt*nh*ns,in);

    /* allocate memory for output */
    tr_out = sf_floatalloc3(nt,nh,ns);

    /* read Green's function (traveltime table) */
    green = sf_input("table");

    if (!sf_histint(green,"n1",&ng)) sf_error("No ng=");
    if (!sf_histfloat(green,"o1",&g0)) sf_error("No g0=");
    if (!sf_histfloat(green,"d1",&dg)) sf_error("No dg=");

    table = sf_floatalloc2(ng,ng);
    sf_floatread(table[0],ng*ng,green);
    sf_fileclose(green);

    /* initialize */
    filt_init(dt,length);

    /* common-shot gather */
#ifdef _OPENMP
#pragma omp parallel for private(c,left,right,ic,cc,tau,dist,shift,it,delta)
#endif
    for (is=0; is < ns; is++) {
	if (verb) sf_warning("Processing common-shot gather %d of %d.",is+1,ns);

	for (ih=0; ih < nh; ih++) {

	    c = (s0+is*ds+h0+ih*dh-g0)/dg+0.5;
	    if (c < 0 || c > ng-1) sf_error("Traveltime table too small.");

	    /* aperture */
	    left  = (ih-aper < 0)?    0:    ih-aper;
	    right = (ih+aper > nh-1)? nh-1: ih+aper;
	    
	    for (ic=left; ic <= right; ic++) {
		
		cc = (s0+is*ds+h0+ic*dh-g0)/dg+0.5;
		if (cc < 0 || cc > ng-1) sf_error("Traveltime table too small.");

		/* time delay */
		tau = table[cc][c];

		/* distance */
		dist = datum*datum+(ic-ih)*dh*(ic-ih)*dh;

		/* filter (tau dependent) */
		filt_set(tau);

		shift = 0;
		delta = 0.;
		for (it=0; it < nt; it++) {
		    if (((float)it)*dt < tau) 
			continue;
		    else if (shift == 0)
			delta = (((float)it*dt)-tau)/dt;

		    tr_out[is][ih][it] += 1./SF_PI
			*dh*datum*tau/dist
			*pick(delta,tr_in[is][ic],shift);
		    shift++;
		}
	    }

	}
    }

    /* zero input */
    for (is=0; is < ns; is++) {
	for (ih=0; ih < nh; ih++) {
	    for (it=0; it < nt; it++) {
		tr_in[is][ih][it] = 0.;
	    }
	}
    }

    /* acquisition */
    s = (ns-1)*ds;
    h = (nh-1)*dh;

    dr = (fabsf(ds)>=fabsf(dh))? fabsf(dh): fabsf(ds);
    nr = (fabsf(s)+fabsf(h))/dr+1.5;

    /* common-receiver gather */
#ifdef _OPENMP
#pragma omp parallel for private(sleft,sright,is,c,ih,left,right,ic,cc,hh,tau,dist,shift,it,delta)
#endif
    for (ir=0; ir < nr; ir++) {
	if (verb) sf_warning("Processing common-receiver gather %d of %d.",ir+1,nr);

	sleft  = (ir*dr+((ds<=0.)?1.:0.)*s+((dh<=0.)?1.:-1.)*h)/ds+0.5;
	sright = (ir*dr+((ds<=0.)?1.:0.)*s)/ds+0.5;

	if (sleft < 0) sleft = 0;
	if (sright > ns-1) sright = ns-1;
	
	for (is=sleft; is <= sright; is++) {
	    
	    c = (s0+is*ds-g0)/dg+0.5;
	    if (c < 0 || c > ng-1) sf_error("Traveltime table too small.");

	    ih = (ir*dr-is*ds+((ds<=0.)?1.:0.)*s+((dh<=0.)?1.:0.)*h)/dh+0.5;
	    
	    /* aperture */
	    left  = (is-aper < sleft)?  sleft:  is-aper;
	    right = (is+aper > sright)? sright: is+aper;
	    
	    for (ic=left; ic <= right; ic++) {
		
		cc = (s0+ic*ds-g0)/dg+0.5;
		if (cc < 0 || cc > ng-1) sf_error("Traveltime table too small.");

		hh = (ir*dr-ic*ds+((ds<=0.)?1.:0.)*s+((dh<=0.)?1.:0.)*h)/dh+0.5;
		
		/* time delay */
		tau = table[cc][c];
		
		/* distance */
		dist = datum*datum+(ic-is)*ds*(ic-is)*ds;
		
		/* filter (tau dependent) */
		filt_set(tau);
		
		shift = 0;
		for (it=0; it < nt; it++) {
		    if (((float)it)*dt < tau) 
			continue;
		    else if (shift == 0)
			delta = (((float)it*dt)-tau)/dt;
		    
		    tr_in[is][ih][it] += 1./SF_PI
			*ds*datum*tau/dist
			*pick(delta,tr_out[ic][hh],shift);
		    shift++;
		}
	    }
	}	
    }

    /* write output */
    sf_floatwrite(tr_in[0][0],nt*nh*ns,out);

    exit(0);
}
