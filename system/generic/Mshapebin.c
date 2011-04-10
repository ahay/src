/* Data binning in 2-D slices by inverse interpolation. */
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

#include <float.h>
#include <math.h>

#include <rsf.h>

#include "laplac2.h"
#include "gaussshape2.h"

int main (int argc, char* argv[])
{
    bool gauss, shape, verb;
    int id, nk, nd, im, nm, nt, it, nx, ny, n2, xkey, ykey;
    int interp, niter, nliter;
    float *pp, *mm, *mm0=NULL, *dd, **xy, *hdr, filt1, filt2, a[3];
    float x0, y0, dx, dy, xmin, xmax, ymin, ymax, f, dt, t0, eps;
    char *header;
    sf_file in, out, head, pattern=NULL, pin=NULL, pout=NULL;
    sf_operator shaping=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nd)) sf_error("Need n1= in in");
    if (!sf_histint(in,"n2",&nt)) sf_error("Need n2= in in");
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_getint("xkey",&xkey)) sf_error("Need xkey=");
    /* x key number */
    if (!sf_getint("ykey",&ykey)) sf_error("Need ykey=");
    /* y key number */

    /* create coordinates */
    xy = sf_floatalloc2(2,nd);

    header = sf_getstring("head");
    if (NULL == header) { 
	header = sf_histstring(in,"head");
	if (NULL == header) sf_error("Need head=");
    }

    head = sf_input(header);

    if (SF_FLOAT != sf_gettype(head)) sf_error("Need float header");
    if (!sf_histint(head,"n1",&nk)) sf_error("No n1= in head");
    if (!sf_histint(head,"n2",&n2) || n2 != nd) 
	sf_error("Wrong n2= in head");

    hdr = sf_floatalloc(nk);

    ymin = xmin = +FLT_MAX;
    ymax = xmax = -FLT_MAX;
    for (id=0; id<nd; id++) {	
	sf_floatread (hdr,nk,head);
	f = hdr[xkey]; 
	if (f < xmin) xmin=f;
	if (f > xmax) xmax=f;
	xy[id][0] = f;
	f = hdr[ykey]; 
	if (f < ymin) ymin=f;
	if (f > ymax) ymax=f;
	xy[id][1] = f;
    }

    sf_fileclose (head);

    /* create model */
    if (!sf_getint ("nx",&nx)) sf_error("Need nx=");
    /* Number of bins in x */
    if (!sf_getint ("ny",&ny)) sf_error("Need ny=");
    /* Number of bins in y */

    sf_putint(out,"n1",nx);
    sf_putint(out,"n2",ny);
    sf_putint(out,"n3",nt);

    if (sf_histfloat(in,"o2",&t0)) sf_putfloat(out,"o3",t0);
    if (sf_histfloat(in,"d2",&dt)) sf_putfloat(out,"d3",dt);

    /* let user overwrite */
    sf_getfloat ("xmin",&xmin);
    sf_getfloat ("xmax",&xmax);
    sf_getfloat ("ymin",&ymin);
    /* Grid dimensions */
    sf_getfloat ("ymax",&ymax);

    if (xmax <= xmin) sf_error ("xmax=%f <= xmin=%f",xmax,xmin);
    if (ymax <= ymin) sf_error ("ymax=%f <= ymin=%f",xmax,xmin);

    if (!sf_getfloat("x0",&x0)) x0=xmin; 
    if (!sf_getfloat("y0",&y0)) y0=ymin; 
    /* grid origin */

    sf_putfloat (out,"o1",x0);
    sf_putfloat (out,"o2",y0);

    if (!sf_getfloat("dx",&dx)) {
	/* bin size in x */
	if (1 >= nx) sf_error("Need dx=");
	dx = (xmax-xmin)/(nx-1);
    }

    if (!sf_getfloat("dy",&dy)) {
	/* bin size in y */
	if (1 >= nx) {
	    dy = dx;
	} else {
	    dy = (ymax-ymin)/(ny-1);
	}
    }

    sf_putfloat (out,"d1",dx);
    sf_putfloat (out,"d2",dy);
    
    /* initialize interpolation */
    if (!sf_getint("interp",&interp)) interp=2;
    /* interpolation length */

    sf_int2_init (xy, x0,y0, dx,dy, nx,ny, sf_spline_int, interp, nd);

    nm = nx*ny;
    mm = sf_floatalloc(nm);
    dd = sf_floatalloc(nd);

    if (!sf_getint("niter",&niter)) niter=nm;
    /* number of iterations */
    if (!sf_getint("nliter",&nliter)) nliter=1;
    /* number of reweighting iterations */

    if (!sf_getfloat("eps",&eps)) eps=1./nd;
    /* regularization parameter */
    if (!sf_getbool("gauss",&gauss)) gauss=false;
    /* if y, use gaussian shaping (estimated from the data) */
    if (!sf_getbool("shape",&shape)) shape=true;
    /* if y, use shaping regularization */
    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */

    if (shape) {
	if (!sf_getfloat("filt1",&filt1)) filt1=3.;
	if (!sf_getfloat("filt2",&filt2)) filt2=filt1;
	/* smoothing length for shaping */

	if (filt1 < 1. || filt2 < 1.) 
	    sf_error("wrong filt1=%g or filt2=%g",filt1,filt2);

	if (gauss) {
	    if (!sf_getfloat("a0",a))   a[0] = (filt1*filt1-1.)/12.;
	    if (!sf_getfloat("b0",a+1)) a[1] = 0.;
	    if (!sf_getfloat("c0",a+2)) a[2] = (filt2*filt2-1.)/12.;
	    /* initial Gaussian shape parameters */

	    gaussshape2_init(nx,ny);

	    if (sf_getstring("pattern")) {
		pattern = sf_input("pattern");
	    } else { /* predefined Gaussian */
		pattern = NULL;
		gaussshape2_set2(a);
	    }
	    shaping = sf_freqfilt2_lop;
	} else {
	    sf_triangle2_init((int) filt1, (int) filt2, nx, ny, 1);
	    shaping = sf_triangle2_lop;
	}

	pp = sf_floatalloc(nm);

	if (NULL != sf_getstring("pin")) {
	    pin = sf_input("pin");
	    mm0 = sf_floatalloc(nm);
	} 
	
	if (NULL != sf_getstring("pout")) {
	    pout = sf_output("pout");
	    sf_putint(pout,"n1",nx);
	    sf_putint(pout,"n2",ny);
	    sf_putint(pout,"n3",nt);
	    sf_putfloat (pout,"o1",x0); sf_putfloat (pout,"d1",dx);
	    sf_putfloat (pout,"o2",y0); sf_putfloat (pout,"d2",dy);
	    if (sf_histfloat(in,"o2",&t0)) sf_putfloat(pout,"o3",t0);
	    if (sf_histfloat(in,"d2",&dt)) sf_putfloat(pout,"d3",dt);
	}
    } else {
	laplac2_init(nx,ny);
	pp = NULL;
    }

    sf_conjgrad_init(nm, nm, nd, nd, eps, 1.e-9, verb, false);

    for (it=0; it < nt; it++) { /* loop over time slices */
	sf_warning("slice %d of %d",it+1,nt);
	sf_floatread (dd,nd,in);

	if (shape) {
	    if (gauss) {
		if (NULL != pattern) {
		    /* estimate shaper */
		    sf_floatread (mm,nm,pattern);
		    gaussshape2_set(a, mm, 100, nliter);
		}
	    }

	    if (NULL != pin) {
		sf_floatread(pp,nm,pin);
		for (im=0; im < nm; im++) {
		    pp[im] = -pp[im];
		}
		shaping(false,false,nm,nm,pp,mm0);
		sf_int2_lop(false,true,nm,nd,mm0,dd);
	    }
	    
	    sf_conjgrad(NULL, sf_int2_lop, shaping, pp, mm, dd, niter);

	    if (NULL != pin) {
		for (im=0; im < nm; im++) {
		    mm[im] -= mm0[im];
		}
	    }
	} else {
	    sf_solver_reg(sf_int2_lop,sf_cgstep,laplac2_lop,
			  nm,nm,nd,mm,dd,niter,eps,"end");
	    sf_cgstep_close();
	}

	sf_floatwrite (mm,nm,out);
	if (NULL != pout) sf_floatwrite(pp,nm,pout);
    }


    exit(0);
}

/* 	$Id$	 */
