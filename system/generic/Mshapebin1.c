/* 1-D inverse interpolation with shaping regularization. */
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

#include "gauss.h"
#include "monofshape.h"

int main (int argc, char* argv[])
{
    bool pef, gauss, verb;
    int id, nd, nt, it, nx, interp, niter, pad, i3, n3;
    float *pp, *mm, *dd, *offset, x0, dx, xmin, xmax, f, filt, eps;
    char *header;
    sf_file in, out, head;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&nd)) nd=1;
    if (!sf_histint(in,"n2",&nt)) nt=1;
    n3 = sf_leftsize(in,2);

    /* create coordinates */
    offset = sf_floatalloc(nd);

    header = sf_getstring("head");
    if (NULL == header) { 
	header = sf_histstring(in,"head");
	if (NULL == header) sf_error("Need head=");
    }

    head = sf_input(header);
    if (SF_FLOAT != sf_gettype(head)) sf_error("Need float head");

    if (1==n3) {
	sf_floatread (offset,nd,head);
	
	xmin = +FLT_MAX;
	xmax = -FLT_MAX;
	for (id=0; id<nd; id++) {	
	    f = offset[id]; 
	    if (f < xmin) xmin=f;
	    if (f > xmax) xmax=f;
	}
    }

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */
 
    /* create model */
    if (!sf_getint ("nx",&nx)) sf_error("Need nx=");
    /* number of bins */

    sf_putint(out,"n1",nx);

    /* let user overwrite */
    sf_getfloat ("xmin",&xmin);
    /* grid size */
    sf_getfloat ("xmax",&xmax);
    if (xmax <= xmin) sf_error ("xmax=%f <= xmin=%f",xmax,xmin);

    if (!sf_getfloat("x0",&x0)) x0=xmin; 
    /* grid origin */
    sf_putfloat (out,"o1",x0);

    if (!sf_getfloat("dx",&dx)) {
	/* grid sampling */
	if (1 >= nx) sf_error("Need dx=");
	dx = (xmax-xmin)/(nx-1);
    }
    sf_putfloat (out,"d1",dx);
    
    /* initialize interpolation */
    if (!sf_getint("interp",&interp)) interp=2;
    /* [1,2] forward interpolation method, 1: nearest neighbor, 2: linear */

    if (!sf_getfloat("filter",&filt)) filt=3.;
    /* smoothing length */

    if (filt < 1.) sf_error("filt must be >= 1, got filt=%g",filt);

    pp = sf_floatalloc(nx);
    mm = sf_floatalloc(nx);
    dd = sf_floatalloc(nd);

    if (!sf_getint("niter",&niter)) niter=nx;
    /* number of conjugate-gradient iterations */
    if (!sf_getfloat("eps",&eps)) eps=1./nd;
    /* regularization parameter */


    if (!sf_getbool("pef",&pef)) pef=false;
    /* if y, use monofrequency regularization */
    if (!sf_getbool("gauss",&gauss)) gauss=false;
    /* if y, use Gaussian shaping */
    if (!sf_getint("pad",&pad)) pad=0;
    /* padding for Gaussian shaping */

    if (gauss) {
	gauss_init (nx+pad, filt);
    } else {
	sf_triangle1_init ((int) filt, nx);
    }
    sf_conjgrad_init(nx, nx, nd, nd, eps, 1.e-10, verb, false);

    if (pef) monofshape_init(nx+pad);

    for (i3=0; i3 < n3; i3++) {
	if (!verb) sf_warning("%d of %d;",i3+1,n3);

	if (n3 > 1) sf_floatread (offset,nd,head);

	switch (interp) {
	    case 1:
		sf_int1_init (offset, x0,dx,nx, sf_bin_int, 1, nd);
		break;
	    case 2:
		sf_int1_init (offset, x0,dx,nx, sf_lin_int, 2, nd);
		break;
	    default:
		sf_error("Unsupported interp=%d",interp);
		break;
	}


	for (it=0; it < nt; it++) { /* loop over time slices */
	    sf_floatread (dd,nd,in);
	    
	    if (gauss) {
		sf_conjgrad(NULL, sf_int1_lop, sf_freqfilt_lop, 
			    pp, mm, dd, niter);
	    } else {
		sf_conjgrad(NULL, sf_int1_lop, sf_triangle1_lop, 
			    pp, mm, dd, niter);
	    }

	    if (pef) {
		monofshape_set(0.1,nx,mm,100);
		sf_conjgrad(NULL, sf_int1_lop, monofshape_lop, pp, mm, dd, niter);
	    }
	
	    sf_floatwrite (mm,nx,out);
	}
    }
    if (!verb) sf_warning(".");

    exit(0);
}

/* 	$Id: Mshapebin1.c 7203 2011-05-03 23:49:03Z sfomel $	 */
