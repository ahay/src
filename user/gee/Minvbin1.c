/* 1-D inverse interpolation. */
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

#include "lint1.h"
#include "tcai1.h"
#include "pef.h"
#include "bound.h"

int main (int argc, char* argv[])
{
    bool prec, pef, verb;
    int id, nd, nt, it, nx, filt, niter, three=3, i;
    float *mm, *dd, *offset, *aa, x0, dx, xmin, xmax, f, eps;
    char *header;
    sf_filter bb;
    sf_file in, out, head;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&nd)) nd=1;
    nt = sf_leftsize(in,1);

    /* create coordinates */
    offset = sf_floatalloc(nd);

    header = sf_getstring("head");
    if (NULL == header) { 
	header = sf_histstring(in,"head");
	if (NULL == header) sf_error("Need head=");
    }

    head = sf_input(header);
    if (SF_FLOAT != sf_gettype(head)) sf_error("Need float head");
    sf_floatread (offset,nd,head);
    sf_fileclose (head);

    xmin = +FLT_MAX;
    xmax = -FLT_MAX;
    for (id=0; id<nd; id++) {	
	f = offset[id]; 
	if (f < xmin) xmin=f;
	if (f > xmax) xmax=f;
    }
 
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
    lint1_init (x0,dx,offset);

    if (!sf_getint("filter",&filt)) filt=1;
    /* filter type */
    if (!sf_getbool("prec",&prec)) prec=true;
    /* if y, use preconditioning */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (!prec) {
	filt++;
	sf_warning("Using differential regularization");	    
	aa = sf_floatalloc(filt);
	tcai1_init(filt,aa);
	if (filt==2) {
	    aa[0] = 1.;
	    aa[1] = -1.;
	} else {
	    aa[0] = 1.;
	    aa[1] = -2.;
	    aa[2] = 1.;
	}
    } else if (filt > 1) {
	sf_triangle1_init (filt, nx);
    }

    mm = sf_floatalloc(nx);
    dd = sf_floatalloc(nd);

    if (!sf_getint("niter",&niter)) niter=nx;
    /* number of conjugate-gradient iterations */
    if (!sf_getfloat("eps",&eps)) eps=0.2;
    /* regularization parameter */
    if (!sf_getbool("pef",&pef)) pef=false;
    /* if y, use PEF for regularization */
 
    for (it=0; it < nt; it++) { /* loop over time slices */
	sf_floatread (dd,nd,in);
	if (prec) {
	    if (filt==1) {
		sf_solver_prec(lint1_lop, sf_cgstep, sf_causint_lop, nx, nx, 
			       nd, mm, dd, niter, eps, "verb", verb, "end");
	    } else {
		sf_solver_prec(lint1_lop, sf_cgstep, sf_triangle1_lop, nx, nx, 
			       nd, mm, dd, niter, eps, "verb", verb, "end");
	    }
	} else {
	    sf_solver_reg(lint1_lop, sf_cgstep, tcai1_lop, nx+filt, nx, nd, 
			  mm, dd, niter, eps, "end");
	}
	sf_cgstep_close();

	if (pef) {
	    bb = sf_allocatehelix (2);
	    aa = sf_floatalloc (3);
	    aa[0] = 1.;
	    bb->lag[0] = 1;
	    bb->lag[1] = 2;
	    bound(1, false, &nx, &nx, &three, bb); 
	    for (i=0; i < 3; i++) {
		find_pef (nx, mm, bb, 3);
		aa[1] = bb->flt[0];
		aa[2] = bb->flt[1];
		if (prec) {
		    sf_polydiv_init (nx, bb);
		    sf_solver_prec(lint1_lop, sf_cgstep, sf_polydiv_lop, 
				   nx, nx, nd, mm, dd, niter, eps, "end");
		} else {
		    tcai1_init (3, aa);
		    sf_solver_reg(lint1_lop, sf_cgstep, tcai1_lop, nx+filt,
				  nx, nd, mm, dd, niter, eps, "end");
		}
		sf_cgstep_close();
	    }
	}
	
    	sf_floatwrite (mm,nx,out);
    }


    exit(0);
}

/* 	$Id$	 */
