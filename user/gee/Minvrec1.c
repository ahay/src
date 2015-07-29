/* 1-D inverse interpolation with recursive filtering. */
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
    bool spline, verb, movie;
    int id, nd, nt, it, nx, ia, na, niter, nw;
    float *mm, *dd, **mmov, *offset, x0, dx, xmin, xmax, f, eps, a0;
    char *header, *lagfile;
    sf_filter aa;
    sf_file in, out, head, filt, lag;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nd)) nd=1;
    if (!sf_histint(in,"n2",&nt)) nt=1;
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

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
    
    if (!sf_getbool("spline",&spline)) spline=false;
    /* if use spline interpolation */

    filt = sf_input("filt");
    if (!sf_histint(filt,"n1",&na)) sf_error("No n1= in filt");
    aa = sf_allocatehelix (na);

    if (!sf_histfloat(filt,"a0",&a0)) a0=1.;
 
    if (NULL == (lagfile = sf_histstring(filt,"lag"))) {
	if (NULL == (lagfile = sf_getstring("lag"))) {
	    /* optional input file with filter lags */
	    for (ia=0; ia < na; ia++) {
		aa->lag[ia]=ia+1;
	    }
	    lag = NULL;
	} else {
	    lag = sf_input("lag");
	}
    } else {
	lag = sf_input(lagfile); 
    }

    if (NULL != lag) {
	if (SF_INT != sf_gettype(lag)) 
	    sf_error("Need int data in lag file '%s'",lagfile);
	sf_intread(aa->lag,na,lag);
	sf_fileclose(lag);
    }

    sf_floatread(aa->flt,na,filt);
    sf_fileclose(filt);
    
    for (ia=0; ia < na; ia++) {
	aa->flt[ia] /= a0;
    }

    if (!sf_getint("nw",&nw)) nw=2;
    /* interpolator size */

    /* initialize interpolation */
    sf_int1_init (offset,x0,dx,nx,spline? sf_spline_int: sf_lg_int,nw,nd);

    /* initialize recursive filtering */
    sf_polydiv_init (nx, aa);

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    mm = sf_floatalloc(nx);
    dd = sf_floatalloc(nd);

    if (!sf_getint("niter",&niter)) niter=nx;
    /* number of conjugate-gradient iterations */
    if (!sf_getfloat("eps",&eps)) eps=0.2;
    /* regularization parameter */
 
    if (!sf_getbool("movie",&movie)) movie=false;
    /* verbosity flag */

    if (movie) {
	mmov = sf_floatalloc2(nx,niter);
	sf_putint(out,"n2",niter);
	sf_putfloat(out,"d2",1);
	sf_putfloat(out,"o2",1);
	sf_shiftdim(in, out, 2);
    } else {
	mmov = NULL;
    }

    for (it=0; it < nt; it++) { /* loop over time slices */
	sf_floatread (dd,nd,in);
	sf_cdstep_init();
	sf_solver_prec(sf_int1_lop, sf_cdstep, sf_polydiv_lop, nx, nx, 
		       nd, mm, dd, niter, eps, "verb", verb, "xmov",mmov,
		       "end");
	sf_cdstep_close();

	if (movie) {
	    sf_floatwrite (mmov[0],nx*niter,out);
	} else {
	    sf_floatwrite (mm,nx,out);
	}
    }


    exit(0);
}

/* 	$Id: Minvbin1.c 2521 2007-02-02 00:25:42Z sfomel $	 */
