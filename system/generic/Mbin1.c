/* Data binning in 1-D slices. */
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

int main (int argc, char* argv[])
{
    int id, nd, nt, it, ix, nx, interp;
    float *mm, *count, *dd, *offset, x0, dx, xmin, xmax, f, clip;
    char *header;
    sf_file in, out, head, fold, pattern;

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

    /* let user overwrite */
    sf_getfloat ("xmin",&xmin);
    /* grid dimensions */
    sf_getfloat ("xmax",&xmax);

    if (xmax <= xmin) sf_error ("xmax=%f <= xmin=%f",xmax,xmin);

    /* create model */
    if (NULL != sf_getstring("pattern")) {
	pattern = sf_input("pattern");

	if (!sf_histint (pattern,"n1",&nx)) sf_error("Need nx=");
	if (!sf_histfloat(pattern,"o1",&x0)) x0=xmin; 
	if (!sf_histfloat(pattern,"d1",&dx)) {
	    /* grid spacing */
	    if (1 >= nx) sf_error("Need dx=");
	    dx = (xmax-xmin)/(nx-1);
	}

	sf_fileclose(pattern);
    } else {
	if (!sf_getint ("nx",&nx)) sf_error("Need nx=");
	/* Number of bins */
	
	if (!sf_getfloat("x0",&x0)) x0=xmin; 
	/* grid origin */

	if (!sf_getfloat("dx",&dx)) {
	    /* grid spacing */
	    if (1 >= nx) sf_error("Need dx=");
	    dx = (xmax-xmin)/(nx-1);
	}	
    }

    sf_putint(out,"n1",nx);
    sf_putfloat (out,"o1",x0);
    sf_putfloat (out,"d1",dx);
    
    /* initialize interpolation */
    if (!sf_getint("interp",&interp)) interp=1;
    /* [1,2] interpolation method, 1: nearest neighbor, 2: linear */

    switch (interp) {
	case 1:
	    sf_int1_init (offset, x0,dx,nx, sf_bin_int, 1, nd, 0.0);
	    sf_warning("Using nearest-neighbor interpolation");
	    break;
	case 2:
	    sf_int1_init (offset, x0,dx,nx, sf_lin_int, 2, nd, 0.0);
	    sf_warning("Using linear interpolation");
	    break;
	default:
	    sf_error("Unsupported interp=%d",interp);
	    break;
    }

    mm = sf_floatalloc(nx);
    dd = sf_floatalloc(nd);
    count  = sf_floatalloc(nx);

    /* compute fold */
    for (id=0; id<nd; id++) {
	dd[id]=1.;
    }

    sf_int1_lop (true, false,nx,nd,count,dd);
 
    if (NULL != sf_getstring("fold")) {
	/* output fold file (optional) */
	fold = sf_output("fold");
	sf_putint(fold,"n1",nx);
	sf_putint(fold,"n2",1);
	sf_putfloat(fold,"o1",x0);
	sf_putfloat(fold,"d1",dx);
	sf_floatwrite (count,nx,fold);
	sf_fileclose (fold);
    }

    if (!sf_getfloat("clip",&clip)) clip = FLT_EPSILON;
    /* clip for fold normalization */

    for (ix=0; ix<nx; ix++) {
	if (clip < count[ix]) count[ix]=1./fabsf(count[ix]);
	else                  count[ix]=0.;
    }

    for (it=0; it < nt; it++) { /* loop over time slices */
	sf_floatread (dd,nd,in);
	sf_int1_lop (true,false,nx,nd,mm,dd);
	for (ix=0; ix<nx; ix++) {
	    mm[ix] *= count[ix];
	}
	sf_floatwrite (mm,nx,out);
    }


    exit(0);
}

/* 	$Id$	 */
