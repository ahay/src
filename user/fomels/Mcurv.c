/* 2-D curvature estimation by plane wave destruction. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

#include "dip3.h"
#include "mask6.h"

void reverse(float *input, float *output, int *n); 
int main (int argc, char *argv[])
{
    int n123, niter, order, dorder, nj1,nj2, i,j, liter, dim;
    int n[SF_MAX_DIM], rect[3], n4, nr, ir, i1, i2, i3; 
    float p0, q0, *u, *p, *q, *m, *rp, *ru, *rm, *inp, *outp, *der;
    float pmin, pmax, qmin, qmax, d1, d2, d3;
    char key[4];
    bool verb, *m1, *m2, dip;
    sf_file in, out, mask, idip0, xdip0;

    sf_init(argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;
    if (!sf_histfloat(in,"d3",&d3)) d3=1.;

    dim = sf_filedims(in,n);
    if (dim < 2) n[1]=1;
    if (dim < 3) n[2]=1;
    n123 = n[0]*n[1]*n[2];
    nr = 1;
    for (j=3; j < dim; j++) {
	nr *= n[j];
    }

    if (1 == n[2]) {
	n4=0; 
    } else {
	if (!sf_getint("n4",&n4)) n4=2;
	/* what to compute in 3-D. 0: in-line, 1: cross-line, 2: both */ 
	if (n4 > 2) n4=2;
	if (2==n4) {
	    sf_putint(out,"n4",n4);
	    for (j=3; j < dim; j++) {
		snprintf(key,4,"n%d",j+2);
		sf_putint(out,key,n[j]);
	    }
	}
    }

    if (!sf_getint("niter",&niter)) niter=5;
    /* number of iterations */
    if (!sf_getint("liter",&liter)) liter=20;
    /* number of linear iterations */

    if (!sf_getint("rect1",&rect[0])) rect[0]=1;
    /* dip smoothness on 1st axis */
    if (!sf_getint("rect2",&rect[1])) rect[1]=1;
    /* dip smoothness on 2nd axis */
    if (!sf_getint("rect3",&rect[2])) rect[2]=1;
    /* dip smoothness on 3rd axuis */

    if (!sf_getfloat("p0",&p0)) p0=0.;
    /* initial in-line dip */
    if (!sf_getfloat("q0",&q0)) q0=0.;
    /* initial cross-line dip */

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */
    if (!sf_getint("dorder",&dorder)) dorder=6;
    /* Derivative filter order */
    if (!sf_getint("nj1",&nj1)) nj1=1;
    /* in-line antialiasing */
    if (!sf_getint("nj2",&nj2)) nj2=1;
    /* cross-line antialiasing */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    if (!sf_getbool("dip",&dip)) dip = false;
    /* If y, output dip, otherwise output curvature */
    if (!sf_getfloat("pmin",&pmin)) pmin = -FLT_MAX;
    /* minimum inline dip */
    if (!sf_getfloat("pmax",&pmax)) pmax = +FLT_MAX;
    /* maximum inline dip */
    if (!sf_getfloat("qmin",&qmin)) qmin = -FLT_MAX;
    /* minimum cross-line dip */
    if (!sf_getfloat("qmax",&qmax)) qmax = +FLT_MAX;
    /* maximum cross-line dip */

    /* initialize dip estimation */
    dip3_init(n[0], n[1], n[2], rect, liter, true);

    u = sf_floatalloc(n123);
    p = sf_floatalloc(n123);
    m = sf_floatalloc(n123);
    inp = sf_floatalloc(n[0]);
    outp = sf_floatalloc(n[0]);
    der = sf_floatalloc(n123);
    ru = sf_floatalloc(n123);
    rp = sf_floatalloc(n123);
    rm = sf_floatalloc(n123);
    q = sf_floatalloc(n123);

    if (NULL != sf_getstring("mask")) {
	m1 = sf_boolalloc(n123);
	m2 = sf_boolalloc(n123);
	mask = sf_input("mask");
    } else {
	m1 = NULL;
	m2 = NULL;
	mask = NULL;
    }

    if (NULL != sf_getstring("idip")) {
	/* initial in-line dip */
	idip0 = sf_input("idip");
    } else {
	idip0 = NULL;
    }

    if (NULL != sf_getstring("xdip")) {
	/* initial cross-line dip */
	xdip0 = sf_input("xdip");
    } else {
	xdip0 = NULL;
    }

    sf_deriv_init(n[0], dorder, 0.);

    for (ir=0; ir < nr; ir++) {
    	if (NULL != mask) {
	    sf_floatread(m,n123,mask);
	    reverse(m,rm,n);
	    mask32 (order, nj1, nj2, n[0], n[1], n[2], m, m1, m2);
	}

	/* read data */
	sf_floatread(u,n123,in);
	reverse(u,ru,n);
	
	if (1 != n4) {
	    /* initialize t-x dip */
	    if (NULL != idip0) {
		sf_floatread(p,n123,idip0);
		reverse(p,rp,n);
	    } else {
		for(i=0; i < n123; i++) {
		    p[i] = p0;
		    rp[i] = p0;
		}
	    }
	    
	    /* estimate t-x dip */
	    dip3(1, niter, order, nj1, verb, u, p, m1, pmin, pmax);

	    /* reverse dip estimation */
	    if (NULL != mask) {
		mask32 (order, nj1, nj2, n[0], n[1], n[2], rm, m1, m2);
	    }
	    dip3(1, niter, order, nj1, verb, ru, rp, m1, pmin, pmax);
	    reverse(rp,q,n);

	    if (dip) {
		for (i=0; i < n123; i++) {
		    p[i] = 0.5*(p[i]-q[i])*d1/(d2+FLT_EPSILON);
		}
	    } else {
		for (i=0; i < n123; i++) {
		    der[i] = 0.5*(p[i]-q[i]);
		}
		for (i3=0; i3 < n[2]; i3++) {
		    for (i2=0; i2 < n[1]; i2++) {
			for (i1=0; i1 < n[0]; i1++) {
			    inp[i1] = der[i3*n[1]*n[0]+i2*n[0]+i1];
			}
			sf_deriv(inp,outp);
			for (i1=0; i1 < n[0]; i1++) {
			    der[i3*n[1]*n[0]+i2*n[0]+i1] *= outp[i1];
			    p[i3*n[1]*n[0]+i2*n[0]+i1] = 
				(p[i3*n[1]*n[0]+i2*n[0]+i1] +
				 q[i3*n[1]*n[0]+i2*n[0]+i1] +
				 der[i3*n[1]*n[0]+i2*n[0]+i1]) *
				d1/(d2*d2+FLT_EPSILON);
			}
		    }
		}
	    }
	    /* write t-x dip */
	    sf_floatwrite(p,n123,out);
	}

	if (0 == n4) continue; /* done if only t-x dip */

	/* initialize t-y dip */
	if (NULL != xdip0) {
	    sf_floatread(p,n123,xdip0);
	    reverse(p,rp,n);
	} else {
	    for(i=0; i < n123; i++) {
		p[i] = q0;
		rp[i] = q0;
	    }
	}	
	if (NULL != mask) {
	    mask32 (order, nj1, nj2, n[0], n[1], n[2], m, m1, m2);
	}
	
	/* estimate t-y dip */
	dip3(2, niter, order, nj2, verb, u, p, m2, pmin, pmax);

	/* reverse dip estimation */
	if (NULL != mask) {
	    mask32 (order, nj1, nj2, n[0], n[1], n[2], rm, m1, m2);
	}
	dip3(1, niter, order, nj1, verb, ru, rp, m1, pmin, pmax);
	reverse(rp,q,n);

	if (!dip) {
	    for (i=0; i < n123; i++) {
		p[i] = 0.5*(p[i]-q[i])*d1/(d3+FLT_EPSILON);
	    }
	} else {
	    for (i=0; i < n123; i++) {
		der[i] = 0.5*(p[i]-q[i]);
	    } 
	    for (i3=0; i3 < n[2]; i3++) {
		for (i2=0; i2 < n[1]; i2++) {
		    for (i1=0; i1 < n[0]; i1++) {
			inp[i1] = der[i3*n[1]*n[0]+i2*n[0]+i1];
		    }
		    sf_deriv(inp,outp);
		    for (i1=0; i1 < n[0]; i1++) {
			der[i3*n[1]*n[0]+i2*n[0]+i1] *= outp[i1];
			p[i3*n[1]*n[0]+i2*n[0]+i1] = 
			    (p[i3*n[1]*n[0]+i2*n[0]+i1] +
			     q[i3*n[1]*n[0]+i2*n[0]+i1] +
			     der[i3*n[1]*n[0]+i2*n[0]+i1]) *
			    d1/(d3*d3+FLT_EPSILON);
		    }
		}
	    }
	}	
	/* write t-y dip */
	sf_floatwrite(p,n123,out);
    }

    exit (0);
}

void reverse(float *input, float *output, int *n) 
{
    int i1, i2, i3;

    for (i3=0; i3 < n[2]; i3++) {
	for (i2=0; i2 < n[1]; i2++) {
	    for (i1=0; i1 < n[0]; i1++) {
		output[(n[2]-1-i3)*n[1]*n[0]+(n[1]-1-i2)*n[0]+i1] =
		    input[i3*n[1]*n[0]+i2*n[0]+i1];
	    }
	}
    }
}


/* 	$Id$	 */
