/* 3-D dip estimation by plane wave destruction. 

The output is dimensionless (stepout in time measured in time samples). 

June 2012 program of the month:
http://ahay.org/blog/2012/06/02/program-of-the-month-sfdip/
*/
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

#include "dip3.h"
#include "mask6.h"

int main (int argc, char *argv[])
{
    int n123, niter, order, nj1,nj2, i,j, liter, dim;
    int n[SF_MAX_DIM], rect[3], n4, nr, ir; 
    float p0, q0, *u, *p, *pi=NULL, *qi=NULL;
    float pmin, pmax, qmin, qmax, eps;
    char key[4];
    bool verb, both, **mm, drift;
    sf_file in, out, mask, idip0, xdip0;

    sf_init(argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

    dim = sf_filedims(in,n);
    if (dim < 2) n[1]=1;
    if (dim < 3) n[2]=1;
    n123 = n[0]*n[1]*n[2];
    nr = 1;
    for (j=3; j < dim; j++) {
	nr *= n[j];
    }

    if (!sf_getbool("both",&both)) both=false;
    /* if y, compute both left and right predictions */

    if (1 == n[2]) {
	n4=0;
	if (both) sf_putint(out,"n3",2);
    } else {
	if (!sf_getint("n4",&n4)) n4=2;
	/* what to compute in 3-D. 0: in-line, 1: cross-line, 2: both */ 
	if (n4 > 2) n4=2;
	if (2==n4) {
	    sf_putint(out,"n4",both? 4:2);
	    for (j=3; j < dim; j++) {
		snprintf(key,4,"n%d",both? j+4:j+2);
		sf_putint(out,key,n[j]);
	    }
	} else if (both) {
	    sf_putint(out,"n4",2);
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
    if (!sf_getint("nj1",&nj1)) nj1=1;
    /* in-line antialiasing */
    if (!sf_getint("nj2",&nj2)) nj2=1;
    /* cross-line antialiasing */

    if (!sf_getbool("drift",&drift)) drift=false;
    /* if shift filter */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    if (!sf_getfloat("pmin",&pmin)) pmin = -FLT_MAX;
    /* minimum inline dip */
    if (!sf_getfloat("pmax",&pmax)) pmax = +FLT_MAX;
    /* maximum inline dip */
    if (!sf_getfloat("qmin",&qmin)) qmin = -FLT_MAX;
    /* minimum cross-line dip */
    if (!sf_getfloat("qmax",&qmax)) qmax = +FLT_MAX;
    /* maximum cross-line dip */

    if (!sf_getfloat("eps",&eps)) eps=0.0f;
    /* regularization */

    /* initialize dip estimation */
    dip3_init(n[0], n[1], n[2], rect, liter, eps, verb);

    u = sf_floatalloc(n123);
    p = sf_floatalloc(n123);

    if (NULL != sf_getstring("mask")) {
	mm = sf_boolalloc2(n123,both? 4:2);
	mask = sf_input("mask");
    } else {
	mm = (bool**) sf_alloc(4,sizeof(bool*));
	mm[0] = mm[1] = mm[2] = mm[3] = NULL;
	mask = NULL;
    }

    if (NULL != sf_getstring("idip")) {
	/* initial in-line dip */
	idip0 = sf_input("idip");
	if (both) pi = sf_floatalloc(n123);
    } else {
	idip0 = NULL;
    }

    if (NULL != sf_getstring("xdip")) {
	/* initial cross-line dip */
	xdip0 = sf_input("xdip");
	if (both) qi = sf_floatalloc(n123);
    } else {
	xdip0 = NULL;
    }

    for (ir=0; ir < nr; ir++) {
    	if (NULL != mask) {
	    sf_floatread(u,n123,mask);
	    mask32 (both, order, nj1, nj2, n[0], n[1], n[2], u, mm);
	}

	/* read data */
	sf_floatread(u,n123,in);
	
	if (1 != n4) {
	    /* initialize t-x dip */
	    if (NULL != idip0) {
		if (both) {
		    sf_floatread(pi,n123,idip0);
		    for(i=0; i < n123; i++) {
			p[i] = pi[i];
		    }
		} else {
		    sf_floatread(p,n123,idip0);
		}
	    } else {
		for(i=0; i < n123; i++) {
		    p[i] = p0;
		}
	    }
	    
	    /* estimate t-x dip */
	    dip3(false, 1, niter, order, nj1, drift, u, p, mm[0], pmin, pmax);
	    
	    /* write t-x dip */
	    sf_floatwrite(p,n123,out);
	}

	if (0 != n4) {
	    /* initialize t-y dip */
	    if (NULL != xdip0) {
		if (both) {
		    sf_floatread(qi,n123,xdip0);
		    for(i=0; i < n123; i++) {
			p[i] = qi[i];
		    }
		} else {
		    sf_floatread(p,n123,xdip0);
		}
	    } else {
		for(i=0; i < n123; i++) {
		    p[i] = q0;
		}
	    }	
	    
	    /* estimate t-y dip */
	    dip3(false, 2, niter, order, nj2, drift, u, p, mm[1], qmin, qmax);
	    
	    /* write t-y dip */
	    sf_floatwrite(p,n123,out);
	}

	if (!both) continue;

	if (1 != n4) {
	    /* initialize t-x dip */
	    if (NULL != idip0) {
		for(i=0; i < n123; i++) {
		    p[i] = -pi[i];
		}
	    } else {
		for(i=0; i < n123; i++) {
		    p[i] = -p0;
		}
	    }
	    
	    /* estimate t-x dip */
	    dip3(true, 1, niter, order, nj1, drift, u, p, mm[2], -pmax, -pmin);
	    
	    /* write t-x dip */
	    sf_floatwrite(p,n123,out);
	}

	if (0 != n4) {
	    /* initialize t-y dip */
	    if (NULL != xdip0) {
		for(i=0; i < n123; i++) {
		    p[i] = -qi[i];
		}
	    } else {
		for(i=0; i < n123; i++) {
		    p[i] = -q0;
		}
	    }	
	    
	    /* estimate t-y dip */
	    dip3(true, 2, niter, order, nj2, drift, u, p, mm[3], -qmax, -qmin);
	    
	    /* write t-y dip */
	    sf_floatwrite(p,n123,out);
	}	
    }

    exit (0);
}

/* 	$Id$	 */
