/* 3-D robust dip estimation by plane wave destruction with non-stationary smoothing and forward-backward space derivative calculation. 
The output is dimensionless (stepout in time measured in time samples). 
also see sfdip
*/
/*
  Copyright (C) 2021 Yangtze University and University of Texas at Austin
  
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

#include "diputil.h"

int main (int argc, char *argv[])
{
    int   *sft[SF_MAX_DIM];	/* storing non-stationary shifting size */
    float *rct[SF_MAX_DIM]; /* storing non-stationary smoothing radii */
    int n123, niter, order, nj1,nj2, i, j, liter, dim, b;
    int n[SF_MAX_DIM], n4, nr, ir; 
    float p0, q0, *u, *p, *pi=NULL, *qi=NULL;
    float pmin, pmax, qmin, qmax, eps;
    char key[8];
    int box[SF_MAX_DIM];
	/*box means maximum (edge) padding for triangle smoothing, box[i]=max(rect[i])*/
    bool verb, both, **mm, drift;
    sf_file in, out, mask, idip0, xdip0;
	sf_file rect[SF_MAX_DIM], shift[SF_MAX_DIM]; 
	
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

		/*Calculate dim1*/
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (NULL != sf_getstring(key)) {
	    /*( rect# size of the smoothing stencil in #-th dimension /auxiliary input file/ )*/
	    rect[i] = sf_input(key);
	    if (SF_FLOAT != sf_gettype(rect[i])) sf_error("Need float %s",key);
	    snprintf(key,8,"shift%d",i+1);
	    if (NULL != sf_getstring(key)) {
		/*( shift# shifting of the smoothing stencil in #-th dimension /auxiliary input file/ )*/
		shift[i] = sf_input(key);
		if (SF_INT != sf_gettype(shift[i])) sf_error("Need int %s",key);
	    } else {
		shift[i] = NULL;
	    }
	} else {
	    rect[i] = NULL;
	    shift[i] = NULL;
	}
    }
	
	/*reading the non-stationary smoothing radii*/
    for (i=0; i < dim; i++) {
	box[i] = 1;
	if (NULL != rect[i]) {
	    rct[i] = sf_floatalloc (n123);
	    sft[i] = sf_intalloc (n123);
		/* non-stationary dip smoothness on 1st/2nd/3rd axis */

	    sf_floatread(rct[i],n123,rect[i]);
	    sf_fileclose(rect[i]);

	    if (NULL != shift[i]) {
		sf_intread(sft[i],n123,shift[i]);
		sf_fileclose(shift[i]);
	    } else {
		for (j=0; j < n123; j++) {
		    sft[i][j] = 0;
		}
	    }

		
	    for (j=0; j < n123; j++) {
		b = ceilf(rct[i][j])+SF_ABS(sft[i][j]);
		if (b > box[i]) box[i] = b;
	    }	    
	} else {
	    rct[i] = NULL;
	    sft[i] = NULL;
	}
    }
    
	sf_warning("dim=%d, n123=%d",dim,n123);
	sf_warning("n[0]=%d, n[1]=%d, n[2]=%d",n[0],n[1],n[2]);	
	
	
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
    dip3n_init(n[0], n[1], n[2], box, rct, sft, liter, eps, verb);

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
	    dip3n_fb(false, 1, niter, order, nj1, drift, u, p, mm[0], pmin, pmax);
	    
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
	    dip3n_fb(false, 2, niter, order, nj2, drift, u, p, mm[1], qmin, qmax);
	    
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
	    dip3n_fb(true, 1, niter, order, nj1, drift, u, p, mm[2], -pmax, -pmin);
	    
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
	    dip3n_fb(true, 2, niter, order, nj2, drift, u, p, mm[3], -qmax, -qmin);
	    
	    /* write t-y dip */
	    sf_floatwrite(p,n123,out);
	}	
    }

    exit (0);
}

