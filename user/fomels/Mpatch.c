/* Patching (N-dimensional). 

w is window size (defaults to n1,n2,...)
p is number of patches in different dimensions (defaults to 1,1,...)

If inv=n, the number of output dimensions is twice the number of input dimensions.
If inv=y, the number of output dimensions is half the number of input dimensions.

September 2013 program of the month:
http://ahay.org/rsflog/index.php?/archives/357-Program-of-the-month-sfpatch.html
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
#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>

#include <rsf.h>

#include "ocpatch.h"
#include "oc.h"

int main (int argc, char *argv[])
{
    bool inv, verb, weight;
    int dim, dim0, w12,p12, j, ip, iw;
    int n[SF_MAX_DIM], p[SF_MAX_DIM], w[SF_MAX_DIM]; 
    off_t nall, n12;
    float *u, *t, *r;
    char key[4], *tmpname, *whtname;
    FILE *tmp=NULL, *wht=NULL;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

    if (!sf_getbool("inv",&inv)) inv=false;
    /* inverse or forward operation */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (!sf_getbool("weight",&weight)) weight=false;
    /* if y, apply weighting to each patch */

    if (inv) {
	dim0 = sf_filedims(in,w)/2;

	if (!sf_getint("dim",&dim)) dim=dim0; 
	for (j=0; j < dim; j++) {
	    p[j]=w[dim+j];
	}
	if (!sf_getints("n0",n,dim) && !sf_histints(in,"n0",n,dim))
	    sf_error("Need n0=");
	/* data dimensions (for inv=y) */

	for (j=0; j < dim; j++) {
	    snprintf(key,4,"n%d",j+1);
	    sf_putint(out,key,n[j]);
	    snprintf(key,4,"n%d",dim+j+1);
	    sf_putint(out,key,1);
	}	
    } else {
	dim = sf_filedims(in,n);
	sf_putints(out,"n0",n,dim);

	if (!sf_getints("w",w,dim)) {
	    /* window size */
	    for (j=0; j < dim; j++) {
		w[j] = n[j];
	    }
	} else {
	    for (j=0; j < dim; j++) {
		if (w[j] > n[j]) w[j] = n[j];
	    }
	}

	for (j=0; j < dim; j++) {
	    snprintf(key,4,"n%d",j+1);
	    sf_putint(out,key,w[j]);
	}

	if (!sf_getints("p",p,dim)) {
	    /* number of windows */
	    for (j=0; j < dim; j++) {
		if (n[j] > w[j]) {
		    p[j] = 1 + 1.5*n[j]/w[j]; /* 50% overlap */
		} else {
		    p[j] = 1;
		}
	    }	
	}

	for (j=0; j < dim; j++) {
	    snprintf(key,4,"n%d",dim+j+1);
	    sf_putint(out,key,p[j]);
	}
    }
	
    n12 = w12 = p12 = 1;    
    for (j=0; j < dim; j++) {
	if (verb) 
	    sf_warning("n[%d]=%d\tw[%d]=%d\tp[%d]=%d",
		       j,n[j],j,w[j],j,p[j]);
	    
	n12 *= n[j];
	w12 *= w[j];
	p12 *= p[j];
    }    

    u = sf_floatalloc(w12);
    nall = n12*sizeof(float);
    
    if (inv) {
	tmp = sf_tempfile(&tmpname,"w+b");
	oc_zero(nall,tmp);

	t = sf_floatalloc(w12);
    } else {
	sf_unpipe(in,nall);
	t = NULL;
    }

    if (weight) {
	r = sf_floatalloc(w12);
	sf_tent2 (dim,w,r);
	
	if (inv) {
	    wht = sf_tempfile(&whtname,"w+b");
	    oc_zero(nall,wht);
	}
    } else {
	r = NULL;
    }
	

    ocpatch_init(dim,w12,p12,p,n,w);
    
    for (ip=0; ip < p12; ip++) {	
	if (verb) sf_warning("patch %d of %d;",ip+1,p12);

	if (inv) {
	    sf_floatread (u,w12,in);
	    ocpatch_lop (ip,false,tmp,t);

	    if (weight) {
		for (iw=0; iw < w12; iw++) {
		    t[iw] += r[iw]*u[iw];
		}
	    } else {
		for (iw=0; iw < w12; iw++) {
		    t[iw] += u[iw];
		}
	    }

	    ocpatch_lop (ip,true,tmp,t);

	    if (weight) {
		ocpatch_lop (ip,false,wht,t);
		for (iw=0; iw < w12; iw++) {
		    t[iw] += r[iw];
		}
		ocpatch_lop (ip,true,wht,t);
	    }

	} else {

	    ocpatch_flop (ip,false,in,u);
	    if (NULL != r) {
		for (iw=0; iw < w12; iw++) {
		    u[iw] *= r[iw];
		}
	    }
	    sf_floatwrite (u,w12,out);

	}
    }	
    if (verb) sf_warning(".");
    
    if (inv) {
	if (weight) {
	    oc_divide(nall,tmp,wht,out);
	    unlink(whtname);
	} else {
	    oc_dump(nall,tmp,out);
	}
	unlink(tmpname);
    }

    exit (0);
}


