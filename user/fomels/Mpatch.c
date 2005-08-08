/* Patching (N-dimensional). 

w is window size (defaults to n1,n2,...)
p is number of patches in different dimensions (defaults to 1,1,...)

If inv=n, the number of output dimensions is twice the number of input dimensions.
If inv=y, the number of output dimensions is half the number of input dimensions.
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
    bool inv, verb;
    int dim,w12,p12, j, ip;
    int n[SF_MAX_DIM], p[SF_MAX_DIM], w[SF_MAX_DIM]; 
    off_t nall, n12;
    float *u;
    char key[4];
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

    if (!sf_getbool("inv",&inv)) inv=false;
    /* inverse or forward operation */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    dim = sf_filedims(in,n);
    
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

    n12 = w12 = p12 = 1;    
    for (j=0; j < dim; j++) {
	if (verb) sf_warning("n[%d]=%d\tw[%d]=%d\tp[%d]=%d",j,n[j],j,w[j],j,p[j]);

	n12 *= n[j];
	w12 *= w[j];
	p12 *= p[j];

	snprintf(key,4,"n%d",j+1);
	sf_putint(out,key,w[j]);
	snprintf(key,4,"n%d",dim+j+1);
	sf_putint(out,key,p[j]);
    }

    u = sf_floatalloc(w12);
    nall = n12*sizeof(float);
    
    sf_unpipe(in,nall);
    
    ocpatch_init(dim,w12,p12,p,n,w);
    
    /* loop over patches */
    for (ip=0; ip < p12; ip++) {
	/* read data */
	ocpatch_flop (ip,false,in,u);
		
	/* write data */
	sf_floatwrite (u,w12,out);
    }		
    
    exit (0);
}

/* 	$Id: Mdip.c 1071 2005-03-20 20:08:23Z fomels $	 */
