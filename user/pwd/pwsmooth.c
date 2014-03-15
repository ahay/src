/* Structure-oriented smoothing */
/*
  Copyright (C) 2014 University of Texas at Austin
  
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

#include "pwspray.h"
#include "pwsmooth.h"

static int ns2, nu;
static float **u, *w, *w1;

void pwsmooth_init(int ns      /* spray radius */,
		   int n1      /* trace length */,
		   int n2      /* number of traces */,
		   int order   /* PWD order */,
		   float eps   /* regularization */,
		   float **dip /* local slope */)
/*< initialize >*/
{
    int is, i1, n12;
    float *t;

    ns2 = pwspray_init(ns,n1,n2,order,eps,dip);
    n12 = n1*n2;
    nu = ns2*n12;
    u = sf_floatalloc2(ns2,n12);
    w = sf_floatalloc(ns2);
    w1 = sf_floatalloc(n12);

    for (is=0; is < ns2; is++) {
	w[is]=ns+1-SF_ABS(is-ns);
    }

    /* Normalization */
    t = sf_floatalloc(n12);

    for (i1=0; i1 < n12; i1++) {
	w1[i1]=1.0f;
    }

    pwsmooth_lop(false,false,n12,n12,w1,t);

    for (i1=0; i1 < n12; i1++) {
	if (0.0f != t[i1]) {
	    w1[i1]=1.0/t[i1];
	} else {
	    w1[i1]=0.0f;
	}
    }

    free(t);    
}

void pwsmooth_close(void)
/*< free allocated storage >*/
{
    free(*u);
    free(u);
    free(w);
    free(w1);
}

void pwsmooth_lop(bool adj, bool add, 
		  int n1, int n2, float* trace, float *smooth)
/*< linear operator >*/
{
    int i1, is;
    float ws;

    if (n1 != n2) sf_error("%s: wrong size %d != %d",__FILE__,n1,n2);

    sf_adjnull(adj,add,n1,n2,trace,smooth);

    if (adj) {
	for (i1=0; i1 < n1; i1++) {
	    ws=w1[i1]; 
	    for (is=0; is < ns2; is++) {
		u[i1][is] = smooth[i1]*w[is]*ws;
	    }
	}

	pwspray_lop(true, add,    n1, nu, trace, u[0]);
    } else {
	pwspray_lop(false, false, n1, nu, trace, u[0]);

	for (i1=0; i1 < n1; i1++) {
	    ws=w1[i1]; 
	    for (is=0; is < ns2; is++) {
		smooth[i1] += u[i1][is]*w[is]*ws;
	    }
	}
    }
}
