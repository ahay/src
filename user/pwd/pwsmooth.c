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

static int n1, n2, ns2;
static float ***u, *w, **w1;

void pwsmooth_init(int ns      /* spray radius */,
		   int m1      /* trace length */,
		   int m2      /* number of traces */,
		   int order   /* PWD order */,
		   float eps   /* regularization */,
		   float **dip /* local slope */)
/*< initialize >*/
{
    int is, i1, n12;
    float *t;

    n1 = m1;
    n2 = m2;
    n12 = n1*n2;

    ns2 = pwspray_init(ns,n1,n2,order,eps,dip);

    u = sf_floatalloc3(n1,ns2,n2);
    w = sf_floatalloc(ns2);
    w1 = sf_floatalloc2(n1,n2);

    for (is=0; is < ns2; is++) {
	w[is]=ns+1-SF_ABS(is-ns);
    }

    /* Normalization */
    t = sf_floatalloc(n12);

    for (i1=0; i1 < n12; i1++) {
	w1[0][i1]=1.0f;
    }

    pwsmooth_lop(false,false,n12,n12,w1[0],t);

    for (i1=0; i1 < n12; i1++) {
	if (0.0f != t[i1]) {
	    w1[0][i1]=1.0/t[i1];
	} else {
	    w1[0][i1]=0.0f;
	}
    }

    free(t);    
}

void pwsmooth_close(void)
/*< free allocated storage >*/
{
    free(**u);
    free(*u);
    free(u);
    free(w);
    free(*w1);
    free(w1);
}

void pwsmooth_lop(bool adj, bool add, 
		  int nin, int nout, float* trace, float *smooth)
/*< linear operator >*/
{
    int i1, i2, is;
    float ws;

    if (nin != nout || nin != n1*n2) 
	sf_error("%s: wrong size %d != %d",__FILE__,n1,n2);
    
    sf_adjnull(adj,add,nin,nout,trace,smooth);

    if (adj) {
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		ws=w1[i2][i1]; 
		for (is=0; is < ns2; is++) {
		    u[i2][is][i1] = smooth[i2*n1+i1]*w[is]*ws;
		}
	    }
	}

	pwspray_lop(true,  true,  nin, nin*ns2, trace, u[0][0]);
    } else {
	pwspray_lop(false, false, nin, nin*ns2, trace, u[0][0]);

	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		ws=w1[i2][i1]; 
		for (is=0; is < ns2; is++) {
		    smooth[i2*n1+i1] += u[i2][is][i1]*w[is]*ws;
		}
	    }
	}
    }
}
