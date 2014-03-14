/* Smoothing by spraying */
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

#include "spray.h"
#include "smspray.h"

static int ns2, nu;
static float **u, *w, *w1;

void smspray_init(int n1    /* data size */, 
		  int ns    /* spray radius */,
		  char type /* weight type */)
/*< initialize >*/
{
    int is, i1;
    float a, *t;

    ns2 = spray_init(ns);
    nu = ns2*n1;
    u = sf_floatalloc2(ns2,n1);
    w = sf_floatalloc(ns2);
    w1 = sf_floatalloc(n1);
    a = 3.0f/(ns*(ns+2));

    for (is=0; is < ns2; is++) {
	switch(type) {
	    case 't':
		w[is]=ns+1-SF_ABS(is-ns);
		break;
	    case 'g':
		w[is]=expf(-a*(is-ns)*(is-ns));
		break;
	}
    }

    /* Normalization */
    t = sf_floatalloc(n1);

    for (i1=0; i1 < n1; i1++) {
	w1[i1]=1.0f;
    }

    smspray_lop(false,false,n1,n1,w1,t);

    for (i1=0; i1 < n1; i1++) {
	if (0.0f != t[i1]) {
	    w1[i1]=1.0/t[i1];
	} else {
	    w1[i1]=0.0f;
	}
    }

    free(t);    
}

void smspray_close(void)
/*< free allocated storage >*/
{
    free(u);
    free(w);
    free(w1);
}

void smspray_lop(bool adj, bool add, int n1, int n2, float* trace, float *smooth)
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

	spray_lop(true, add,    n1, nu, trace, u[0]);
    } else {
	spray_lop(false, false, n1, nu, trace, u[0]);

	for (i1=0; i1 < n1; i1++) {
	    ws=w1[i1]; 
	    for (is=0; is < ns2; is++) {
		smooth[i1] += u[i1][is]*w[is]*ws;
	    }
	}
    }
}
