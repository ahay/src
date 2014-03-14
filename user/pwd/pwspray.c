/* Plane-wave spraying */
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

#include "predict.h"

static int n1, n2, ns, ns2;
static float *trace, **p;

int pwspray_init(int nr      /* spray radius */, 
		 int nt      /* trace length */, 
		 int n       /* number of traces */,
		 int order   /* PWD order */,
		 float eps   /* regularization */,
		 float **dip /* local slope */)
/*< initialize >*/
{
    n1=nt;
    n2=n;

    ns=nr;
    ns2=2*ns+1;

    predict_init (n1, n2, eps*eps, order, 1, false);
    trace = sf_floatalloc(n1);
    p = dip;

    return ns2;
}

void pwspray_close(void)
/*< free allocated storage >*/
{
    predict_close();
    free(trace);
}

void pwspray_lop(bool adj, bool add, int n, int nu, float* u1, float *u)
/*< linear operator >*/
{
    int i, is, ip, j, i1;

    if (n  != n1*n2) sf_error("%s: wrong size %d != %d*%d",__FILE__,n, n1,n2);
    if (nu != n*ns2) sf_error("%s: wrong size %d != %d*%d",__FILE__,nu,n,ns2);

    sf_adjnull(adj,add,n,nu,u1,u);

    for (i=0; i < n2; i++) { 	
	if (adj) {
	    for (i1=0; i1 < n1; i1++) {
		trace[i1] = u[(i*ns2+ns)*n1+i1];
	    }
	    for (i1=0; i1 < n1; i1++) {
		u1[i*n1+i1] += trace[i1];
	    }

	    /* predict backward */
	    for (is=ns-1; is >= 0; is--) {
		ip = i-is-1;
		if (ip < 0) continue;
		j = ip*ns2+ns-is-1;
		for (i1=0; i1 < n1; i1++) {
		    trace[i1] += u[j*n1+i1];
		}
		predict_step(false,false,trace,p[ip]);
	    }
	    

	} else {
	    for (i1=0; i1 < n1; i1++) {
		trace[i1] = u1[i*n1+i1];
	    }
	    for (i1=0; i1 < n1; i1++) {
		u[(i*ns2+ns)*n1+i1] += trace[i1];
	    }

            /* predict forward */
	    for (is=0; is < ns; is++) {
		ip = i-is-1;
		if (ip < 0) break;
		j = ip*ns2+ns-is-1;
		predict_step(false,false,trace,p[ip]);
		for (i1=0; i1 < n1; i1++) {
		    u[j*n1+i1] += trace[i1];
		}
	    }
	}
    }
}

