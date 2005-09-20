/* Seislet transform */
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

#include "seislet.h"
#include "predict.h"

static int nt, n;
static bool inv;
static float **t, **d;

static void haar(bool adj) 
/* Lifting Haar transform in place */
{
    int i, j, i1, i2;

    if (adj) {
	for (j=1; j <= nt/2; j *= 2) {
	    for (i=0; i < nt-j; i += 2*j) {
		if (inv) {
		    for (i2=i; i2 < i+j; i2++) {
			predict_step(false,true,t[i],d[i2]);
		    }
		    for (i1=0; i1 < n; i1++) {
			t[i+j][i1] -= t[i][i1];     /* d = o - P[e] */
			t[i][i1]   += t[i+j][i1]/2; /* s = e + U[d] */
		    }
		    for (i2=i+j-1; i2 >= i; i2--) {
			predict_step(false,false,t[i],d[i2]);
		    }
		} else {
		    for (i2=i; i2 < i+j; i2++) {
			predict_step(true,false,t[i],d[i2]);
		    }
		    for (i1=0; i1 < n; i1++) {
			t[i][i1]   += t[i+j][i1];
			t[i+j][i1] -= t[i][i1]/2;
		    }
		    for (i2=i+j-1; i2 >= i; i2--) {
			predict_step(true,true,t[i],d[i2]);
		    }
		}
	    }
	}
    } else {
	for (j=nt/2; j >= 1; j /= 2) {
	    for (i=0; i < nt-j; i += 2*j) {
		for (i2=i; i2 < i+j; i2++) {
		    predict_step(false,true,t[i],d[i2]);
		}
		for (i1=0; i1 < n; i1++) {
		    t[i][i1]   -= t[i+j][i1]/2;
		    t[i+j][i1] += t[i][i1];
		}
		for (i2=i+j-1; i2 >= i; i2--) {
		    predict_step(false,false,t[i],d[i2]);
		}
	    }
	}
    } 
}

void seislet_init(int n1      /* trace length */, 
		  int n2      /* number of traces */, 
		  bool inv1   /* inversion flag */, 
		  float eps   /* regularization parameter */,
		  float **dip /* local slope */) 
/*< allocate space >*/
{
    inv = inv1;
    n = n1;
    for (nt=1; nt < n2; nt *= 2) ;
    t = sf_floatalloc2(n,nt);
    d = dip;
    predict_init (n, nt, eps*eps, 1);
}

void seislet_close(void) 
/*< deallocate space >*/
{
    free (*t);
    free (t);
    predict_close();
}

void seislet_lop(bool adj, bool add, int nx, int ny, float *x, float *y)
/*< linear operator >*/
{
    int it, i, j, i1;

    sf_adjnull (adj,add,nx,ny,x,y);

    if (adj) {
	for (it=0; it < nx; it++) {
	    t[0][it] = y[it];
	}
	for (it=nx; it < n*nt; it++) {
	    t[0][it] = 0.;
	}
    } else {
	for (i1=0; i1 < n; i1++) {
	    t[0][i1] = x[i1];
	}
	it = n;
	for (j=nt/2; j >= 1; j /= 2) {
	    for (i=0; i < nt-j; i += 2*j) {
		for (i1=0; i1 < n; i1++) {
		    if (it < ny) {
			t[i+j][i1]=x[it];
			it++;
		    } else {
			t[i+j][i1]=0.;
		    }
		}
	    }	    	    
	}
    }

    haar(adj);

    if (adj) {
	for (i1=0; i1 < n; i1++) {
	    x[i1] += t[0][i1];
	}
	it = n;
	for (j=nt/2; j >= 1; j /= 2) {
	    for (i=0; i < nt-j; i += 2*j) {
		for (i1=0; i1 < n; i1++) {
		    x[it] += t[i+j][i1];
		    it++;
		    if (it >= ny) return;
		}
	    }	    	    
	}
    } else {
	for (it=0; it < nx; it++) {
	    y[it] += t[0][it];
	}
    }
}
