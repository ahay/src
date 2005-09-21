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
static float **t, **d, **t2=NULL;
static void (*transform)(bool);

static void predict_forw(bool adj, float *tt, int i, int j)
/* Predict forward */
{
    int i2;

    for (i2=i; i2 < i+j; i2++) {
	predict_step(adj,true,tt,d[i2]);
    }
}

static void predict_back(bool adj, float *tt, int i, int j)    
/* Predict backward */
{
    int i2;

    for (i2=i+j-1; i2 >= i; i2--) {
	predict_step(adj,false,tt,d[i2]);
    }
}

static void haar(bool adj) 
/* Lifting Haar transform in place */
{
    int i, j, i1;

    if (adj) {
	for (j=1; j <= nt/2; j *= 2) {
	    for (i=0; i < nt-j; i += 2*j) {
		if (inv) {
		    predict_forw(false,t[i],i,j);
		    for (i1=0; i1 < n; i1++) {
			t[i+j][i1] -= t[i][i1];     /* d = o - P[e] */
			t[i][i1]   += t[i+j][i1]/2; /* s = e + U[d] */
		    }
		    predict_back(false,t[i],i,j);
		} else {
		    predict_forw(true,t[i],i,j);
		    for (i1=0; i1 < n; i1++) {
			t[i][i1]   += t[i+j][i1];
			t[i+j][i1] -= t[i][i1]/2;
		    }
		    predict_back(true,t[i],i,j);
		}
	    }
	}
    } else {
	for (j=nt/2; j >= 1; j /= 2) {
	    for (i=0; i < nt-j; i += 2*j) {
		predict_forw(false,t[i],i,j);
		for (i1=0; i1 < n; i1++) {
		    t[i][i1]   -= t[i+j][i1]/2;
		    t[i+j][i1] += t[i][i1];
		}
		predict_back(false,t[i],i,j);
	    }
	}
    } 
}

static void linear(bool adj)
{
    int i, j, i1;
    float **t1;

    t1=t;
    for (i=0; i < nt; i++) {
	for (i1=0; i1 < n; i1++) {
	    t2[i][i1] = t[i][i1];
	}
    }
    
    if (adj) {
	for (j=1; j <= nt/2; j *= 2) {
	    if (inv) {
		for (i=0; i < nt-2*j; i += 2*j) {
		    predict_forw(false,t1[i],i,j);
		    predict_back(false,t2[i+2*j],i+j,j);
		    for (i1=0; i1 < n; i1++) {
			t[i+j][i1] -= (t1[i][i1]+t2[i+2*j][i1])/2;
			t1[i][i1] += t[i+j][i1]/2;
			t2[i][i1] += t[i+j][i1]/2;
		    }
		    predict_back(false,t1[i],i,j);
		    predict_forw(false,t2[i+2*j],i+j,j);
		}
		if (i+j < nt) {
		    predict_forw(false,t1[i],i,j);
		    for (i1=0; i1 < n; i1++) {
			t[i+j][i1] -= t1[i][i1];
			t1[i][i1] += t[i+j][i1]/2;
		    }		    
		    predict_back(false,t1[i],i,j);
		}
		for (i1=0; i1 < n; i1++) {
		    t[0][i1] = t2[j][i1];
		}
		for (i=2*j; i < nt-j; i += 2*j) {
		    for (i1=0; i1 < n; i1++) {
			t[i][i1] = (t1[i-j][i1]+t2[i+j][i1])/2;
		    }
		}
	    }
	}
    } else {	
    }
}

void seislet_init(int n1      /* trace length */, 
		  int n2      /* number of traces */, 
		  bool inv1   /* inversion flag */, 
		  float eps   /* regularization parameter */,
		  float **dip /* local slope */,
		  char type   /* transform type */) 
/*< allocate space >*/
{
    inv = inv1;
    n = n1;
    for (nt=1; nt < n2; nt *= 2) ;
    t = sf_floatalloc2(n,nt);
    d = dip;
    predict_init (n, nt, eps*eps, 1);

    switch(type) {
	case 'h': 
	    transform = haar;
	    break;
	case 'l':
	    transform = linear;

	    t2 = sf_floatalloc2(n,nt);
	    break;
	default:
	    sf_error("Unknown wavelet type=%c",type);
	    break;
    }
}

void seislet_close(void) 
/*< deallocate space >*/
{
    free (*t);
    free (t);
    if (NULL != t2) {
	free (*t2);
	free (t2);
    }
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

    transform(adj);

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
