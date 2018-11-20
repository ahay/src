/* Chain of PWD operators */
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
#include "allp3.h"

static int n, nc;
static float *x0, *xn, *pn, *tmp;
static allpass *ap;

void pwdchain_init(int n1,
		   int n2     /* data size */, 
		   int nw     /* filter order */,
		   int nc1    /* number of chain elements */,		   
		   float *x1  /* [n] input */,
		   float *pn1 /* [n*nc] dips */,
		   float *xn1 /* [n*(nc-1)] layers */)
/*< initialize >*/
{
    int ic;

    n = n1*n2;
    nc = nc1;
    x0 = x1;
    xn = xn1;
    pn = pn1;

    tmp = sf_floatalloc(n);

    ap = (allpass*) sf_alloc(nc,sizeof(*ap));
    for (ic=0; ic < nc; ic++) {
	ap[ic] = allpass_init (nw,1,n1,n2,1,pn+ic*n);
    }
}

void pwdchain_close(void)
/*< free allocated storage >*/
{
    int ic;

    free(tmp);

    for (ic=0; ic < nc; ic++) {
	allpass_close(ap[ic]);
    }
    free(ap);
}

void pwdchain_apply(const float *x2, float *y)
/*< apply the matrix operator >*/
{
    int ic, i;
    float *xc, *yc;

    if (nc > 1) {
	allpass1(false, false, ap[0], xn, y);
	for (i=0; i < n; i++) {
	    y[i] = x2[i] - y[i];
	}
	for (ic=1; ic < nc-1; ic++) {
	    yc = y+ic*n;
	    xc = xn+ic*n;
	    allpass1(false, false, ap[ic], xc, yc);
	    for (i=0; i < n; i++) {
		yc[i] = xc[i-n] - yc[i];
	    }
	}
	yc = y+(nc-1)*n;
	xc = xn+(nc-1)*n;
	allpass1(false, false, ap[nc-1], x0, yc);
	for (i=0; i < n; i++) {
	    yc[i] = xc[i-n] - yc[i];
	}
    } else {
	allpass1(false, false, ap[0], x0, y);
	for (i=0; i < n; i++) {
	    y[i] = x2[i] - y[i];
	}
    }
}

void pwdchain(float *y)
/*< apply the chain >*/
{
    int i, ic;

    for (i=0; i < n; i++) {
	tmp[i] = x0[i];
    }
    for (ic=nc-1; ic >= 0; ic--) {
	allpass1(false, false, ap[ic], tmp, y);
	for (i=0; i < n; i++) {
	    tmp[i] = y[i];
	}
    }
}

void pwdchain_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
/*< linear operator >*/
{
    int ic, i, j, k;

    if (nx != (2*nc-1)*n || ny != nc*n) sf_error("%s: Wrong size",__FILE__);

    sf_adjnull(adj,add,nx,ny,x,y);

    if (adj) {
	if (nc > 1) {
	    allpass1t(false, false, ap[0], tmp, y);
	    for (i=0; i < n; i++) {
		x[i+n*nc] -= tmp[i];
	    }
	    for (ic=1; ic < nc-1; ic++) {
		allpass1t(false, false, ap[ic], tmp, y);
		for (i=0; i < n; i++) {		    
		    j = ic*n+i;
		    k = j+n*nc;
		    x[k-n] += y[j];
		    x[k] -= tmp[i];
		}
	    } 
	    for (i=0; i < n; i++) {
		j = (nc-1)*n+i;
		k = j+n*nc;
		x[k-n] += y[j];
	    } 
	} 
	/* delta A */ 
	for (ic=0; ic < nc-1; ic++) {
	    allpass1(false, true, ap[ic], xn+ic*n, tmp);
	    for (i=0; i < n; i++) {
		j = ic*n+i;
		x[j] -= y[j]*tmp[i];
	    }
	} 
	allpass1(false, true, ap[nc-1], x0, tmp);
	for (i=0; i < n; i++) {
	    j = (nc-1)*n+i;
	    x[j] -= y[j]*tmp[i];
	} 
    } else {
	if (nc > 1) {
	    allpass1(false, false, ap[0], x+nc*n, tmp);
	    for (i=0; i < n; i++) {
		y[i] -= tmp[i];
	    }
	    for (ic=1; ic < nc-1; ic++) {
		allpass1(false, false, ap[ic], x+(nc+ic)*n, tmp);
		for (i=0; i < n; i++) {
		    j = ic*n+i;
		    k = j+n*nc;
		    y[j]  += x[k-n] - tmp[i];
		}
	    } 
	    for (i=0; i < n; i++) {
		j = (nc-1)*n+i;
		k = j+n*nc;
		y[j] += x[k-n];
	    } 
	} 
	/* delta A */
	for (ic=0; ic < nc-1; ic++) {
	    allpass1(false, true, ap[ic], xn+ic*n, tmp);	    
	    for (i=0; i < n; i++) {
		j = ic*n+i;
		y[j] -= x[j]*tmp[i];
	    }
	} 
	allpass1(false, true, ap[nc-1], x0, tmp);
	for (i=0; i < n; i++) {
	    j = (nc-1)*n+i;
	    y[j] -= x[j]*tmp[i];
	} 
    }
} 
