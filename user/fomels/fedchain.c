/* Fast explicit diffusion as a chain of operators */
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

static int n, nc;
static float *x0, *xn, *w;

void fedchain_init(int n1  /* data size */, 
		   int nc1 /* number of chain elements */,
		   float *x1  /* [n] input */,
		   float *w1 /* [n] weights */,
		   float *xn1 /* [n*(nc-1)] layers */)
/*< initialize >*/
{
    n = n1;
    nc = nc1;
    x0 = x1;
    xn = xn1;
    w = w1;
}

void fedchain_apply(const float *x2, float *y)
/*< apply the matrix operator >*/
{
    int ir, i, j;
    float r;
   
    for (ir=0; ir < nc; ir++) {
	r = 2*sinf(SF_PI*(ir+1)/nc);
	r = 1/(r*r);
	
	if (ir < nc-1) {	    
	    j = ir*n;	    
	    y[j] = xn[j]+ r*(xn[j+1]-xn[j])*w[0];
	    for (i=1; i < n-1; i++) {
		j = ir*n+i;
		y[j] = xn[j]+ r*(xn[j-1] + xn[j+1] - 2*xn[j])*w[i];
	    }
	    j = ir*n+n-1;
	    y[j] = xn[j]+ r*(xn[j-1]-xn[j])*w[n-1];
	} else {
	    j = ir*n;	    
	    y[j] = x0[0]+ r*(x0[1]-x0[0])*w[0];
	    for (i=1; i < n-1; i++) {
		j = ir*n+i;
		y[j] = x0[i]+ r*(x0[i-1] + x0[i+1] - 2*x0[i])*w[i];
	    }
	    j = ir*n+n-1;
	    y[j] = x0[n-1]+ r*(x0[n-2]-x0[n-1])*w[n-1];
		
	}
	
	for (i=0; i < n; i++) {
	    j = ir*n+i;
	    if (0==ir) {
		y[j] = x2[i] - y[j];
	    } else {
		y[j] = xn[j-n] - y[j];
	    }
	}    
    }
}

void fedchain(float *y)
/*< apply the chain >*/
{
    int i, ir;
    float *tmp, r;

    tmp = sf_floatalloc(n);

    for (i=0; i < n; i++) {
	tmp[i] = x0[i];
    }

    for (ir=0; ir < nc; ir++) {
	r = 2*sinf(SF_PI*(ir+1)/nc);
	r = 1/(r*r);

	y[0] = tmp[0]+ r*(tmp[1]-tmp[0])*w[0];
	for (i=1; i < n-1; i++) {
	    y[i] = tmp[i]+ r*(tmp[i-1] + tmp[i+1] - 2*tmp[i])*w[i];
	}
	y[n-1] = tmp[n-1]+ r*(tmp[n-2]-tmp[n-1])*w[n-1];

        for (i=0; i < n; i++) {
	    tmp[i] = y[i];
	}
    }
    
    free(tmp);
}

#ifdef aedfadrgwfgw

void fedchain_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
/*< linear operator >*/
{
    int ic, i, j, k;
    float r;

    if (nx != nc*n || ny != nc*n) sf_error("%s: Wrong size",__FILE__);

    sf_adjnull(adj,add,nx,ny,x,y);

    if (adj) {
	/* delta A */
    } else {
	/* delta A */

    }
}

