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
static float *x0, *xn, *w, *rc, *t;

void fedchain_init(int n1  /* data size */, 
		   int nc1 /* number of chain elements */,
		   float *x1  /* [n] input */,
		   float *w1 /* [n] weights */,
		   float *xn1 /* [n*(nc-1)] layers */)
/*< initialize >*/
{
    int ir;
    float r;
    
    n = n1;
    nc = nc1;
    x0 = x1;
    xn = xn1;
    w = w1;

    t = sf_floatalloc(n1);
    rc = sf_floatalloc(nc);
    for (ir=0; ir < nc; ir++) {
	r = 2*sinf(SF_PI*(ir+1)/(nc+1));
	rc[ir] = 1/(r*r);
    }
}

void fedchain_close(void)
/*< free allocated storage >*/
{
    free(t);
    free(rc);
}

void fedchain_apply(const float *x2, float *y)
/*< apply the matrix operator >*/
{
    int ir, i, j;
    float r;
   
    for (ir=0; ir < nc; ir++) {
	r = rc[ir];
	
	if (ir < nc-1) {	    
	    j = ir*n;	    
	    t[0] = xn[j]+ r*(xn[j+1]-xn[j])*w[0];
	    for (i=1; i < n-1; i++) {
		j = ir*n+i;
		t[i] = xn[j]+ r*(xn[j-1] + xn[j+1] - 2*xn[j])*w[i];
	    }
	    j = ir*n+n-1;
	    t[n-1] = xn[j]+ r*(xn[j-1]-xn[j])*w[n-1];
	} else {	    
	    t[0] = x0[0]+ r*(x0[1]-x0[0])*w[0];
	    for (i=1; i < n-1; i++) {
		t[i] = x0[i]+ r*(x0[i-1] + x0[i+1] - 2*x0[i])*w[i];
	    }
	    t[n-1] = x0[n-1]+ r*(x0[n-2]-x0[n-1])*w[n-1];
		
	}
	
	for (i=0; i < n; i++) {
	    j = ir*n+i;
	    if (0==ir) {
		y[j] = x2[i] - t[i];
	    } else {
		y[j] = xn[j-n] - t[i];
	    }
	}    
    }
}

void fedchain(float *y)
/*< apply the chain >*/
{
    int i, ir;
    float r;

    for (i=0; i < n; i++) {
	t[i] = x0[i];
    }

    for (ir=0; ir < nc; ir++) {
	r = rc[ir];

	y[0] = t[0]+ r*(t[1]-t[0])*w[0];
	for (i=1; i < n-1; i++) {
	    y[i] = t[i]+ r*(t[i-1] + t[i+1] - 2*t[i])*w[i];
	}
	y[n-1] = t[n-1]+ r*(t[n-2]-t[n-1])*w[n-1];

        for (i=0; i < n; i++) {
	    t[i] = y[i];
	}
    }
}

void fedchain_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
/*< linear operator >*/
{
    int ir, i, j;
    float r;

    if (nx != nc*n || ny != nc*n) sf_error("%s: Wrong size",__FILE__);

    sf_adjnull(adj,add,nx,ny,x,y);

    for (ir=0; ir < nc; ir++) {
	r = rc[ir];

	if (adj) {
	    for (i=0; i < n; i++) {
		j = ir*n+i;
		t[i] = y[j];
		if (0 < ir) {
		    x[j] -= y[j];
		}
	    }
	    if (ir < nc-1) {	    
		j = (ir+1)*n;
		x[j] += t[0] - r*w[0]*t[0];
		x[j+1] += r*w[0]*t[0];
		for (i=1; i < n-1; i++) {
		    j = (ir+1)*n+i;
		    x[j] += t[i] - 2*r*w[i]*t[i];
		    x[j-1] += r*w[i]*t[i];
		    x[j+1] += r*w[i]*t[i];
		}
		j = (ir+1)*n+n-1;
		x[j] += t[n-1] - r*w[n-1]*t[n-1];
		x[j-1] += r*w[n-1]*t[n-1];
		} 
	    	    
	    /* delta A */
	    if (ir < nc-1) {	    
		j = ir*n;	    
		x[0]  += r*(xn[j+1]-xn[j])*y[j];
		for (i=1; i < n-1; i++) {
		    j = ir*n+i;
		    x[i] += r*(xn[j-1] + xn[j+1] - 2*xn[j])*y[j];
		}
		j = ir*n+n-1;
		x[n-1] += r*(xn[j-1]-xn[j])*y[j];
	    } else {
		j = ir*n;	    
		x[0] += r*(x0[1]-x0[0])*y[j];
		for (i=1; i < n-1; i++) {
		    j = ir*n+i;
		    x[i] += r*(x0[i-1] + x0[i+1] - 2*x0[i])*y[j];
		}
		j = ir*n+n-1;
		x[n-1] += r*(x0[n-2]-x0[n-1])*y[j];
	    }
	} else {
	    if (ir < nc-1) {	    
		j = (ir+1)*n;
		t[0] = x[j]+ r*(x[j+1]-x[j])*w[0];
		for (i=1; i < n-1; i++) {
		    j = (ir+1)*n+i;
		    t[i] = x[j]+ r*(x[j-1] + x[j+1] - 2*x[j])*w[i];
		}
		j = (ir+1)*n+n-1;
		t[n-1] = x[j]+ r*(x[j-1]-x[j])*w[n-1];
	    } else {
		 for (i=0; i < n; i++) {
		     t[i] = 0.0f;
		 }
	    }		    
	    for (i=0; i < n; i++) {
		j = ir*n+i;
		if (0==ir) {
		    y[j] += t[i];
		} else {
		    y[j] += t[i] - x[j];
		}
	    }    
	    
	    /* delta A */
	    if (ir < nc-1) {	    
		j = ir*n;	    
		y[j] += r*(xn[j+1]-xn[j])*x[0];
		for (i=1; i < n-1; i++) {
		    j = ir*n+i;
		    y[j] += r*(xn[j-1] + xn[j+1] - 2*xn[j])*x[i];
		}
		j = ir*n+n-1;
		y[j] += r*(xn[j-1]-xn[j])*x[n-1];
	    } else {
		j = ir*n;	    
		y[j] += r*(x0[1]-x0[0])*x[0];
		for (i=1; i < n-1; i++) {
		    j = ir*n+i;
		    y[j] += r*(x0[i-1] + x0[i+1] - 2*x0[i])*x[i];
		}
		j = ir*n+n-1;
		y[j] += r*(x0[n-2]-x0[n-1])*x[n-1];
	    }
	}
    }
}

void fedchainx_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
/*< linear operator for intermediate results >*/
{
    int ir, i, j;
    float r;

    if (nx != (nc-1)*n || ny != nc*n) sf_error("%s: Wrong size",__FILE__);

    sf_adjnull(adj,add,nx,ny,x,y);

    for (ir=0; ir < nc; ir++) {
	r = rc[ir];

	if (adj) {
	    for (i=0; i < n; i++) {
		j = ir*n+i;
		t[i] = y[j];
		if (0 < ir) {
		    x[j-n] -= y[j];
		}
	    }
	    if (ir < nc-1) {	    
		j = ir*n;
		x[j] += t[0] - r*w[0]*t[0];
		x[j+1] += r*w[0]*t[0];
		for (i=1; i < n-1; i++) {
		    j = ir*n+i;
		    x[j] += t[i] - 2*r*w[i]*t[i];
		    x[j-1] += r*w[i]*t[i];
		    x[j+1] += r*w[i]*t[i];
		}
		j = ir*n+n-1;
		x[j] += t[n-1] - r*w[n-1]*t[n-1];
		x[j-1] += r*w[n-1]*t[n-1];
		} 
	} else {
	    if (ir < nc-1) {	    
		j = ir*n;
		t[0] = x[j]+ r*(x[j+1]-x[j])*w[0];
		for (i=1; i < n-1; i++) {
		    j = ir*n+i;
		    t[i] = x[j]+ r*(x[j-1] + x[j+1] - 2*x[j])*w[i];
		}
		j = ir*n+n-1;
		t[n-1] = x[j]+ r*(x[j-1]-x[j])*w[n-1];
	    } else {
		 for (i=0; i < n; i++) {
		     t[i] = 0.0f;
		 }
	    }		    
	    for (i=0; i < n; i++) {
		j = ir*n+i;
		if (0==ir) {
		    y[j] += t[i];
		} else {
		    y[j] += t[i] - x[j-n];
		}
	    }    
	}
    }
}
