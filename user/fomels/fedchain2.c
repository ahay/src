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

static int m1, m2, n, nc;
static float *x0, *xn, *w, *rc, *t;

void fedchain2_init(int n1  /* data size */,
		    int n2,
		    int nc1 /* number of chain elements */,
		    float *x1  /* [n] input */,
		    float *w1 /* [n] weights */,
		    float *xn1 /* [n*(nc-1)] layers */)
/*< initialize >*/
{
    int ir;
    float r;

    m1 = n1;
    m2 = n2;
    n = n1*n2;
    nc = nc1;
    x0 = x1;
    xn = xn1;
    w = w1;

    t = sf_floatalloc(n);
    rc = sf_floatalloc(nc);
    for (ir=0; ir < nc; ir++) {
	r = 2*sinf(SF_PI*(ir+1)/(nc+1));
	rc[ir] = 0.5/(r*r);
    }
}

void fedchain2_close(void)
/*< free allocated storage >*/
{
    free(t);
    free(rc);
}

void fedchain2_apply(const float *x2, float *y)
/*< apply the matrix operator >*/
{
    int ir, i, j, i2, i1;
    float r, ti;
   
    for (ir=0; ir < nc; ir++) {
	r = rc[ir];

	for (i2=0; i2 < m2; i2++) {
	    for (i1=0; i1 < m1; i1++) {
		i = i2*m1+i1;
		if (ir < nc-1) {	    
		    j = ir*n+i;
		    
		    ti = 0.0f;
		    if (i1 < m1-1) {
			ti += xn[j+1]-xn[j];
		    }
		    if (i1 > 0) {
			ti += xn[j-1]-xn[j];
		    }
		    if (i2 < m2-1) {
			ti += xn[j+m1]-xn[j];
		    }
		    if (i2 > 0) {
			ti += xn[j-m1]-xn[j];
		    }
		
		    t[i] = xn[j] + r*w[i]*ti;
		} else {
		    ti = 0.0f;
		    if (i1 < m1-1) {
			ti += x0[i+1]-x0[i];
		    }
		    if (i1 > 0) {
			ti += x0[i-1]-x0[i];
		    }
		    if (i2 < m2-1) {
			ti += x0[i+m1]-x0[i];
		    }
		    if (i2 > 0) {
			ti += x0[i-m1]-x0[i];
		    }
		    
		    t[i] = x0[i] + r*w[i]*ti;
		}
	    }
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

void fedchain2(float *y)
/*< apply the chain >*/
{
    int i, ir, i1, i2;
    float r, ti;

    for (i=0; i < n; i++) {
	t[i] = x0[i];
    }

    for (ir=0; ir < nc; ir++) {
	r = rc[ir];
	
	for (i2=0; i2 < m2; i2++) {
	    for (i1=0; i1 < m1; i1++) {
		i = i2*m1+i1;

		ti = 0.0f;
		if (i1 < m1-1) {
		    ti += t[i+1]-t[i];
		}
		if (i1 > 0) {
		    ti += t[i-1]-t[i];
		}
		if (i2 < m2-1) {
		    ti += t[i+m1]-t[i];
		}
		if (i2 > 0) {
		    ti += t[i-m1]-t[i];
		}
		    
		y[i] = t[i] + r*w[i]*ti;
	    }
	}
	
        for (i=0; i < n; i++) {
	    t[i] = y[i];
	}
    }
}

void fedchain2_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
/*< linear operator >*/
{
    int ir, i, j, i1, i2;
    float r, ti;

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
		for (i2=0; i2 < m2; i2++) {
		    for (i1=0; i1 < m1; i1++) {
			i = i2*m1+i1;
			j = (ir+1)*n+i;
			
			ti = r*w[i]*t[i];

			if (i1 < m1-1) {
			    x[j+1] += ti;
			    x[j]   -= ti;
			}
			if (i1 > 0) {
			    x[j-1] += ti;
			    x[j]   -= ti; 
			}
			if (i2 < m2-1) {
			    x[j+m1] += ti;
			    x[j]    -= ti;
			}
			if (i2 > 0) {
			    x[j-m1] += ti;
			    x[j]    -= ti;
			}
			
			x[j] += t[i];
		    } 
		}	    
	    }

	    /* delta A */

	    for (i2=0; i2 < m2; i2++) {
		for (i1=0; i1 < m1; i1++) {
		    i = i2*m1+i1;
		    j = ir*n+i;
		    
		    if (ir < nc-1) {	    
			ti = 0.0f;
			if (i1 < m1-1) {
			    ti += xn[j+1]-xn[j];
			}
			if (i1 > 0) {
			    ti += xn[j-1]-xn[j];
			}
			if (i2 < m2-1) {
			    ti += xn[j+m1]-xn[j];
			}
			if (i2 > 0) {
			    ti += xn[j-m1]-xn[j];
			}
		    } else {
			ti = 0.0f;
			if (i1 < m1-1) {
			    ti += x0[i+1]-x0[i];
			}
			if (i1 > 0) {
			    ti += x0[i-1]-x0[i];
			}
			if (i2 < m2-1) {
			    ti += x0[i+m1]-x0[i];
			}
			if (i2 > 0) {
			    ti += x0[i-m1]-x0[i];
			}
		    }
		    
		    x[i] += r*y[j]*ti;
		}
	    }
	} else {

	    if (ir < nc-1) {	    
		for (i2=0; i2 < m2; i2++) {
		    for (i1=0; i1 < m1; i1++) {
			i = i2*m1+i1;
			j = (ir+1)*n+i;
		    
			ti = 0.0f;
			if (i1 < m1-1) {
			    ti += x[j+1]-x[j];
			}
			if (i1 > 0) {
			    ti += x[j-1]-x[j];
			}
			if (i2 < m2-1) {
			    ti += x[j+m1]-x[j];
			}
			if (i2 > 0) {
			    ti += x[j-m1]-x[j];
			}
		
			t[i] = x[j] + r*w[i]*ti;
		    } 
		}
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

	    for (i2=0; i2 < m2; i2++) {
		for (i1=0; i1 < m1; i1++) {
		    i = i2*m1+i1;
		    j = ir*n+i;

		    if (ir < nc-1) {	    
			ti = 0.0f;
			if (i1 < m1-1) {
			    ti += xn[j+1]-xn[j];
			}
			if (i1 > 0) {
			    ti += xn[j-1]-xn[j];
			}
			if (i2 < m2-1) {
			    ti += xn[j+m1]-xn[j];
			}
			if (i2 > 0) {
			    ti += xn[j-m1]-xn[j];
			}
		    } else {
			ti = 0.0f;
			if (i1 < m1-1) {
			    ti += x0[i+1]-x0[i];
			}
			if (i1 > 0) {
			    ti += x0[i-1]-x0[i];
			}
			if (i2 < m2-1) {
			    ti += x0[i+m1]-x0[i];
			}
			if (i2 > 0) {
			    ti += x0[i-m1]-x0[i];
			}
		    }
		    
		    y[j] += r*x[i]*ti;
		}
	    }
	}
    }
}
