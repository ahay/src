/* Derivative for two-frequency estimation */
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
/*^*/

#include "expder2.h"

static int n1, n2;
static float *a, *b, *c, *d;

void expder2_init(int m1, int m2 /* data size */, 
		  float **aa     /* filter [4][m1*m2] */)
/*< initialize >*/
{
    n1 = m1;
    n2 = m2;

    a = aa[0];
    b = aa[1];
    c = aa[2];
    d = aa[3];
}

void expder2_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy)
/*< linear operator >*/
{
    
    int i, j, k;

    if (ny != 4*nx) sf_error("%s: size error: %d != 4*%d",__FILE__,ny,nx);

    sf_adjnull(adj,add,nx,ny,xx,yy);

    for (j=0; j < n2; j++) {
	for (i=4; i < n1; i++) {
	    k = i + j*n1;
	
	    if (adj) {
		xx[k-1] += b[k]*yy[k] + a[k]*yy[k+nx] + 
		    d[k]*yy[k+2*nx] + c[k]*yy[k+3*nx];
		xx[k-2] -= (2.*b[k]*d[k]*(c[k]*yy[k] + a[k]*yy[k+2*nx]) +
			    (b[k]+2.*a[k]*c[k]*d[k])*yy[k+nx] +
			    (d[k]+2.*a[k]*c[k]*b[k])*yy[k+3*nx]);
		xx[k-3] += b[k]*d[k]*(d[k]*yy[k]+b[k]*yy[k+2*nx]) +
		    d[k]*(a[k]*d[k]+2.*b[k]*c[k])*yy[k+nx] +
		    b[k]*(c[k]*b[k]+2.*d[k]*a[k])*yy[k+3*nx];
		xx[k-4] -= b[k]*d[k]*(d[k]*yy[k+nx] + b[k]*yy[k+3*nx]);
	    } else {
		yy[k] += b[k]*(xx[k-1] + 
			       d[k]*(d[k]*xx[k-3] - 2.*c[k]*xx[k-2])); 
		yy[k+nx] += a[k]*xx[k-1]-(b[k]+2.*a[k]*c[k]*d[k])*xx[k-2] +
		    d[k]*((a[k]*d[k]+2.*b[k]*c[k])*xx[k-3] - 
			  d[k]*b[k]*xx[k-4]);
		yy[k+2*nx] += d[k]*(xx[k-1] + 
				    b[k]*(b[k]*xx[k-3] - 2.*a[k]*xx[k-2]));
		yy[k+3*nx] += c[k]*xx[k-1]-(d[k]+2.*a[k]*c[k]*b[k])*xx[k-2] +
		    b[k]*((c[k]*b[k]+2.*d[k]*a[k])*xx[k-3] - 
			  d[k]*b[k]*xx[k-4]);
	    }
	}
    }
}
