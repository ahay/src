/* Frequency filtering, 2-D */
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

#include "expont2.h"

static int n1, n2;
static float *a, *b, *c, *d;

void expont2_init(int m1, int m2 /* data dimensions */, 
		  float **aa     /* filter[4][m1*m2] */)
/*< initialize >*/
{
    n1 = m1;
    n2 = m2;

    a = aa[0];
    b = aa[1];
    c = aa[2];
    d = aa[3];
}

void expont2_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy)
/*< linear operator >*/
{
    
    int i, j, k;
    float ab, cd, b2, d2;

    if (ny != nx) sf_error("%s: size error: %d != %d",__FILE__,ny,nx);

    sf_adjnull(adj,add,nx,ny,xx,yy);

    for (j=0; j < n2; j++) {
	for (i=4; i < n1; i++) {
	    k = i + j*n1;

	    ab = -2.*a[k]*b[k];
	    cd = -2.*c[k]*d[k];
	    b2 = b[k]*b[k];
	    d2 = d[k]*d[k];
	
	    if (adj) {
		xx[k  ] += 0.5*yy[k];
		xx[k-1] += 0.5*(ab+cd)*yy[k];
		xx[k-2] += 0.5*(b2+d2+ab*cd)*yy[k];
		xx[k-3] += 0.5*(ab*d2+cd*b2)*yy[k];
		xx[k-4] += 0.5*b2*d2*yy[k];
	    } else {
		yy[k] += 0.5*(xx[k] + 
			      (ab+cd)*xx[k-1] + 
			      (b2+d2+ab*cd)*xx[k-2] +
			      (ab*d2+cd*b2)*xx[k-3] + 
			      b2*d2*xx[k-4]);
	    }
	}
    }
}
