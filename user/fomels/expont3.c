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

#include "expont3.h"

static int n1, n2, n3;
static float *a, *b, *c, *d, *aq, *bq, *cq, *dq;

void expont3_init(int m1, int m2, int m3 /* data dimensions */, 
		  float **aa     /* filter[4][m1*m2] ,*/
                  /* float **aaq   filter[4][m1*m3] */)
/*< initialize >*/
{
    n1 = m1;
    n2 = m2;
    n3 = m3;

    a = aa[0];
    b = aa[1];
    c = aa[2];
    d = aa[3];

    aq = aa[4];
    bq = aa[5];
    cq = aa[6];
    dq = aa[7];

}

void expont3_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy)
/*< linear operator >*/
{
    
    int i, j, k, l;
    float ab, cd, b2, d2;

    if (ny != 2*nx) sf_error("%s: size error: %d != 2*%d",__FILE__,ny,(2*nx));
    if (n1*n2*n3 != nx) sf_error("%s: size mismatch",__FILE__);    

    sf_adjnull(adj,add,nx,ny,xx,yy);

    for (l=0; l < n3; l++) {
	for (j=0; j < n2; j++) {
            for (i=4; i<n1; i++) {
	    k = i + n1*(j+n2*l);

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

    for (l=0; l < n2; l++) {
        for (j=0; j < n3; j++) {
            for (i=4; i<n1; i++) {
            k = i + n1*(j+n3*l);

            ab = -2.*aq[k]*bq[k];
            cd = -2.*cq[k]*dq[k];
            b2 = bq[k]*bq[k];
            d2 = dq[k]*dq[k];

            if (adj) {
                xx[k ] += 0.5*yy[k+nx];
                xx[k-1+n1*n2 ] += 0.5*(ab+cd)*yy[k+nx];
                xx[k-2+n1*n2] += 0.5*(b2+d2+ab*cd)*yy[k+nx];
                xx[k-3+n1*n2] += 0.5*(ab*d2+cd*b2)*yy[k+nx];
                xx[k-4+n1*n2] += 0.5*b2*d2*yy[k+nx];
            } else {
                yy[k+nx] += 0.5*(xx[k] +
                              (ab+cd)*xx[k-1+n1*n2] +
                              (b2+d2+ab*cd)*xx[k-2+n1*n2] +
                              (ab*d2+cd*b2)*xx[k-3+n1*n2] +
                              b2*d2*xx[k-4+n1*n2]);
            }
        }
    }
    }
}
