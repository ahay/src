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

#include "expder3.h"

static int n1, n2, n3;
static float *a, *b, *c, *d, *aq, *bq, *cq, *dq;

void expder3_init(int m1, int m2, int m3 /* data size */, 
		  float **aa     /* filter [4][m1*m2] , */
                  /* float **aaq    filter [4][m1*m3] */)
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

void expder3_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy)
/*< linear operator >*/
{
    
    int i, j, k, l;

    if (ny != 8*nx) sf_error("%s: size error: %d != 8*%d",__FILE__,ny,nx);
    
    if (n1*n2*n3 != nx) sf_error("%s: size mismatch",__FILE__);

    sf_adjnull(adj,add,nx,ny,xx,yy);

    for (l=0; l<n3; l++) {
    for (j=0; j < n2; j++) {
	for (i=4; i < n1; i++) {
	    k = i + n1*(j+n2*l);
	
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

    for (l=0; l<n2; l++) {
    for (j=0; j < n3; j++) {
        for (i=4; i < n1; i++) {
            k = i + n1*(j+n3*l);

            if (adj) {
                xx[k-1] += bq[k]*yy[k+4*nx] + aq[k]*yy[k+4*nx] +
                    dq[k]*yy[k+2*nx+4*nx] + cq[k]*yy[k+3*nx+4*nx];
                xx[k-2+n2*n1] -= (2.*bq[k]*dq[k]*(cq[k]*yy[k+4*nx] + aq[k]*yy[k+2*nx+4*nx]) +
                            (bq[k]+2.*aq[k]*cq[k]*dq[k])*yy[k+nx+4*nx] +
                            (dq[k]+2.*aq[k]*cq[k]*bq[k])*yy[k+3*nx+4*nx]);
                xx[k-3+n2*n1] += bq[k]*dq[k]*(dq[k]*yy[k+4*nx]+bq[k]*yy[k+2*nx+4*nx]) +
                    dq[k]*(aq[k]*dq[k]+2.*bq[k]*cq[k])*yy[k+nx+4*nx] +
                    bq[k]*(cq[k]*bq[k]+2.*dq[k]*aq[k])*yy[k+3*nx+4*nx];
                xx[k-4+n2*n1] -= bq[k]*dq[k]*(dq[k]*yy[k+nx+4*nx] + bq[k]*yy[k+3*nx+4*nx]);
            } else {
                yy[k+4*nx] += bq[k]*(xx[k-1] +
                               dq[k]*(dq[k]*xx[k-3+n2*n1] - 2.*cq[k]*xx[k-2+n2*n1]));
                yy[k+nx+4*nx] += aq[k]*xx[k-1+n2*n1]-(bq[k]+2.*aq[k]*cq[k]*dq[k])*xx[k-2+n2*n1] +
                    dq[k]*((aq[k]*dq[k]+2.*bq[k]*cq[k])*xx[k-3+n2*n1] -
                          dq[k]*bq[k]*xx[k-4+n2*n1]);
                yy[k+2*nx+4*nx] += dq[k]*(xx[k-1] +
                                    bq[k]*(bq[k]*xx[k-3+n2*n1] - 2.*aq[k]*xx[k-2+n2*n1]));
                yy[k+3*nx+4*nx] += cq[k]*xx[k-1]-(dq[k]+2.*aq[k]*cq[k]*bq[k])*xx[k-2+n2*n1] +
                    bq[k]*((cq[k]*bq[k]+2.*dq[k]*aq[k])*xx[k-3+n2*n1] -
                          dq[k]*bq[k]*xx[k-4+n2*n1]);
            }
        }
    }
    }
}
