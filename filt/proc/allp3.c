/* 3-D Plane-wave destruction filter. */
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
#include "apfilt.h"

#ifndef _allp3_h

typedef struct Allpass *allpass;
/* abstract data type */
/*^*/

#endif

struct Allpass {
    int nx, ny, nz, nw, nj;
    float*** pp;
};

static allpass ap1, ap2;

allpass allpass_init(int nw                 /* filter size (1,2,3) */, 
		     int nj                 /* filter step */, 
		     int nx, int ny, int nz /* data size */, 
		     float ***pp            /* data [nz][ny][nx] */)
/*< Initialize >*/
{
    allpass ap;

    ap = (allpass) sf_alloc(1,sizeof(*ap));

    ap->nw = nw;
    ap->nj = nj;
    ap->nx = nx;
    ap->ny = ny;
    ap->nz = nz;
    ap->pp = pp;

    return ap;
}

void allpass1 (bool der         /* derivative flag */, 
	       const allpass ap /* PWD object */, 
	       float*** xx      /* input */, 
	       float*** yy      /* output */)
/*< in-line plane-wave destruction >*/
{
    int ix, iy, iz, iw, is;
    float a[7];

    for (iz=0; iz < ap->nz; iz++) {
	for (iy=0; iy < ap->ny; iy++) {
	    for (ix=0; ix < ap->nx; ix++) {
		yy[iz][iy][ix] = 0.;
	    }
	}
    }
  
    for (iz=0; iz < ap->nz; iz++) {
	for (iy=0; iy < ap->ny-1; iy++) {
	    for (ix = ap->nw*ap->nj; ix < ap->nx-ap->nw*ap->nj; ix++) {
		if (der) {
		    aderfilter(ap->nw, ap->pp[iz][iy][ix], a);
		} else {
		    passfilter(ap->nw, ap->pp[iz][iy][ix], a);
		}
	      
		for (iw = 0; iw <= 2*ap->nw; iw++) {
		    is = (iw-ap->nw)*ap->nj;
		  
		    yy[iz][iy][ix] += (xx[iz][iy+1][ix+is] - 
				       xx[iz][iy  ][ix-is]) * a[iw];
		}
	    }
	}
    }
}

void allpass2 (bool der         /* derivative flag */, 
	       const allpass ap /* PWD object */, 
	       float*** xx      /* input */, 
	       float*** yy      /* output */)
/*< cross-line plane-wave destruction >*/
{
    int ix, iy, iz, iw, is;
    float a[7];
    
    for (iz=0; iz < ap->nz; iz++) {
	for (iy=0; iy < ap->ny; iy++) {
	    for (ix=0; ix < ap->nx; ix++) {
		yy[iz][iy][ix] = 0.;
	    }
	}
    }
    
    for (iz=0; iz < ap->nz-1; iz++) {
	for (iy=0; iy < ap->ny; iy++) {
	    for (ix = ap->nw*ap->nj; ix < ap->nx-ap->nw*ap->nj; ix++) {
		if (der) {
		    aderfilter(ap->nw, ap->pp[iz][iy][ix], a);
		} else {
		    passfilter(ap->nw, ap->pp[iz][iy][ix], a);
		}
		
		for (iw = 0; iw <= 2*ap->nw; iw++) {
		    is = (iw-ap->nw)*ap->nj;
		    
		    yy[iz][iy][ix] += (xx[iz+1][iy][ix+is] - 
				       xx[iz  ][iy][ix-is]) * a[iw];
		}
	    }
	}
    }
}

void allpass3_init (allpass ap, allpass aq)
/*< Initialize linear operator >*/
{
    ap1 = ap;
    ap2 = aq;
}

void allpass3_lop (bool adj, bool add, int n1, int n2, float* xx, float* yy)
/*< PWD as linear operator >*/
{
    int i, ix, iy, iz, iw, is, nx, ny, nz, nw, nj;
    float a[7];

    if (n2 != 2*n1) sf_error("%s: size mismatch: %d != 2*%d",__FILE__,n2,n1);

    sf_adjnull(adj, add, n1, n2, xx, yy);

    nx = ap1->nx;
    ny = ap1->ny;
    nz = ap1->nz;
    nw = ap1->nw;
    nj = ap1->nj;

    if (nx*ny*nz != n1) sf_error("%s: size mismatch",__FILE__);
    
    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny-1; iy++) {
	    for (ix = nw*nj; ix < nx-nw*nj; ix++) {
		i = ix + nx*(iy + ny*iz);

		passfilter(nw, ap1->pp[iz][iy][ix], a);
	      
		for (iw = 0; iw <= 2*nw; iw++) {
		    is = (iw-nw)*nj;
	
		    if (adj) {
			xx[i+nx+is] += yy[i] * a[iw];
			xx[i-is]    -= yy[i] * a[iw];
		    } else {
			yy[i] += (xx[i+nx+is] - xx[i-is]) * a[iw];
		    }
		}
	    }
	}
    }

    nx = ap2->nx;
    ny = ap2->ny;
    nz = ap2->nz;
    nw = ap2->nw;
    nj = ap2->nj;

    if (nx*ny*nz != n1) sf_error("%s: size mismatch",__FILE__);
    
    for (iz=0; iz < nz-1; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix = nw*nj; ix < nx-nw*nj; ix++) {
		i = ix + nx*(iy + ny*iz);

		passfilter(nw, ap2->pp[iz][iy][ix], a);
		
		for (iw = 0; iw <= 2*nw; iw++) {
		    is = (iw-nw)*nj;
		    
		    if (adj) {
			xx[i+nx*ny+is] += yy[i+n1] * a[iw];
			xx[i-is]       -= yy[i+n1] * a[iw];
		    } else {
			yy[i+n1] += (xx[i+nx*ny+is] - xx[i-is]) * a[iw];
		    }
		}
	    }
	}
    }
}

/* 	$Id$	 */
