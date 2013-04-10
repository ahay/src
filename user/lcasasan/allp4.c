/* 3-D Plane-wave destruction filter + anisotropy scaling factor */
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
#include <rsfpwd.h>

#include "allp4.h"

#ifndef _allp4_h

typedef struct Allpass4 *allpass4;
/* abstract data type */
/*^*/

#endif

struct Allpass4 {
    int nx, ny, nz, nw, nj;
    float *flt, *pp, *ani;
};

static allpass4 ap1, ap2;

allpass4 allpass4_init(int nw          /* filter size (1,2,3) */, 
		     int nj                  /* filter step */, 
		     int nx, int ny, int nz  /* data size */, 
		     float *pp               /* data [nz*ny*nx] */,
			 float *ani              /* anisotropy [nz*ny*nx] */)
/*< Initialize >*/
{
    allpass4 ap;

    ap = (allpass4) sf_alloc(1,sizeof(*ap));

    ap->nw = nw;
    ap->nj = nj;
    ap->nx = nx;
    ap->ny = ny;
    ap->nz = nz;
    ap->pp = pp;
	ap->ani = ani;

    ap->flt = sf_floatalloc(2*nw+1);
    apfilt_init(nw);

    return ap;
}

void allpass4_close(allpass4 ap)
/*< free allocated storage >*/
{
    apfilt_close();
    free(ap->flt);
    free(ap);
}

void allpass41 (bool der         /* derivative flag */, 
	       const allpass4 ap /* PWD object */, 
	       float* xx        /* input */, 
	       float* yy        /* output */)
/*< in-line plane-wave destruction >*/
{
    int ix, iy, iz, iw, is, i, nx, ny, nz;

    nx = ap->nx;
    ny = ap->ny;
    nz = ap->nz;


    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		i = ix + nx * (iy + ny * iz);
		yy[i] = 0.;
	    }
	}
    }
  
    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny-1; iy++) {
	    for (ix = ap->nw*ap->nj; ix < nx-ap->nw*ap->nj; ix++) {
		i = ix + nx * (iy + ny * iz);

		if (der) {
		    aderfilter(ap->pp[i], ap->flt);
		} else {
		    passfilter(ap->pp[i], ap->flt);
		}
	      
		for (iw = 0; iw <= 2*ap->nw; iw++) {
		    is = (iw-ap->nw)*ap->nj;
		  
		    yy[i] += (xx[i+is+nx] - xx[i-is]) * ap->flt[iw]*ap->ani[i];
		}
	    }
	}
    }
}

void allpass42 (bool der         /* derivative flag */, 
	       const allpass4 ap /* PWD object */, 
	       float* xx        /* input */, 
	       float* yy        /* output */)
/*< cross-line plane-wave destruction >*/
{
    int ix, iy, iz, iw, is, i, nx, ny, nz;

    nx = ap->nx;
    ny = ap->ny;
    nz = ap->nz;
    
    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		i = ix + nx * (iy + ny * iz);
		yy[i] = 0.;
	    }
	}
    }
    
    for (iz=0; iz < nz-1; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix = ap->nw*ap->nj; ix < nx-ap->nw*ap->nj; ix++) {
		i = ix + nx * (iy + ny * iz);
		
		if (der) {
		    aderfilter(ap->pp[i], ap->flt);
		} else {
		    passfilter(ap->pp[i], ap->flt);
		}
		
		for (iw = 0; iw <= 2*ap->nw; iw++) {
		    is = (iw-ap->nw)*ap->nj;
		    
		    yy[i] += (xx[i+is+nx*ny] - xx[i-is]) * ap->flt[iw]*ap->ani[i];;
		}
	    }
	}
    }
}

void allpass43_init (allpass4 ap, allpass4 aq)
/*< Initialize linear operator >*/
{
    ap1 = ap;
    ap2 = aq;
}

void allpass43_lop (bool adj, bool add, int n1, int n2, float* xx, float* yy)
/*< PWD as linear operator >*/
{
    int i, ix, iy, iz, iw, is, nx, ny, nz, nw, nj;

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

		passfilter(ap1->pp[i], ap1->flt);
	      
		for (iw = 0; iw <= 2*nw; iw++) {
		    is = (iw-nw)*nj;
	
		    if (adj) {
			xx[i+nx+is] += yy[i] * ap1->flt[iw];
			xx[i-is]    -= yy[i] * ap1->flt[iw];
		    } else {
			yy[i] += (xx[i+nx+is] - xx[i-is]) * ap1->flt[iw]*ap1->ani[i];
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

		passfilter(ap2->pp[i], ap2->flt);
		
		for (iw = 0; iw <= 2*nw; iw++) {
		    is = (iw-nw)*nj;
		    
		    if (adj) {
			xx[i+nx*ny+is] += yy[i+n1] * ap2->flt[iw];
			xx[i-is]       -= yy[i+n1] * ap2->flt[iw];
		    } else {
			yy[i+n1] += (xx[i+nx*ny+is] - xx[i-is]) * ap2->flt[iw]*ap2->ani[i];
		    }
		}
	    }
	}
    }
}

/* 	$Id: allp3.c 4781 2009-09-25 04:01:38Z sfomel $	 */
