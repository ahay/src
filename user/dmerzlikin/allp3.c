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
    float *flt, *pp;
};

static allpass ap1, ap2;

allpass allpass_init(int nw                 /* filter size */, 
		     int nj                 /* filter step */, 
		     int nx, int ny, int nz /* data size */, 
		     float *pp              /* dip [nz*ny*nx] */)
/*< Initialize >*/
{
    allpass ap;

    ap = (allpass) sf_alloc(1,sizeof(*ap));

//^^^ you might need to address nz here 

    ap->nw = nw;
    ap->nj = nj;
    ap->nx = nx;
    ap->ny = ny;
    ap->nz = nz;
    ap->pp = pp;

    ap->flt = sf_floatalloc(2*nw+1);
    apfilt_init(nw);

    return ap;
}

void allpass_close(allpass ap)
/*< free allocated storage >*/
{
    apfilt_close();
    free(ap->flt);
    free(ap);
}

void allpass1 (bool left        /* left or right prediction */,
	       bool der         /* derivative flag */, 
	       const allpass ap /* PWD object */, 
	       float* xx        /* input */, 
	       float* yy        /* output */)
/*< in-line plane-wave destruction >*/
{
    int ix, iy, iz, iw, is, i, nx, ny, nz, i1, i2, ip;

    nx = ap->nx;
    ny = ap->ny;
    nz = ap->nz;

    if (left) {
	i1=1; i2=ny;   ip=-nx;
    } else {
	i1=0; i2=ny-1; ip=nx;
    }

    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		i = ix + nx * (iy + ny * iz);
		yy[i] = 0.;
	    }
	}
    }
  
    for (iz=0; iz < nz; iz++) {
	for (iy=i1; iy < i2; iy++) {
	    for (ix = ap->nw*ap->nj; ix < nx-ap->nw*ap->nj; ix++) {
		i = ix + nx * (iy + ny * iz);

		if (der) {
		    aderfilter(ap->pp[i], ap->flt);
		} else {
		    passfilter(ap->pp[i], ap->flt);
		}
	      
		for (iw = 0; iw <= 2*ap->nw; iw++) {
		    is = (iw-ap->nw)*ap->nj;
		  
		    yy[i] += (xx[i+is+ip] - xx[i-is]) * ap->flt[iw];
		}
	    }
	}
    }
}

void left1 (bool left        /* left or right prediction */,
	       bool der         /* derivative flag */, 
	       const allpass ap /* PWD object */, 
	       float* xx        /* input */, 
	       float* yy        /* output */)
/*< left part of in-line plane-wave destruction >*/
{
    int ix, iy, iz, iw, is, i, nx, ny, nz, i1, i2, ip;

    nx = ap->nx;
    ny = ap->ny;
    nz = ap->nz;

    if (left) {
	i1=1; i2=ny;   ip=-nx;
    } else {
	i1=0; i2=ny-1; ip=nx;
    }

    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		i = ix + nx * (iy + ny * iz);
		yy[i] = 0.;
	    }
	}
    }
  
    for (iz=0; iz < nz; iz++) {
	for (iy=i1; iy < i2; iy++) {
	    for (ix = ap->nw*ap->nj; ix < nx-ap->nw*ap->nj; ix++) {
		i = ix + nx * (iy + ny * iz);

		if (der) {
		    aderfilter(ap->pp[i], ap->flt);
		} else {
		    passfilter(ap->pp[i], ap->flt);
		}
	      
		for (iw = 0; iw <= 2*ap->nw; iw++) {
		    is = (iw-ap->nw)*ap->nj;
		  
		    yy[i] += xx[i+is+ip] * ap->flt[iw];
		}
	    }
	}
    }
}

void right1 (bool left        /* left or right prediction */,
	       bool der         /* derivative flag */, 
	       const allpass ap /* PWD object */, 
	       float* xx        /* input */, 
	       float* yy        /* output */)
/*< right part of in-line plane-wave destruction >*/
{
    int ix, iy, iz, iw, is, i, nx, ny, nz, i1, i2;

    nx = ap->nx;
    ny = ap->ny;
    nz = ap->nz;

    if (left) {
	i1=1; i2=ny;   
    } else {
	i1=0; i2=ny-1;
    }

    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		i = ix + nx * (iy + ny * iz);
		yy[i] = 0.;
	    }
	}
    }
  
    for (iz=0; iz < nz; iz++) {
	for (iy=i1; iy < i2; iy++) {
	    for (ix = ap->nw*ap->nj; ix < nx-ap->nw*ap->nj; ix++) {
		i = ix + nx * (iy + ny * iz);

		if (der) {
		    aderfilter(ap->pp[i], ap->flt);
		} else {
		    passfilter(ap->pp[i], ap->flt);
		}
	      
		for (iw = 0; iw <= 2*ap->nw; iw++) {
		    is = (iw-ap->nw)*ap->nj;
		  
		    yy[i] += xx[i-is] * ap->flt[iw];
		}
	    }
	}
    }
}

void allpass2 (bool left        /* left or right prediction */,
	       bool der         /* derivative flag */, 
	       const allpass ap /* PWD object */, 
	       float* xx        /* input */, 
	       float* yy        /* output */)
/*< cross-line plane-wave destruction >*/
{
    int ix, iy, iz, iw, is, i, nx, ny, nz, i1, i2, ip;

    nx = ap->nx;
    ny = ap->ny;
    nz = ap->nz;

    if (left) {
	i1=1; i2=nz;   ip=-nx*ny;
    } else {
	i1=0; i2=nz-1; ip=nx*ny;
    }
    
    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		i = ix + nx * (iy + ny * iz);
		yy[i] = 0.;
	    }
	}
    }
    
    for (iz=i1; iz < i2; iz++) {
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
		    
		    yy[i] += (xx[i+is+ip] - xx[i-is]) * ap->flt[iw];
		}
	    }
	}
    }
}

void allpass32d_init (allpass ap)
/*< Initialize linear operator >*/
{
    //^^^ deleted second allpass filter structure
    ap1 = ap;
}

void allpass32d_lop (bool adj, bool add, int n1, int n2, float* xx, float* yy)
/*< PWD as linear operator >*/
{
    int i, ix, iy, iz, iw, is, nx, ny, nz, nw, nj;

    if (n2 != n1) sf_error("%s: size mismatch: %d != 2*%d",__FILE__,n2,n1);

    sf_adjnull(adj, add, n1, n2, xx, yy);

    nx = ap1->nx;
    ny = ap1->ny;
    nz = 1;
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
			yy[i] += (xx[i+nx+is] - xx[i-is]) * ap1->flt[iw];
		    }//else
		}//for iw
	    }//for ix
	}//for iy
    }//for iz

	//^^^ deleted second crossline loop
    
}

/* 	$Id$	 */
