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
#include <rsfpwd.h>

#include "allp3norm.h"

#ifndef _allp3norm_h

typedef struct Allpassn *allpassn;
/* abstract data type */
/*^*/

#endif

struct Allpassn {
    int nx, ny, nz, nw, nj;
    float *flt, *pp;
};



allpassn allpassn_init(int nw         /* filter size (1,2,3) */, 
		     int nj                 /* filter step */, 
		     int nx, int ny, int nz /* data size */, 
		     float *pp              /* data [nz*ny*nx] */)
/*< Initialize >*/
{
    allpassn ap;

    ap = (allpassn) sf_alloc(1,sizeof(*ap));

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

void allpassn_close(allpassn ap)
/*< free allocated storage >*/
{
    apfilt_close();
    free(ap->flt);
    free(ap);
}

void allpassn1 (bool norm         /* normalization flag */, 
	       const allpassn ap /* PWD object */, 
	       float* xx        /* input */, 
	       float* yy        /* output */)
/*< in-line plane-wave destruction >*/
{
    int ix, iy, iz, iw, is, i, nx, ny, nz;
	float E;

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
	
		passfilter(ap->pp[i], ap->flt);
		E=0.0;	      
		/* filter energy normalization*/	
		if (norm) {			
			for (iw = 0; iw <= 2*ap->nw; iw++) 
		    	E+= 2.0*ap->flt[iw]* ap->flt[iw];

		//sf_warning("P=%f,ap[1]=%f,E=%f,ap_N[1]=%f",ap->pp[i],ap->flt[1],E,ap->flt[1] / sqrtf(E) );	


			for (iw = 0; iw <= 2*ap->nw; iw++) 
		    	ap->flt[iw]=ap->flt[iw]/sqrtf(E);
		}

		//sf_warning("k=%f",sqrtf(E));
		for (iw = 0; iw <= 2*ap->nw; iw++) {
		    is = (iw-ap->nw)*ap->nj;
		  
		    yy[i] += (xx[i+is+nx] - xx[i-is]) 
					* ap->flt[iw];
			}
		//sf_warning("y[%d]=%f",i,yy[i]);
	    }
	}
    }
}

void allpassn2 (bool norm         /* normliaztion flag */, 
	       const allpassn ap /* PWD object */, 
	       float* xx        /* input */, 
	       float* yy        /* output */)
/*< cross-line plane-wave destruction >*/
{
    int ix, iy, iz, iw, is, i, nx, ny, nz;
	float E;

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
	
		passfilter(ap->pp[i], ap->flt);
		E=0.0;

		/* filter energy normalization*/	
		if (norm) {			
			for (iw = 0; iw <= 2*ap->nw; iw++) 
		    	E+= 2.0*ap->flt[iw]* ap->flt[iw];
			
			for (iw = 0; iw <= 2*ap->nw; iw++) 
		    	ap->flt[iw]=ap->flt[iw]/sqrtf(E);
		}
		
		for (iw = 0; iw <= 2*ap->nw; iw++) {
		    is = (iw-ap->nw)*ap->nj;
		    
		    yy[i] += (xx[i+is+nx*ny] - xx[i-is]) 
					* ap->flt[iw];
		}
	    }
	}
    }
}

