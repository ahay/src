/* 2-D omnidirectional plane-wave destruction filter. */
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

#include "opwd2.h"
#include "apfilt.h"

#ifndef _opwd2_h

typedef struct Omni2 *omni2;
/* abstract data type */
/*^*/

#endif

struct Omni2 {
    int nx, ny, nw;
    float *flt, **t1, **t2, *p1, *p2;
};

static omni2 ap2;

omni2 opwd2_init(int nw               /* filter order */, 
		 int nx, int ny       /* data size */, 
		 float *p1, float *p2 /* dip [ny][nx] */) 
/*< Initialize >*/
{
    omni2 ap;
    
    ap = (omni2) sf_alloc(1,sizeof(*ap));
    
    ap->nw = nw;
    ap->nx = nx;
    ap->ny = ny;
    ap->p1 = p1;
    ap->p2 = p2;

    ap->flt = sf_floatalloc(2*nw+1);
    ap->t1 = sf_floatalloc2(nx,ny);
    ap->t2 = sf_floatalloc2(nx,ny);
    
    apfilt_init(nw);
    return ap;
}

void opwd2_close(omni2 ap)
/*< free allocated storage >*/
{
    apfilt_close();
    free(ap->flt);
    free(*(ap->t1)); free(ap->t1);
    free(*(ap->t2)); free(ap->t2);
    free(ap);
}

void opwd22_init (omni2 ap1)
/*< Initialize linear operator >*/
{
    ap2 = ap1;
}

void opwd21 (bool der1, bool der2 /* derivative flags */, 
	     const omni2 ap /* OPWD object */, 
	     float* xx     /* input */, 
 	     float* yy     /* output */)
/*< plane-wave destruction >*/
{
    int ix, iy, iw, is, i;

    for (iy=0; iy < ap->ny; iy++) {
	for (ix=0; ix < ap->nx; ix++) {
	    i = ix + ap->nx * iy;
	    ap->t1[iy][ix] = 0.0f;
	    ap->t2[iy][ix] = 0.0f;
	    yy[i] = 0.0f;
	}
    }
  
    for (iy=0; iy < ap->ny-1; iy++) {
	for (ix = ap->nw; ix < ap->nx-ap->nw; ix++) {
	    i = ix + ap->nx * iy;
	    if (der1) {
		aderfilter(ap->p1[i], ap->flt);
	    } else {
		passfilter(ap->p1[i], ap->flt);
	    }
	      
	    for (iw = 0; iw <= 2*ap->nw; iw++) {
		is = iw-ap->nw;
		  
		ap->t1[iy][ix] += xx[i-is] * ap->flt[iw];
		ap->t2[iy][ix] += xx[i+is] * ap->flt[iw];
	    }
	}
    }

    for (iy = ap->nw; iy < ap->ny-ap->nw; iy++) {
	for (ix=0; ix < ap->nx-1; ix++) {
	    i = ix + ap->nx * iy;
	    if (der2) {
		aderfilter(ap->p2[i], ap->flt);
	    } else {
		passfilter(ap->p2[i], ap->flt);
	    }
	      
	    for (iw = 0; iw <= 2*ap->nw; iw++) {
		is = iw-ap->nw;
		  
		yy[i] += (ap->t2[iy+is][ix] -
			  ap->t1[iy-is][ix]) * ap->flt[iw];
	    }
	}
    }

}

void opwd12 (bool der1, bool der2 /* derivative flags */, 
	     const omni2 ap /* OPWD object */, 
	     float* xx     /* input */, 
 	     float* yy     /* output */)
/*< plane-wave destruction >*/
{
    int ix, iy, iw, is, i;

    for (iy=0; iy < ap->ny; iy++) {
	for (ix=0; ix < ap->nx; ix++) {
	    i = ix + ap->nx * iy;
	    ap->t1[iy][ix] = 0.0f;
	    ap->t2[iy][ix] = 0.0f;
	    yy[i] = 0.0f;
	}
    }
  
    for (iy = ap->nw; iy < ap->ny-ap->nw; iy++) {
	for (ix=0; ix < ap->nx-1; ix++) {
	    i = ix + ap->nx * iy;
	    if (der2) {
		aderfilter(ap->p2[i], ap->flt);
	    } else {
		passfilter(ap->p2[i], ap->flt);
	    }
	      
	    for (iw = 0; iw <= 2*ap->nw; iw++) {
		is = iw-ap->nw;
		  
		ap->t1[iy][ix] += xx[i-is*ap->nx] * ap->flt[iw];
		ap->t2[iy][ix] += xx[i+is*ap->nx] * ap->flt[iw];
	    }
	}
    }

    for (iy=0; iy < ap->ny-1; iy++) {
	for (ix = ap->nw; ix < ap->nx-ap->nw; ix++) {
	    i = ix + ap->nx * iy;
	    if (der1) {
		aderfilter(ap->p1[i], ap->flt);
	    } else {
		passfilter(ap->p1[i], ap->flt);
	    }
	      
	    for (iw = 0; iw <= 2*ap->nw; iw++) {
		is = iw-ap->nw;
		  
		yy[i] += (ap->t2[iy][ix+is] -
			  ap->t1[iy][ix-is]) * ap->flt[iw];
	    }
	}
    }
}
