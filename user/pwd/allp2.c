/* 2-D Plane-wave destruction filter. */
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

#include "allp2.h"
#include "apfilt.h"

#ifndef _allp2_h

typedef struct Allpass2 *allpas2;
/* abstract data type */
/*^*/

#endif

struct Allpass2 {
    int nx, ny, nw, nj;
    float *flt, **pp;
};

static allpas2 ap2;

allpas2 allpass2_init(int nw         /* filter order */, 
		       int nj         /* filter step */, 
		       int nx, int ny /* data size */, 
		       float **pp     /* dip [ny][nx] */) 
/*< Initialize >*/
{
    allpas2 ap;
    
    ap = (allpas2) sf_alloc(1,sizeof(*ap));
    
    ap->nw = nw;
    ap->nj = nj;
    ap->nx = nx;
    ap->ny = ny;
    ap->pp = pp;

    ap->flt = sf_floatalloc(2*nw+1);
    
    apfilt_init(nw);
    return ap;
}

void allpass2_close(allpas2 ap)
/*< free allocated storage >*/
{
    apfilt_close();
    free(ap->flt);
    free(ap);
}

void allpass22_init (allpas2 ap1)
/*< Initialize linear operator >*/
{
    ap2 = ap1;
}

void allpass21_lop (bool adj, bool add, int n1, int n2, float* xx, float* yy)
/*< PWD as linear operator >*/
{
    int i, ix, iy, iw, is, nx, ny;

    sf_adjnull(adj,add,n1,n2,xx,yy);
  
    nx = ap2->nx;
    ny = ap2->ny;

    for (iy=0; iy < ny-1; iy++) {
	for (ix = ap2->nw*ap2->nj; ix < nx-ap2->nw*ap2->nj; ix++) {
	    passfilter(ap2->pp[iy][ix], ap2->flt);
	    i = ix + iy*nx;
	      
	    for (iw = 0; iw <= 2*ap2->nw; iw++) {
		is = (iw-ap2->nw)*ap2->nj;
		  
		if (adj) {
		    xx[i+is+nx] += yy[i]*ap2->flt[iw];
		    xx[i-is]    -= yy[i]*ap2->flt[iw];
		} else {
		    yy[i] += (xx[i+is+nx] - xx[i-is]) * ap2->flt[iw];
		}
	    }
	}
    }
}

void allpass21 (bool der          /* derivative flag */, 
		const allpas2 ap /* PWD object */, 
		float** xx        /* input */, 
		float** yy        /* output */)
/*< plane-wave destruction >*/
{
    int ix, iy, iw, is;

    for (iy=0; iy < ap->ny; iy++) {
	for (ix=0; ix < ap->nx; ix++) {
	    yy[iy][ix] = 0.;
	}
    }
  
    for (iy=0; iy < ap->ny-1; iy++) {
	for (ix = ap->nw*ap->nj; ix < ap->nx-ap->nw*ap->nj; ix++) {
	    if (der) {
		aderfilter(ap->pp[iy][ix], ap->flt);
	    } else {
		passfilter(ap->pp[iy][ix], ap->flt);
	    }
	      
	    for (iw = 0; iw <= 2*ap->nw; iw++) {
		is = (iw-ap->nw)*ap->nj;
		  
		yy[iy][ix] += (xx[iy+1][ix+is] - 
			       xx[iy  ][ix-is]) * ap->flt[iw];
	    }
	}
    }
}

/* 	$Id: allp2.c 10160 2013-04-10 20:11:15Z sfomel $	 */
