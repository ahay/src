/* allp2 for constant dips */
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

#include "callp2.h"
#include "apfilt.h"

#ifndef _callp2_h

typedef struct Callpass2 *callpass2;
/* abstract data type */
/*^*/

#endif

struct Callpass2 {
    int nx, ny, nw, nj;
    float *a, *d;
};

callpass2 callpass2_init(int nw          /* filter size (1,2,3) */, 
			 int nj          /* filter step */, 
			 int nx , int ny /* data size */)
/*< Initialize >*/
{
    callpass2 ap;
    
    ap = (callpass2) sf_alloc(1,sizeof(*ap));
    
    ap->nw = nw;
    ap->nj = nj;
    ap->nx = nx;
    ap->ny = ny;
   
    ap->a = sf_floatalloc(2*nw+1);
    ap->d = sf_floatalloc(2*nw+1);
    apfilt_init(nw);

    return ap;
}

void callpass2_close(callpass2 ap)
/*< Free allocated storage >*/
{
    apfilt_close();
    free(ap->a);
    free(ap->d);
    free(ap);
}

void callpass21_set (callpass2 ap, float p /* constant dip */)
/*< set dip >*/
{
    aderfilter(p, ap->d);
    passfilter(p, ap->a);
}

void callpass21 (bool der           /* derivative flag */, 
		 const callpass2 ap /* PWD object */, 
		 float** xx         /* input */, 
		 float** yy         /* output */)
/*< plane-wave destruction >*/
{
    int ix, iy, iw, is;
    float *b;

    for (iy=0; iy < ap->ny; iy++) {
	for (ix=0; ix < ap->nx; ix++) {
	    yy[iy][ix] = 0.;
	}
    }
  
    b = der? ap->d: ap->a;

    for (iy=0; iy < ap->ny-1; iy++) {
	for (ix = ap->nw*ap->nj; ix < ap->nx-ap->nw*ap->nj; ix++) {
	    for (iw = 0; iw <= 2*ap->nw; iw++) {
		is = (iw-ap->nw)*ap->nj;
		  
		yy[iy][ix] += (xx[iy+1][ix+is] - 
			       xx[iy  ][ix-is]) * b[iw];
	    }
	}
    }
}

/* 	$Id: callp2.c 704 2004-07-13 18:22:06Z fomels $	 */
