/*  Tent-like window weighting function */
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

#include <math.h>

#include <rsf.h>

#include "tent.h"

void tent (int dim           /* number of dimensions */, 
	   const int* nwind  /* window size [dim] */, 
	   const int* center /* filter center [dim] */, 
	   const int* a      /* filter size [dim] */, 
	   float* windwt     /* output weight */)
/*< compute weight >*/
{
    int i, j, nw, start[SF_MAX_DIM], end[SF_MAX_DIM], x[SF_MAX_DIM];
    float w, mid[SF_MAX_DIM], wid[SF_MAX_DIM];

    nw = 1;
    for (j=0; j < dim; j++) {
	start[j] = a[j]-center[j];
	end[j] = nwind[j]-center[j];
	mid[j]= (end[j]+start[j])/2.;
	wid[j]= (end[j]-start[j]+1.)/2.;
	nw *= nwind[j]; /* compute window size */
    }

    /* loop in the window */
    for (i=0; i < nw; i++) {
	sf_line2cart(dim, nwind, i, x);
    
	windwt[i] = 1.;
	for (j=0; j < dim; j++) {
	    if (x[j] >= start[j] && x[j] <= end[j]) {
		w = (x[j]-mid[j])/wid[j];
		windwt[i] *= SF_MAX(0.,1.-fabs(w));
	    }	else {
		windwt[i] = 0.;
	    }
	}
    }
}

