/* Cosine window weighting function */
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

#include "tent2.h"

#include "_defs.h"
#include "file.h"
#include "decart.h"

void sf_tent2 (int dim          /* number of dimensions */, 
	       const int* nwind /* window size [dim] */, 
	       float* windwt    /* window weight */)
/*< define window weight >*/
{
    int i, j, nw, x[SF_MAX_DIM];
    double w;

    /* compute window size */
    nw = 1;
    for (j=0; j < dim; j++) {
	nw *= nwind[j];
    }

    /* loop inside the window */
    for (i=0; i < nw; i++) { 
	sf_line2cart(dim, nwind, i, x);
    
	windwt[i] = 1.;
	for (j=0; j < dim; j++) {
	    if (nwind[j] > 1) {
		w = cosf(2.*SF_PI*(x[j]+1.)/(nwind[j] + 1.));
		w = 0.5*(1.-w);
		if (w > 0.) 
		    windwt[i] *= w;
		else
		    windwt[i] = 0.;
	    }
	}
    }
}
