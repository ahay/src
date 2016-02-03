/* Weighting with multiple components */
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

#include "_bool.h"
/*^*/

#include "weight2.h"
#include "alloc.h"
#include "error.h"
#include "adjnull.h"

static int nw;
static float **w;

void sf_weight2_init(int nw1   /* number of components */, 
		     int n     /* model size */, 
		     float *ww /* weight [nw*n] */)
/*< initialize >*/
{
    int iw;

    nw = nw1;
    w = (float**) sf_alloc(nw,sizeof(float*));

    for (iw=0; iw < nw; iw++) {
	w[iw] = ww+iw*n;
    }
}

void sf_weight2_close(void)
/*< free allocated storage >*/
{
    free(w);
}

void sf_weight2_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy)
/*< linear operator >*/
{
    int i, iw;

    if (nw*ny != nx) sf_error("%s: size mismatch: %d*%d != %d",
			      __FILE__,nw,ny,nx);

    sf_adjnull (adj, add, nx, ny, xx, yy);
  
    if (adj) {
        for (iw=0; iw < nw; iw++) {
	    for (i=0; i < ny; i++) {
	        xx[i+iw*ny] += yy[i] * w[iw][i];
	    }
	}
    } else {
        for (iw=0; iw < nw; iw++) {
	    for (i=0; i < ny; i++) {
	        yy[i] += xx[i+iw*ny] * w[iw][i];
	    }
	}
    }
}

