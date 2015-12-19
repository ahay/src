/* 2 x 2 matrix multiplications */
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

#include "twobytwo.h"

static float **w;

void twobytwo_init(float **ww)
/*< initialize >*/
{
    w = ww;
}

void twobytwo_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy)
/*< linear operator >*/
{
    int i, n;

    if (ny != nx) sf_error("%s: size mismatch: %d != %d",__FILE__,ny,nx);

    n = nx/2;
    
    sf_adjnull (adj, add, nx, ny, xx, yy);
  
    if (adj) {
	for (i=0; i < n; i++) {
	    xx[i]   += yy[i] * w[0][i] + yy[i+n] * w[2][i];
	    xx[i+n] += yy[i] * w[1][i] + yy[i+n] * w[3][i];
	}
    } else {
	for (i=0; i < n; i++) {
	    yy[i]   += xx[i] * w[0][i] + xx[i+n] * w[1][i];
	    yy[i+n] += xx[i] * w[2][i] + xx[i+n] * w[3][i];
	}
    }
}

