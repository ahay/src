/* Transient convolution, adjoint is filter */
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

#include "tcaf1.h"

static int nx;
static float *xx;

void tcaf1_init(int ny    /* data size */, 
		float* yy /* data [ny] */)
/*< initialize >*/
{
    nx = ny;
    xx = yy;
}

void tcaf1_lop(bool adj, bool add, int nb, int ny, float *bb, float *yy)
/*< linear operator >*/
{
    int x, b, y;

    if(ny < nx+nb-1) sf_error("%s: size problem: %d < %d+%d-1",
			      __FILE__,ny,nx,nb);
    sf_adjnull (adj, add, nb, ny, bb, yy);

    for (b=0; b < nb; b++) {
	for (x=0; x < nx; x++) { y = x + b;
	    if( adj) bb[b] += yy[y] * xx[x];
	    else     yy[y] += bb[b] * xx[x];
        }
    }
}
