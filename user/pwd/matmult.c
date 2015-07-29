/* Simple matrix multiplication operator */
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

#include "matmult.h"

static float** B;

void matmult_init (float** bb) 
/*< initialize with a pointer to a matrix >*/
{
    B = bb;
}

void matmult_lop (bool adj, bool add, 
		  int nx, int ny, float* x, float*y) 
/*< linear operator >*/
{
    int ix, iy;
    sf_adjnull (adj, add, nx, ny, x, y);
    for (ix = 0; ix < nx; ix++) {
		for (iy = 0; iy < ny; iy++) {
			if (adj) x[ix] += B[iy][ix] * y[iy];
			else     y[iy] += B[iy][ix] * x[ix];
		}
    }
}

