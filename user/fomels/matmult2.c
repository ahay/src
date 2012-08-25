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

#include "matmult2.h"

static float** B;

void matmult2_init (float** bb) 
/*< initialize with a pointer to a matrix >*/
{
    B = bb;
}

void matmult2_lop (bool adj, bool add, 
		  int nx, int ny, float* x, float*y) 
/*< linear operator, no adjoint >*/
{
    int ix, iy;

    for (iy = 0; iy < ny; iy++) {
	y[iy] = 0.;
	for (ix = 0; ix < nx; ix++) {
	    y[iy] += B[iy][ix] * x[ix];
	}
    }
}

void matmult2 (int nx, const float* x, float* y, void* mat) 
/*< square linear operator, no adjoint >*/
{
    float** A;
    int ix, iy;

    A = (float**) mat;

    for (iy = 0; iy < nx; iy++) {
	y[iy] = 0.;
	for (ix = 0; ix < nx; ix++) {
	    y[iy] += A[iy][ix] * x[ix];
	}
    }
}
