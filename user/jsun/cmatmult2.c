/* Simple complex matrix multiplication operator */
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

#include "cmatmult2.h"

static sf_complex** B;

void cmatmult2_init (sf_complex** bb) 
/*< initialize with a pointer to a matrix >*/
{
    B = bb;
}

void cmatmult2_lop (bool adj, bool add, 
		  int nx, int ny, sf_complex* x, sf_complex*y)
/*< linear operator, no adjoint >*/
{
    int ix, iy;

    for (iy = 0; iy < ny; iy++) {
        y[iy] = sf_cmplx(0.,0.);
	for (ix = 0; ix < nx; ix++) {
	    y[iy] += B[iy][ix] * x[ix];
	}
    }
}

void cmatmult2 (int nx, sf_complex* x, sf_complex* y, void* mat) 
/*< square linear operator, no adjoint >*/
{
    sf_complex** A;
    int ix, iy;

    A = (sf_complex**) mat;

    for (iy = 0; iy < nx; iy++) {
        y[iy] = sf_cmplx(0.,0.);
	for (ix = 0; ix < nx; ix++) {
	    y[iy] += A[iy][ix] * x[ix];
	}
    }
}
