/* Riesz transform FIR filter. */
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

/* From "Closed-form design of maximally flat FIR Hilbert transformers,
 * differentiators, and fractional delayers by power series expansion" by
 * Soo-Chang Pei and Peng-Hua Wang, IEEE Trans. on Circuits and Systems - Part
 * I: Fundamental theory and applications, v. 48, No. 4, 2001, 389-398. */

#include <rsf.h>

#include "riesz.h"

static float c, c2;
static int n;

void riesz_init(int n1   /* transform length */, 
		float c1 /* filter parameter */)
/*< initialize >*/
{
    n = n1;
    c = 1./(2*sqrtf(c1));
    c2 = c*c;
}

void riesz (int nx, int ny,       /* data size */
	    float** d,            /* input (gets destroyed) */
	    float**dx, float **dy /* output */)
/*< transform >*/
{
    int i, ix, iy;
    
    for (iy=0; iy < ny; iy++) {
		for (ix=0; ix < nx; ix++) {
			dx[iy][ix] = d[iy][ix];
			dy[iy][ix] = 0.;
		}
    }
    
    for (i=n; i >= 1; i--) {
		/* Add proper boundary conditions */

		for (iy=2; iy < ny-2; iy++) {
			for (ix=2; ix < nx-2; ix++) {
				dy[iy][ix] = dx[iy][ix]+(dx[iy+2][ix]+dx[iy][ix+2]-4.*dx[iy][ix]+
										 dx[iy-2][ix]+dx[iy][ix-2])*c2;
			}
		}
		for (iy=0; iy < ny; iy++) {
			for (ix=0; ix < nx; ix++) {
				dx[iy][ix] = d[iy][ix] + dy[iy][ix]*(2*i-1)/(2*i);
			}
		}
    }

    for (iy=0; iy < ny; iy++) {
		for (ix=0; ix < nx; ix++) {
			d[iy][ix] = dx[iy][ix];
		}
    }

    for (iy=0; iy < ny; iy++) {
		dx[iy][0] = 2.*(d[iy][0]-d[iy][1])*c;
		for (ix=1; ix < nx-1; ix++) {
			dx[iy][ix] = (d[iy][ix-1]-d[iy][ix+1])*c;
		}
		dx[iy][nx-1] = 2.*(d[iy][nx-2]-d[iy][nx-1])*c;
    }

    for (ix=0; ix < nx; ix++) {
		dy[0][ix] = 2.*(d[0][ix]-d[1][ix])*c;

		for (iy=1; iy < ny-1; iy++) {
			dy[iy][ix] = (d[iy-1][ix]-d[iy+1][ix])*c;
		}
		dy[ny-1][ix] = 2.*(d[ny-2][ix]-d[ny-1][ix])*c;
    }
}

