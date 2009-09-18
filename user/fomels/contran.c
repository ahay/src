/* Non-stationary convolution/deconvolution */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

static float **b;
static int* lg, nb, n;
static bool dec;

void contran_init (int nf, float **filter /* [nd][nf] */, int nd, int *lag, bool decon)
/*< initialize >*/
{
    b = filter;
    nb = nf;
    lg = lag; 
    n = nd; 
    dec = decon;
}

void contran_lop (bool adj, bool add, int nx, int ny, float *x, float *y)
/*< apply >*/
{
    int i, ix, iy;
    float t;

    if (nx != n || ny != n) sf_error("%s: Wrong dimensions",__FILE__);

    sf_adjnull(adj,add,nx,ny,x,y);

    if (adj) {
	for (ix=n-1; ix >= 0; ix--) {
	    t = y[ix];

	    for (i=0; i < nb; i++) {
		iy = ix + lg[i];
		if (iy >= n) break; 
		if (iy < 0) continue;

		if (dec) {
		    t -= b[iy][i] * x[iy];
		} else {
		    t += b[iy][i] * y[iy];
		}
	    }
	    x[ix] += t;
	}
    } else {
	for (iy=0; iy < n; iy++) {
	    t = x[iy];
	    for (i=0; i < nb; i++) {
		ix = iy - lg[i];
		if (ix < 0) break;
		if (ix >= n) continue;

		if (dec) {
		    t -= b[iy][i] * y[ix];
		} else {
		    t += b[iy][i] * x[ix];
		}
	    } 
	    y[iy] += t;
	}
    }
}

