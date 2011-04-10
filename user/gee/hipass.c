/* A simple 1-D high-pass filter */
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
/*^*/

#include "hipass.h"

static float r;

void hipass_init (float eps)
/*< initialize >*/
{
    r = 1+0.5*eps*eps;
    r -= sqrtf(r*r-1.);
}

void hipass_lop(bool adj, bool add, int nx, int ny, float *x, float *y)
/*< linear operator >*/
{
    int i;
    float t;

    sf_adjnull(adj,add,nx,ny,x,y);

    if ( adj) {
	t = y[ny-1];
	x[ny-1] += t;
	for (i=ny-2; i >=0; i--) {
	    t = y[i] - y[i+1] + r*t;
	    x[i] += t;
	}
    } else {
	t = x[0];
	y[0] += t;
	for (i=1; i < nx; i++) {
	    t = x[i] - x[i-1] + r*t;
	    y[i] += t;
	}
    }
}

