/* Define a prediction-error filter for steep-dip deconvolution */
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

#include "steepdip.h"

sf_filter steep(int dim    /* number of dimensions */, 
		int *n     /* data size [dim] */, 
		int *a     /* filter size [dim] */, 
		float *d   /* axis sampling [dim] */, 
		float vel  /* velocity */, 
		float tgap /* time gap */) 
/*< define PEF >*/
{
    int *lag, c[SF_MAX_DIM], i, h, na, j;
    float x, t0;
    sf_filter aa;

    na = 1;
    for (j=0; j < dim; j++) {
	na *= a[j];
    }

    lag = sf_intalloc(na);

    for (h=i=0; i < na; i++) { 
	sf_line2cart(dim, a, i, c);

	for (j=0; j < dim-1; j++) {
	    c[j] -= (a[j]-1)/2;
	}

	t0 = 0.;
	for (j=0; j < dim-1; j++) {
	    x = d[j]*c[j];
	    t0 += x*x;
	}
	t0 = sqrtf(t0)/vel;
	if (t0 < tgap) t0 = tgap;

	x = d[dim-1]*c[dim-1];
	if(x >= t0) {
	    lag[h] = sf_cart2line(dim, n, c);
	    h++;
	}
    }

    aa = sf_allocatehelix(h);

    for (i=0; i < h; i++) {
	aa->lag[i] = lag[i];
	aa->flt[i] = 0.;
    }

    free(lag);

    return aa;
}

