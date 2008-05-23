/* Operations with real sinusoids */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

#include "rsin.h"

static float c0;

void rsin_init(float c1)
/*< initialize >*/
{
    c0 = c1;
}

void rsin_destruct(bool adj, bool add, int nx, int ny, float *xx, float *yy)
/*< destruction operator >*/
{
    int i;

    if (nx != ny) sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull(adj, add, nx, nx, xx, yy);
    
    for (i = 2; i < nx-1; i++) {
	if(adj) {
	    xx[i-1] += yy[i];
	    xx[i+1] += yy[i];
	    xx[i] -= 2*c0*yy[i];
	} else {
	    yy[i] += xx[i-1] + xx[i+1] - 2*c0*xx[i];
	}
    }
}
