/* Frequency filtering, 2-D */
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

#include "expont.h"

static int n1, n2;
static float *a, *b;

void expont_init(int n1_in,int n2_in /* data size */, 
		 float *a1           /* filter component [n1*n2] */, 
		 float *b1           /* filter component [n1*n2] */)
/*< initialize >*/
{
    n1 = n1_in;
    n2 = n2_in;
    a = a1;
    b = b1;
}

void expont_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy)
/*< linear operator >*/
{
    int it, ix, i;

    sf_adjnull(adj,add,nx,ny,xx,yy);

    for (ix=0; ix < n2; ix++) {
	for (it=2; it < n1; it++) {
	    i = it + ix*n1;
	    if (adj) {
		xx[i-1] -= a[i]*b[i]*yy[i];
		xx[i]   += 0.5*yy[i];
		xx[i-2] += 0.5*b[i]*b[i]*yy[i];
	    } else {
		yy[i] += 0.5*(xx[i] + b[i]*(b[i]*xx[i-2] - 2.*a[i]*xx[i-1]));
	    }
	}
    }
}
