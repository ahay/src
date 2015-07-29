/* Transient convolution with complex numbers, adjoint is filter */
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

#include "ctcaf1.h"

static int nx;
static sf_complex *xx;

void ctcaf1_init(int ny         /* data size */, 
		sf_complex* yy /* data [ny] */)
/*< initialize >*/
{
    nx = ny;
    xx = yy;
}

void ctcaf1_lop(bool adj, bool add, int nb, int ny, sf_complex *bb, sf_complex *yy)
/*< linear operator >*/
{
    int x, b, y;

    if(ny < nx+nb-1) sf_error("%s: size problem: %d < %d+%d-1",
			      __FILE__,ny,nx,nb);
    sf_cadjnull (adj, add, nb, ny, bb, yy);

    for (b=0; b < nb; b++) {
	for (x=0; x < nx; x++) { y = x + b;
#ifdef SF_HAS_COMPLEX_H
	    if( adj) bb[b] += yy[y] * conjf(xx[x]);
	    else     yy[y] += bb[b] * xx[x];
#else
	    if( adj) bb[b] = sf_cadd(bb[b],sf_cmul(yy[y],conjf(xx[x])));
	    else     yy[y] = sf_cadd(yy[y],sf_cmul(bb[b],xx[x]));
#endif
        }
    }
}
