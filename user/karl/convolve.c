/* convolution, adjoint is the input */
/*
  Copyright (C) 2014 University of Texas at Austin
  
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

#include "convolve.h"

static int nb;
static const float* bb;

void convolve_init (int nb_in          /* filter size */, 
		    const float* bb_in /* filter [na] */) 
/*< initialize >*/
{
    nb = nb_in;
    bb = bb_in;
}

void convolve_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy) 
/*< linear operator >*/
{
    int b, x, y;

    if(ny < nx+nb-1) sf_error("%s: size problem: %d < %d+%d-1",
			      __FILE__,ny,nx,nb);
    sf_adjnull (adj, add, nx, ny, xx, yy);
    
    for( b=0; b < nb; b++) {
	for( x=0; x < nx; x++) { y = x + b;
	    if( adj) xx[x] += yy[y] * bb[b];
	    else     yy[y] += xx[x] * bb[b];
	}
    }
}

/* 	$Id: convolve.c 7107 2014-12-13 02:04:14Z kschleicher $	 */
