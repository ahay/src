/* Convert a helix filter from one data size to another */
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

#include "regrid.h"

#include "helix.h"
/*^*/

void regrid( int dim         /* number of dimensions */, 
	     const int* nold /* old data size [dim] */, 
	     const int* nnew /* new data size [dim] */, 
	     filter aa       /* modified filter */) 
/*< change data size >*/
{
    int i, h0, h1, h, ii[SF_MAX_DIM];

    for (i=0; i < dim; i++) {
	ii[i] = nold[i]/2-1;
    }
  
    h0 = sf_cart2line( dim, nold, ii); /* lag of near middle point on nold */
    h1 = sf_cart2line( dim, nnew, ii); /* lag                      on nnew */
    for (i=0; i < aa->nh; i++) { /* forall given filter coefficients */
	h = aa->lag[i] + h0;
	sf_line2cart( dim, nold, h, ii);
	aa->lag[i] = sf_cart2line( dim, nnew, ii) - h1;
    }
}

/* 	$Id$	 */

