/* 1-D finite difference */
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
#include "_bool.h"
/*^*/

#include "igrad1.h"
#include "adjnull.h"

void  sf_igrad1_lop(bool adj, bool add, 
		    int nx, int ny, float *xx, float *yy)
/*< linear operator >*/
{
    int i;

    sf_adjnull(adj,add,nx,ny,xx,yy);
    for (i=0; i < nx-1; i++) {
        if (adj) {
	    xx[i+1] += yy[i];
	    xx[i]   -= yy[i];
	} else {
	    yy[i] += xx[i+1] - xx[i];
	}
    }
}

/* 	$Id$	 */


