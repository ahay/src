/* Simple identity (copy) operator */
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

#include "c99.h"
#include "_bool.h"
/*^*/

#include "error.h"
#include "adjnull.h"
#include "komplex.h"

#include "ccopy.h"

void sf_ccopy_lop (bool adj, bool add, int nx, int ny, 
		   sf_complex* xx, sf_complex* yy)
/*< linear operator >*/
{
    int i;
    
    if (ny!=nx) sf_error("%s: size mismatch: %d != %d",__FILE__,ny,nx);

    sf_cadjnull (adj, add, nx, ny, xx, yy);
  
    for (i=0; i < nx; i++) {
	if (adj) {
#ifdef SF_HAS_COMPLEX_H
	    xx[i] += yy[i];
#else
	    xx[i] = sf_cadd(xx[i],yy[i]);
#endif
	} else {
#ifdef SF_HAS_COMPLEX_H
	    yy[i] += xx[i];
#else
	    yy[i] = sf_cadd(yy[i],xx[i]);
#endif
	}
    }
}

/* 	$Id: ccopy.c 4396 2009-04-29 02:16:06Z jennings_jim $	 */
