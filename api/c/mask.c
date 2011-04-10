/*< Simple mask operator >*/
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

#include "mask.h"

#include "_bool.h"
#include "c99.h"
/*^*/

#include "adjnull.h"
#include "komplex.h"
#include "error.h"

static const bool *m;

void sf_mask_init(const bool *m_in)
/*< initialize with mask >*/
{
    m = m_in;
}

void sf_mask_lop(bool adj, bool add, int nx, int ny, float *x, float *y)
/*< linear operator >*/
{
    int ix;

    if (nx != ny) sf_error("%s: wrong size: %d != %d",nx,ny);

    sf_adjnull (adj,add,nx,ny,x,y);

    for (ix=0; ix < nx; ix++) {
	if (m[ix]) {
	    if (adj) x[ix] += y[ix];
	    else     y[ix] += x[ix];
	}
    }
}

void sf_cmask_lop(bool adj, bool add, int nx, int ny, sf_complex *x, sf_complex *y)
/*< linear operator >*/
{
    int ix;

    if (nx != ny) sf_error("%s: wrong size: %d != %d",nx,ny);

    sf_cadjnull (adj,add,nx,ny,x,y);

    for (ix=0; ix < nx; ix++) {
	if (m[ix]) {
#ifdef SF_HAS_COMPLEX_H	
	    if (adj) x[ix] += y[ix];
	    else     y[ix] += x[ix];
#else
	    if (adj) x[ix] = sf_cadd(x[ix],y[ix]);
	    else     y[ix] = sf_cadd(y[ix],x[ix]);
#endif
	} 
    }
}

/* 	$Id$	 */

