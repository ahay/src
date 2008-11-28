/* Helical shaping. */
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
/*^*/

#include "polydiv.h"

static sf_filter aa, bb;
static float* tt;

void hshape_init( int nd       /* data size */,
		  int ns       /* scaling */,
		  sf_filter ff /* filter */) 
/*< initialize >*/
{
    int ia, na;

    aa = ff;

    na = aa->nh;
    bb = sf_allocatehelix(na);
    free (bb->flt);
    bb->flt = aa->flt;
    for (ia=0; ia < na; ia++) {
	bb->lag[ia] = ns * aa->lag[ia];
    }

    tt = sf_floatalloc (nd);
    polydiv_init(nd,aa);
    sf_helicon_init(bb);
}

void hshape_lop( bool adj, bool add, 
		 int nx, int ny, float* xx, float*yy)
/*< linear operator >*/
{
    sf_chain(sf_helicon_lop,polydiv_lop,adj,add,nx,ny,nx,xx,yy,tt);
}

void hshape_close (void) 
/*< free allocated storage >*/
{
    free (tt);
    free (bb->lag);
    free (bb);
    polydiv_close();
}

