/* Double polynomial division with a non-stationary helical filter */
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

#include "npolydiv2.h"
#include "npolydiv.h"

#include "nhelix.h"
/*^*/

static float *tt;

void npolydiv2_init (int nd     /* data size */, 
		     nfilter aa /* non-stationary filter */)
/*< initialize >*/
{
    npolydiv_init (nd, aa);
    tt = sf_floatalloc(nd);
}

void npolydiv2_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy)
/*< linear operator >*/
{
    if (nx != ny) sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull(adj,add,nx,ny,xx,yy);

    if (adj) {
        npolydiv_lop (false,false,nx,ny,yy,tt);
        npolydiv_lop (true,true,nx,ny,xx,tt);
    } else {
        npolydiv_lop (false,false,nx,ny,xx,tt);
        npolydiv_lop (true,true,nx,ny,yy,tt);
    }
}

void npolydiv2_close(void)
/*< free allocated storage >*/
{
    npolydiv_close();
    free(tt);
}
