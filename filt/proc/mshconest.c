/* Multi-scale helix convolution for filter estimation */
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

#include "mshconest.h"
#include "hconest.h"

#include "mshelix.h"
/*^*/

static msfilter msaa;

void mshconest_init(float *x, msfilter msaa_in)
/*< initialize with data and filter >*/
{
    msaa = msaa_in;
    hconest_init(x,msaa->one);
}

void mshconest_lop(bool adj, bool add, int na, int ny, float *a, float *y)
/*< linear operator >*/
{
    int  is, nx;
    
    if (na != msaa->nh) sf_error("%s: Wrong data dimensions",__FILE__);
    nx = ny/msaa->ns;

    sf_adjnull(adj, add, na, ny, a, y);

    for (is=0; is < msaa->ns; is++) {
	onescale(is,msaa);
	hconest_lop(adj,true,na,nx,a,y+is*nx);
    }
}

/* 	$Id$	 */
