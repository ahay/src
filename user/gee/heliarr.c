/* An array of two helix convolutions */	
/*
  Copyright (C) 2006 University of Texas at Austin
  
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

static sf_filter aa1, aa2;

void heliarr_init (sf_filter a1, sf_filter a2)
/*< initialize >*/
{
    aa1 = a1;
    aa2 = a2;
}

void heliarr_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy)
/*< linear operator >*/
{
    if (ny != 2*nx) sf_error("%s: wrong dimensions: %d != 2*%d",
			     __FILE__,ny,nx);

    sf_adjnull(adj,add,nx,ny,xx,yy);

    sf_helicon_init(aa1);
    sf_helicon_lop (adj,true,nx,nx,xx,yy);
    sf_helicon_init(aa2);
    sf_helicon_lop (adj,true,nx,nx,xx,yy+nx); 
}
