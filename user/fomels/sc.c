/* Surface-consistent operations */
/*
  Copyright (C) 2015 University of Texas at Austin
  
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

static int nm, **indx, *size;

void sc_init(int nm1, int **indx1, int *size1)
/*< initialize >*/
{
    nm = nm1;
    indx = indx1;
    size = size1;
}

void sc_lop(bool adj, bool add, int nx, int nd, float* x, float* d)
/*< linear operator >*/
{
    int id, im, ix, sx;
    
    sf_adjnull(adj,add,nx,nd,x,d);

    sx=0;
    for (im=0; im < nm; im++) {
	for (id=0; id < nd; id++) {
	    ix = indx[im][id]+sx;
	    
	    if (adj) {
		x[ix] += d[id];
	    } else {
		d[id] += x[ix];
	    }
	    sx += size[im];
	}
    }
}
