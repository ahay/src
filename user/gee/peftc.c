/* Convolution array */
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

#include "peftc.h"

#include "tcai1.h"
#include "tcaf1.h"

static int n;

void peftc_init (int na    /* filter length */, 
		 int ny    /* data length */, 
		 float *aa /* filter [na] */, 
		 float *yy /* data [ny] */)
/*< initialize >*/
{
    n = ny;
    tcai1_init (na, aa);
    tcaf1_init (ny, yy);
}

void peftc_lop (bool adj, bool add, int nx, int nr, float *x, float *r)
/*< linear operator >*/
{
    sf_adjnull(adj,add,nx,nr,x,r);

    tcai1_lop (adj, true, n,    nr, x,   r);
    tcaf1_lop (adj, true, nx-n, nr, x+n, r);
}
