/* Rational filter by helical convolution and deconvolution. */
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

#include "heldiv.h"
#include "helicon.h"

static int n;
static float *tmp;

void heldiv_init(int nd /* data size */, sf_filter aa, sf_filter bb) 
/*<  Initialized with two filters. >*/
{
    helicon_init(aa);
    sf_polydiv_init(nd,bb);

    n = nd;
    tmp = sf_floatalloc(n);
}

void heldiv_lop(bool adj, bool add, 
		int nx, int ny, float* xx, float*yy) 
/*< linear operator >*/
{
    if (nx != n || ny != n) sf_error("%s: wrong dimensions",__FILE__);
    
    sf_adjnull(adj, add, n, n, xx, yy);
    sf_chain(helicon_lop,sf_polydiv_lop,adj,add,n,n,n,xx,yy,tmp);
}

void heldiv_close(void)
/*< free allocated storage >*/
{
    sf_polydiv_close();
    free(tmp);
}
