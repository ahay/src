/* Data-push binning in 2-D. */
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

static int m1, m2;
static float o1,d1,o2,d2, **xy;

void bin2_init(int m1_in, int m2_in      /* model dimensions */, 
	       float o1_in, float d1_in, 
	       float o2_in, float d2_in  /* model grid */, 
	       float **xy_in             /* data coordinates */)
/*< initialization >*/
{
    m1=m1_in; m2=m2_in;
    o1=o1_in; d1=d1_in;
    o2=o2_in; d2=d2_in;
    xy = xy_in;
}

void bin2_lop (bool adj, bool add, int nm, int nd, float *mm, float *dd)
/*< linear operator >*/
{
    int    i1,i2,im,id;

    if (nm != m1*m2) sf_error("%s: wrong size: %d != %d*%d",
			      __FILE__,nm,m1,m2);

    sf_adjnull(adj,add,nm,nd,mm,dd);

    for (id=0; id < nd; id++) {
        i1 = 0.5 + (xy[0][id]-o1)/d1;
	i2 = 0.5 + (xy[1][id]-o2)/d2;
        if (0<=i1 && i1<m1 &&
            0<=i2 && i2<m2) {
	    im = i1+i2*m1;	    
	    if (adj) mm[im] += dd[id];
	    else     dd[id] += mm[im];
	}
    }
}
