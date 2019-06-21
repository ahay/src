/* Plane-wave destruction with two slopes */
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

#include "twoplane2.h"
#include "allp2.h"

static int nx, ny, nw; 
static allpas2 p, q;
static float *tmp1, *tmp2;

void twoplane2_init (int nw_in            /* filter size */, 
		     int nj1, int nj2     /* dealising stretch */, 
		     int nx_in, int ny_in /* data size */,
		     bool drift           /* if shift filter */,
		     float **pp           /* first slope */, 
		     float **qq           /* second slope */)
/*< initialize >*/
{
    nw = nw_in; nx = nx_in; ny = ny_in;
    p = allpass2_init(nw, nj1, nx, ny, drift, pp);
    q = allpass2_init(nw, nj2, nx, ny, drift, qq);

    nw *= 2*(nj1 > nj2? nj1: nj2);
    
    tmp1 = sf_floatalloc(nx*ny);
    tmp2 = sf_floatalloc(nx*ny);
}

void twoplane2_close(void)
/*< free allocated storage >*/
{
    free (tmp1);
    free (tmp2);
}

void twoplane2_lop (bool adj, bool add, int n1, int n2, float *xx, float *yy)
/*< linear operator >*/
{
    int ix, iy, i;

    if (n1 != n2 || n1 != nx*ny) sf_error("%s: wrong sizes",__FILE__);

    sf_adjnull (adj,add,n1,n2,xx,yy);

    if (adj) {
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		i = ix+iy*nx;

		if (iy < ny-2 && ix >= nw && ix < nx-nw) {
		    tmp1[i] = yy[i];
		} else {
		    tmp1[i] = 0.;
		}
	    }
	}
	
	allpass22_init(q);
	allpass21_lop (true,false,n1,n1,tmp2,tmp1);
	allpass22_init(p);
	allpass21_lop (true,add,n1,n1,xx,tmp2);
    } else {
	allpass22_init(p);
	allpass21_lop (false,false,n1,n1,xx,tmp2);
	allpass22_init(q);
	allpass21_lop (false,false,n1,n1,tmp2,tmp1);

	for (iy=0; iy < ny-2; iy++) {
	    for (ix=nw; ix < nx-nw; ix++) {
		i = ix+iy*nx;
		yy[i] += tmp1[i];
	    }
	}
    }
}
