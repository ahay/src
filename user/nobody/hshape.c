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

#include "conv.h"

static sf_filter aa, bb;
static float *t1, *t2, wt;

void hshape_init( int nd       /* data size */,
		  int ns       /* scaling */,
		  sf_filter ff /* filter */) 
/*< initialize >*/
{
    int is, ia, na;
    sf_filter cc;

    aa = ff;
    na = aa->nh;

    bb =  sf_allocatehelix(na);
    for (ia=0; ia < na; ia++) {
	bb->lag[ia] = aa->lag[ia];
	bb->flt[ia] = aa->flt[ia];
    }

    /* convolve ns times */
    for (is=0; is < ns-1; is++) {
	cc = conv (bb, aa, false);
	sf_deallocatehelix(bb);
	bb = cc;
    }

    if (0==(ns%2)) {
	for (ia=0; ia < bb->nh; ia++) {
	    bb->flt[ia] = - bb->flt[ia];
	}
    }

    wt = 1.0/ns;
    t1 = sf_floatalloc (nd);
    t2 = sf_floatalloc (nd);

    sf_polydiv_init(nd,aa);
    sf_helicon_init(bb);
}

void hshape_lop( bool adj, bool add, 
		 int nx, int ny, float* xx, float*yy)
/*< linear operator >*/
{
    int i;

    if (ny != nx) sf_error("%s: Different size",__FILE__);

    sf_adjnull(adj,add,nx,ny,xx,yy);

    if (adj) {
	for (i=0; i < nx; i++) {
	    t2[i] = wt*yy[i];
	}

	sf_polydiv_lop (true, false, nx, nx, t1, t2);
	sf_helicon_lop(true, true, nx, nx, xx, t1);	
    } else {
	sf_helicon_lop(false, false, nx, nx, xx, t1);
	sf_polydiv_lop (false, false, nx, nx, t1, t2);
	
	for (i=0; i < nx; i++) {
	    yy[i] += wt*t2[i];
	}
    }
}

void hshape_close (void) 
/*< free allocated storage >*/
{
    free (t1);
    free (t2);
    sf_deallocatehelix(bb);
    sf_polydiv_close();
}

