/* Inverse filtering with a non-stationary helix filter */
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

#include "npolydiv.h"

#include "nhelix.h"
/*^*/

static int nd;
static nfilter aa;
static float *tt;

void npolydiv_init (int nd_in     /* data size */, 
		    nfilter aa_in /* filter */)
/*< initialize >*/
{
    nd = nd_in;
    aa = aa_in;
    tt = sf_floatalloc(nd);
}

void npolydiv_close(void)
/*< free allocated storage >*/
{
    free (tt);
}


void npolydiv_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy)
/*< linear operator >*/
{
    int id, ia, na, ix, iy, ip, *lag;
    float *flt;

    sf_adjnull(adj,add,nx,ny,xx,yy);

    for (id=0; id < nd; id++) {
	tt[id] = adj? yy[id]: xx[id];
    }

    if (adj) {
        for (iy=nd-1; iy >= 0; iy--) { 
	    ip = aa->pch[iy];
	    lag = aa->hlx[ip]->lag; 
	    flt = aa->hlx[ip]->flt;
	    na = aa->hlx[ip]->nh;
	    for (ia=0; ia < na; ia++) {
		ix = iy - lag[ia];     
		if (ix < 0)  continue;
		tt[ix] -=  flt[ia] * tt[iy];
	    } 
	}
	for (id=0; id < nd; id++) {
	    xx[id] += tt[id];
	}
    } else { 
        for (iy=0; iy < nd; iy++) { 
	    ip = aa->pch[iy];
	    lag = aa->hlx[ip]->lag; 
	    flt = aa->hlx[ip]->flt;
	    na = aa->hlx[ip]->nh;
	    for (ia=0; ia < na; ia++) {
		ix = iy - lag[ia]; 
		if (ix < 0)  continue;
		tt[iy] -=  flt[ia] * tt[ix];
	    } 
	}
	for (id=0; id < nd; id++) {
	    yy[id] += tt[id];
        }
    }
}

