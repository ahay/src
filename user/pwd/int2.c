/* 2-D Interpolation in time slices */
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

#include "int2.h"

static int nm, nt, nd;
static float *mod, *dat;

void int2_init (float** coord              /* cooordinates [nd] */, 
		float o1, float o2, 
		float d1, float d2,
		int   n1, int   n2         /* axes */,
		int ns                     /* number of slices */,
		sf_interpolator interp     /* interpolation function */, 
		int nf                     /* interpolator length */, 
		int nd_in                  /* number of data points */)
/*< initialize >*/
{
    nm = n1*n2;
    nt = ns;
    nd = nd_in;

    sf_int2_init(coord,o1,o2,d1,d2,n1,n2,interp,nf,nd);

    mod = sf_floatalloc(nm);
    dat = sf_floatalloc(nd);
}

void int2_lop (bool adj, bool add, int nx, int ny, float* x, float* y)
/*< linear operator >*/
{ 
    int id, im, it;
    
    if (ny != nd*nt) sf_error("%s: wrong data size: %d != %d",__FILE__,ny,nd*nt);
    if (nx != nm*nt) sf_error("%s: wrong data size: %d != %d",__FILE__,nx,nm*nt);

    sf_adjnull (adj,add,nx,ny,x,y);
    
    for (it=0; it < nt; it++) { /* loop over time slices */
	if (adj) {
	    for (id=0; id < nd; id++) {
		dat[id] = y[id*nt+it];
	    }
	} else {
	    for (im=0; im < nm; im++) {
		mod[im] = x[im*nt+it];
	    }
	}

	/* apply interpolation */
	sf_int2_lop(adj,false,nm,nd,mod,dat);

	if (adj) {
	    for (im=0; im < nm; im++) {
		x[im*nt+it] += mod[im];
	    }
	} else {
	    for (id=0; id < nd; id++) {
		y[id*nt+it] += dat[id];
	    }
	} 
    }
}

void int2_close (void)
/*< free allocated storage >*/
{
    free(mod);
    free(dat);
    sf_int2_close();
}

/* 	$Id: int1.c 1855 2006-05-22 17:04:00Z fomels $	 */
