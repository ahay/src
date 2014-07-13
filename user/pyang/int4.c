/* 4-D spatial interpolation in time slices */
/*
  Copyright (C) 2014 Xi'an Jiaotong University, UT Austin (Pengliang Yang)
  
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
#include <complex.h>

#include "int4.h"
#include "lint4.h"

static int nm, nt, nd;
static sf_complex *mod, *dat;

void int4_init (float o1, float o2, float o3, float o4,
		float d1, float d2, float d3, float d4,
		int   n1, int   n2, int   n3, float n4,
		int nt_,  int nd_in/*dimension product besides time for model*/,
		float* x, float* y, float *p, float *q /* coordinates */)
/*< initialize >*/
{
    nm = n1*n2*n3*n4;
    nt = nt_;
    nd = nd_in;
    lint4_init(n1,o1,d1,n2,o2,d2,n3,o3,d3,n4,o4,d4,x,y,p,q) ;

    mod = sf_complexalloc(nm);
    dat = sf_complexalloc(nd);
}

void int4_lop (bool adj, bool add, int nx, int ny, sf_complex *x, sf_complex *y)
/*< linear operator >*/
{ 
    int id, im, it;
    
    if (ny != nd*nt) sf_error("%s: wrong data size: %d != %d",__FILE__,ny,nd*nt);
    if (nx != nm*nt) sf_error("%s: wrong data size: %d != %d",__FILE__,nx,nm*nt);

    sf_cadjnull (adj,add,nx,ny,x,y);
    
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
	lint4_lop(adj,false,nm,nd,mod,dat);

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

void int4_close (void)
/*< free allocated storage >*/
{
    free(mod);
    free(dat);
}
