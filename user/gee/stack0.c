/* Nearest-neighbor NMO and stack */
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

#include "nmo0.h"

static int nt, nx;
static float x0, dx;

void stack0_init(float *slow  /* slowness */,
		 int nt_in,
		 float t0,
		 float dt     /* time axis */,
		 int nx_in,
		 float x0_in,
		 float dx_in  /* offset axis */)
/*< initialize >*/
{
    nt = nt_in;
    nx = nx_in;
    x0 = x0_in;
    dx = dx_in;

    nmo0_init(slow,nt,t0,dt);
}

void stack0_lop(bool adj, bool add, int ns, int ng, 
		float *stack, float *gather)
/*< linear operator >*/
{
    int ix;

    if (ns != nt || ng != nt*nx) sf_error("%s: wrong size",__FILE__);

    sf_adjnull(adj, add, ns, ng, stack, gather);

    for (ix=0; ix < nx; ix++) {
        nmo0_set(x0 + dx * ix); /* set offset */
        nmo0_lop( adj, true, nt, nt, stack, gather+ix*nt);
    }
}
