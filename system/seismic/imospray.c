/* Dumb modeling by inverse NMO */
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

#include "hypotenusei.h"

static int nt, nx;
static float x0,dx, t0,dt;

void imospray_init (float slow /* slowness */, 
		    float y0   /* offset origin */, 
		    float dy   /* offset sampling */, 
		    float z0   /* time origin */, 
		    float dz   /* time sampling */, 
		    int nz     /* time samples */, 
		    int ny     /* offset samples */)
/*< initialize >*/
{
    x0 = y0*slow;
    dx = dy*slow;
    t0 = z0;
    dt = dz;
    nt = nz;
    nx = ny;

    hypotenusei_init (nt);
}

void imospray_lop(bool adj, bool add, int n1, int n2, 
		  float *stack, float *gather)
/*< linear operator >*/
{
    int ix;
    float x;

    sf_adjnull(adj,add,n1,n2,stack,gather);

    for (ix=0; ix < nx; ix++) {
	x = x0 + dx*ix;

        hypotenusei_set (t0, dt, x);
	hypotenusei_lop (adj, true, nt, nt, stack, gather+ix*nt);
    }
}

void imospray_close(void)
/*< free allocated storage >*/
{
    hypotenusei_close ();
}
