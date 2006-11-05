/* Slow Kirchhoff from J.F.Claerbout */
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

static float velhalf, t0,dt,dx;
static int nt,nx;

void kirchslow_init(float velhalf_in          /* half velocity */, 
		    float t0_in, float dt_in  /* time axis */,
		    float dx_in               /* midpoint axis */, 
		    int nt_in, int nx_in      /* time-midpoint dimensions */)
{
    velhalf = velhalf_in;
    t0 = t0_in;
    dt = dt_in;
    dx = dx_in;
    nt = nt_in;
    nx = nx_in;
}

void kirchslow(bool adj, bool add, int nm, int nd, float *modl, float *data)
{
    int ix,iy,it,iz, im,id;
    float t,z;

    sf_adjnull(adj,add,nm,nd,modl,data);

    for (ix=0; ix < nx; ix++) { 
	for (iy=0; iy < nx; iy++) { 
	    for (iz=0; iz < nt; iz++) { 
		z = t0+dt*iz; /* travel-time depth */
		t = hypotf(z,(ix-iy)*dx / velhalf);
		it = 0.5 + (t-t0) / dt;
		id = it + iy*nt;
		im = iz + ix*nt;
		
		if( it < nt ) {
		    if( adj) modl[im] += data[id];
		    else     data[id] += modl[im];
		}
	    }
	}
    }
}

