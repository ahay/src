/* Nearest-neighbor NMO with weighting */
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

static float *slow, t0, dt, x;
static int n;

void nmo1_init(float *slow_in /* slowness */,
	       int n_in,
	       float t0_in,
	       float dt_in    /* time axis */)
/*< initialize >*/
{
    slow = slow_in;
    n = n_in;
    t0 = t0_in;
    dt = dt_in;
}

void nmo1_set(float x_in /* offset */)
/*< set offset >*/
{
    x = x_in;
}

void nmo1_lop(bool adj, bool add, int nz, int nt, float *zz,  float *tt)
/*< linear operator >*/
{
    int  it, iz;
    float  xs, t, z, wt;

    if (nt != n || nz != n) sf_error("%s: wrong size",__FILE__);

    sf_adjnull(adj, add, n, n, zz, tt);

    for (iz=0; iz < n; iz++) {   
	z = t0 + dt*iz;            /* Travel-time depth */
	xs= x * slow[iz];
	t = hypotf(z,xs) + 1.e-20; /* Hypotenuse */
	wt = z/t * (1./sqrt(t));   /* Weighting function */
	it = 0.5 + (t - t0) / dt;  /* Round to nearest neighbor. */
	if( it < n ) {
	    if( adj ) zz[iz] += tt[it] * wt;
	    else      tt[it] += zz[iz] * wt;
	}
    }
}
