/* Velocity transform */
/*
  Copyright (C) 2011 University of Texas at Austin
  
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

#include "velxf.h"

static int nt, nx, ns;
static float ot, dt, *x2, *z2, *s;
static bool *mask;

void velxf_init(int nt1, int nx1, int ns1,
		float ot1, float dt1,
		float *xx2, float *zz2, float *ss, bool *mmask)
/*< initialize >*/
{
    nt = nt1;
    nx = nx1;
    ns = ns1;

    ot = ot1;
    dt = dt1;

    x2 = xx2;
    z2 = zz2;
    s = ss;
    mask = mmask;
}


void velxf(bool adj, bool add, int nts, int ntx, float* modl, float* data)
/*< linear operator >*/
{
    int is, ix, iz, it, endt;
    float x2s, t, ft, gt;

    if (nts != nt*ns || ntx != nt*nx) sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull( adj, add, nts, ntx, modl, data);

    for (is=0; is < ns; is++) {
	for (ix=0; ix < nx; ix++) {
	    x2s = x2[ix]*s[is];
	    
	    if ( x2s > z2[nt-1] ) break;
	    if (!mask[ix]) continue;
	    
	    endt = floorf(1. + sqrt( z2[nt-1] - x2s )/dt);
	    
	    for (iz=0; iz < endt; iz++) {
		t = sqrtf( z2[iz] + x2s );
		ft = (t-ot)/dt;

		it = floorf(ft);
		ft -= it;
		gt = 1.-ft;
		
		if (adj) {
		    modl[is*nt+iz] += gt*data[ix*nt+it] + ft*data[ix*nt+it+1];
		} else {
		    data[ix*nt+it]   += gt*modl[is*nt+iz];
		    data[ix*nt+it+1] += ft*modl[is*nt+iz];
		}
	    }
	}
    }
}


