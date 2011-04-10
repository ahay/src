/* Dumb inverse NMO in constant velocity */
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

#include <math.h>

#include <rsf.h>
/*^*/

#include "hypotenusei.h"

static int nt, *iz;

void hypotenusei_init(int nt1)
/*< initialize with trace length >*/
{
    nt = nt1;
    iz = sf_intalloc(nt);
}

void hypotenusei_close(void)
/*< free allocated storage >*/
{
    free(iz);
}

void hypotenusei_set(float t0 /* time origin */, 
		     float dt /* time sampling */, 
		     float xs /* offset times slowness */)
/*< set up >*/
{
    int it;
    float t, z2;

    for (it=0; it < nt; it++) {  
	t = t0 + dt*it;
        z2 =  t * t - xs * xs;
	iz[it] = ( z2 >= 0.)? 0.5 + (sqrtf(z2) - t0) /dt: -1;
    }
}

void hypotenusei_lop(bool adj, bool add, 
		     int n1, int n2, float *zz, float *tt)
/*< linear operator >*/
{
    int  it;
    
    sf_adjnull(adj,add,n1,n2,zz,tt);

    for (it=0; it < nt; it++) {  
	if (iz[it] < 0) continue;

	if (adj) 
	    zz[iz[it]] +=  tt[it];
	else 
	    tt[it] +=  zz[iz[it]];
    }
}

