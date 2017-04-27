/* Inverse moveout and spray into a gather */
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

#include "slant.h"

#include <rsf.h>
/*^*/

static int nt, nx, ns;
static float x0, dx, s0, ds, s1, t0, dt, anti;
static float *tmp, *amp, *str, *tx;
static bool rho;

/*------------------------------------------------------------*/
void slant_init (bool rho1                      /* use rho filter */,
		 float x01, float dx1, int nx1  /* offset axis */, 
		 float s01, float ds1, int ns1  /* slowness axis */, 
		 float t01, float dt1, int nt1  /* time axis */, 
		 float s11                      /* reference slowness */, 
		 float anti1                    /* antialiasing */) 
/*< initialize >*/
{
    rho = rho1;

    x0 = x01; dx = dx1; nx = nx1;
    s0 = s01; ds = ds1; ns = ns1;
    t0 = t01; dt = dt1; nt = nt1; 
    s1 = s11; 
    anti = anti1;

    sf_aastretch_init (false, nt, t0, dt, nt);
    if (rho) {
	sf_halfint_init (true,2*nt,1.-1./nt);
	tmp  = sf_floatalloc(nt);
    }

    amp  = sf_floatalloc(nt);
    str  = sf_floatalloc(nt);
    tx   = sf_floatalloc(nt);
}

/*------------------------------------------------------------*/
void slant_close (void)
/*< free allocated storage >*/
{
    free(amp);
    free(str);
    free(tx);

    sf_aastretch_close();
    if (rho) {
	sf_halfint_close();
	free(tmp);
    }
}

/*------------------------------------------------------------*/
void slant_lop (bool adj, 
		bool add, 
		int   nm, 
		int   nd, 
		float *modl, 
		float *data)
/*< linear operator >*/
{
    int ix, is, it;
    float x, s, sxx, t, z;

    if (nm != nt*ns || nd != nt*nx) sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull (adj, add, nm, nd, modl, data);

    for (is=0; is < ns; is++) { /* slowness */ 
	s = s0 + is*ds;
	for (ix=0; ix < nx; ix++) { /* offset */
	    x = x0 + ix*dx;
	    sxx = s*x;

	    for (it=0; it < nt; it++) { /* time */		
		z = t0 + it*dt;
		t = z + sxx;

		str[it] = t;
		tx[it] = anti*fabsf(s-s1)*dx;
		amp[it] = 1.;
	    } /* it */

	    sf_aastretch_define (str, tx, amp);
	    
	    if (rho) {
		if (adj) {
		    sf_halfint_lop (false, false, nt, nt, data+ix*nt, tmp);
		    sf_aastretch_lop (true, true, nt, nt, modl+is*nt, tmp);
		} else {
		    sf_aastretch_lop (false, false, nt, nt, modl+is*nt, tmp);
		    sf_halfint_lop (true, true, nt, nt, data+ix*nt, tmp);
		}
	    } else {
		sf_aastretch_lop(adj,true,nt,nt,modl+is*nt,data+ix*nt);
	    }
	} /* ix */
    } /* is */
}
