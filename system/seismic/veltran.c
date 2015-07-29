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

#include "veltran.h"

#include <rsf.h>
/*^*/

static int nt, nx, ns, psun1, psun2;
static float x0, dx, s0, ds, s1, t0, dt, anti;
static float *tmp, *amp, *str, *tx;
static bool pull;

void veltran_init (bool pull1                     /* pull or push mode */, 
		   float x01, float dx1, int nx1  /* offset axis */, 
		   float s01, float ds1, int ns1  /* slowness axis */, 
		   float t01, float dt1, int nt1  /* time axis */, 
		   float s11                      /* reference slowness */, 
		   float anti1                    /* antialiasing */, 
		   int psun11, int psun21         /* amplitude case */)
/*< initialize >*/
{
    pull = pull1;
    x0 = x01; dx = dx1; nx = nx1;
    s0 = s01; ds = ds1; ns = ns1;
    t0 = t01; dt = dt1; nt = nt1; 
   
    s1 = s11; psun1 = psun11; psun2 = psun21;
    anti = anti1;

    sf_aastretch_init (false, nt, t0, dt, nt);
    sf_halfint_init (true,2*nt,1.-1./nt);

    amp  = sf_floatalloc(nt);
    str  = sf_floatalloc(nt);
    tx   = sf_floatalloc(nt);
    tmp  = sf_floatalloc(nt);
}

void veltran_close (void)
/*< free allocated storage >*/
{
    free(amp);
    free(str);
    free(tx);
    free(tmp);

    sf_aastretch_close();
    sf_halfint_close();
}

void veltran_lop (bool adj, bool add, int nm, int nd, float *modl, float *data)
/*< linear operator >*/
{
    int ix, is, it;
    float x, s, sx, sxx, t, z, w;

    if (nm != nt*ns || nd != nt*nx) sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull (adj, add, nm, nd, modl, data);

    for (is=0; is < ns; is++) { 
	s = s0 + is*ds;
	for (ix=0; ix < nx; ix++) { 
	    x = x0 + ix*dx;
	    sx = (s-s1)*x*dx;
	    sxx = s*x*x;

	    for (it=0; it < nt; it++) {		
		z = t0 + it*dt;
		t = pull? z*z + sxx: z*z - sxx;

		if (t > 0. && z > 0.) {
		    t = sqrtf(t);
		    str[it] = t;
		    tx[it] = fabsf(anti*sx)/t;
		    switch (adj? psun1: psun2) {
			case 2:
			    w = x*x;
			    break;
			case 1:
			    w = fabsf(x);
			    break;
			default:
			    w = 1.;
			    break;
		    }

		    amp[it] = pull? 
			(z/t) * sqrtf(w/t):
			(t/z) * sqrtf(w/z);
		} else {
		    str[it] = t0 - 2.*dt;
		    tx[it] = 0.;
		    amp[it] = 0.;
		}		
	    } /* it */

	    sf_aastretch_define (str, tx, amp);

	    if (pull) {
		sf_chain(sf_halfint_lop,sf_aastretch_lop,
			 adj,true,nt,nt,nt,modl+is*nt,data+ix*nt,tmp);
	    } else {
		sf_chain(sf_aastretch_lop,sf_halfint_lop,
			 (bool) !adj,true,nt,nt,nt,data+ix*nt,modl+is*nt,tmp);
	    }
	} /* ix */
    } /* is */
}
