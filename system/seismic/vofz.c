/* Traveltime in a V(z) model. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include "kirmod.h"

static float *v, d1;
static int n1, niter;

void vofz_init(int miter /* maximum number of iterations */,
	       int nz    /* number of depth samples */, 
	       float dz  /* depth sampling */,
	       float *vz /* [nz] velocity */)
/*< initialize >*/
{
    v = vz;
    n1 = nz;
    d1 = dz;
    niter = miter;
}

void vofz_table(char   type  /* type of distribution */,
		float z      /* depth */,
		float h      /* offset */,
		ktable table /* output traveltime attributes */)
/*< compute traveltime attributes >*/
{
    int iz, nz, iter;
    float p,vz,t,hp,dhp,zn,vn,cos2,q;

    nz = SF_MIN(floorf(z/d1),n1);
    zn = z-nz*d1;
    vn = (nz < n1)? v[nz]:v[n1-1];

    switch(type) {
	case 'h': /* hyperbolic approximation */


	    t = vz = 0.0;
	    for (iz=0; iz < nz; iz++) {
		vz += v[iz];
		t += 1.0f/v[iz];
	    }
	    if (iz == n1) iz--;

	    vz = vz*d1 + vn*zn; 
	    t  =  t*d1 + zn/vn;
	    
	    table->t = sqrtf(t*(t+h*h/vz));
	    break;
	case 'e': /* exact */
	    p = 0.0;
	    
	    for (iter=0; iter < niter; iter++) {
		hp = 0.0;
		dhp = 0.0;
		
		for (iz=0; iz < nz; iz++) {
		    cos2 = fabsf(1.0f-p*p*v[iz]*v[iz]);
		    q = v[iz]/sqrtf(cos2);

		    hp += q;
		    dhp += q/cos2;
		}
		cos2 = fabsf(1.0f-p*p*vn*vn);
		q = vn/sqrtf(cos2);

		hp = p*(hp*d1+q*zn);
		dhp = dhp*d1+q*zn/cos2;
		
		/* Newton step */			
		p -= (hp-h)/dhp;

		/* Check convergence */
	    }
	    
	    t = 0.0;
	    for (iz=0; iz < nz; iz++) {
		cos2 = fabsf(1.0f-p*p*v[iz]*v[iz]);
		t += 1.0f/(v[iz]*sqrtf(cos2));
	    }
	    cos2 = fabsf(1.0f-p*p*vn*vn);
	    t = t*d1 + zn/(vn*sqrtf(cos2));
	    
	    table->t = t;
	    break;
	default:
	    sf_error("Unknown type");
	    break;
    }
}
