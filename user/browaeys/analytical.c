/* Analytical traveltime */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

float analytical(float x1, float z1, float v1, const float *g1, float *p1,
		 float x2, float z2, float v2, const float *g2, float *p2)
/*< return traveltime >*/
{
    float d[2], a, g[2], v, z, dist, grad, disc, time, num, den, r, r1;
    float gg1, gg2, aa;

    d[0] = z2-z1;
    d[1] = x2-x1;

    a = 1.;

    num = v2-g2[0]*d[0]-g2[1]*d[1] - v1;
    den = v1+g1[0]*d[0]+g1[1]*d[1] - v2;

    if (fabsf(den) < fabsf(num)) {
	r = den/num;
	r1 = 1.+r;
	if (fabsf(r1) > 0.) a = 1./r1;
    } else if (fabsf(num) < fabsf(den)) {
	r = num/den;
	r1 = 1.+r;
	if (fabsf(r1) > 0.) a = r/r1;
    }
    
    if (a < 1. && a >= 0.5) {
	a = 1.;
    } else if (a < 0.5 && a > 0.) {
	a = 0.;
    }

    g[0] = g2[0]+a*(g1[0]-g2[0]);
    g[1] = g2[1]+a*(g1[1]-g2[1]);

    v = v1+v2;
    dist = d[0]*d[0]+d[1]*d[1];
    grad = g[0]*g[0]+g[1]*g[1];

    disc = v*v-dist*grad;
    if (disc > 0.) {
	z = 4*dist/(v+sqrtf(disc));
    } else {
	z = 4*dist/v;
    }

    p1[0] = d[0]-0.25*z*g[0];
    p1[1] = d[1]-0.25*z*g[1];
    disc = 1./hypotf(p1[0],p1[1]);
    p1[0] *= disc;
    p1[1] *= disc;

    p2[0] = d[0]+0.25*z*g[0];
    p2[1] = d[1]+0.25*z*g[1];    
    disc = 1./hypotf(p2[0],p2[1]);
    p2[0] *= disc;
    p2[1] *= disc;
    
    if (z > 0.) {
	time = dist + z*z*grad/48.;

	if (a > 1. && a < 0.) {
	    gg1 = (g1[0]-g2[0])*g1[0]+(g1[1]-g2[1])*g1[1];
	    gg2 = (g1[0]-g2[0])*g2[0]+(g1[1]-g2[1])*g2[1];

	    aa = (a-1.)*a;

	    time += aa*z*z*
		(-2.5*gg2 - 
		 (a*((7  + 2*a*(19 + 5*a*(8*a*(10 + 3*(a-3)*a)-37)))*gg1 + 
		     (13 - 2*a*(64 + 5*a*(8*a*(10 + 3*(a-3)*a)-43)))*gg2))/2. 
		 + 30*aa*aa*
		 (a*(3 - 6*a + 4*a*a)*(gg1 - gg2) + gg2)*logf(a/(a-1)))/120.;
	}
	
	time /= sqrtf(z);
    } else {
	time = 0.;
    }

    return time;
}
