/* Common phase-shift computation. */ 
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

#include "pshift.h"

static bool depth;
static char rule;
static float dz;

void pshift_init(bool depth1 /* depth (or time) */,
		 float dz1 /* depth step */, 
		 char rule1 /* interpolation rule */)
/*< Initialize >*/
{
    depth = depth1;
    rule = rule1;
    dz = dz1;
}

float complex pshift(float complex w2, float k2, float v1, float v2)
/*< phase shift for different rules >*/
{
    float complex cshift, cshift1, cshift2, y;
    float x;

    w2 *= w2;

    switch (rule) {
	case 's': /* simple */			
	    if (depth) {
		w2 = w2 * v1 + k2;
	    } else {
		w2 = w2 + v1 * k2;
	    }
	    cshift = csqrtf(w2);
	    break;
	case 'm': /* midpoint */
	    if (depth) {
		w2 = 0.5 * w2 * (v1+v2) + k2;
	    } else {
		w2 = w2 + 0.5 * (v1+v2) * k2;
	    }
	    cshift = csqrtf(w2);
	    break;
	case 'l': /* linear slowth */ 
	default:
	    if (depth) {
		cshift1 = csqrtf(w2 * v1 + k2);
		cshift2 = csqrtf(w2 * v2 + k2);
			
		cshift = (cshift1 + cshift2 - 
			  1./(1./cshift1+1./cshift2))/1.5;
	    } else {
		cshift1 = csqrtf(w2 + v1 * k2);
		cshift2 = csqrtf(w2 + v2 * k2);
			
		x = 1./(1.+v2/v1);
		cshift = cshift1 + x*(cshift2-cshift1);
		y = x*cshift2/cshift;
		cshift *= (1.-y*(1.-y))/(1.-x*(1.-x));
	    }
	    break;
    }
    return cexpf(-cshift*dz);
}
