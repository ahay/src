/* B-spline erf function for tapering */
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

#include "berf.h"
#include "sinint.h"

float berf(float x)
/*< B-spline Erf >*/
{
    float b, c, s, x2;

    x *= sqrt(6.0);
    x2 = x*x;

    if (fabsf(x) < SF_EPS) {
	b = x*(2-x2/9.0)/sqrt(6.0*SF_PI);
    } else {
	c = cosf(x);
	s = sinf(x);	

	b = (((c-1)*(1 + x*(s+x) + c*(2*x2-1))) + x*x2*(2*sin_integral(2*x)-sin_integral(x)))*4*sqrt(2.0/(3.0*SF_PI))/(3.*x*x2);
    }

    return b;
}
	
