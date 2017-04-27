/*
 Copyright (C) 2009 University of Texas at Austin
 
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

#include "newton.h"
/*^*/

static func1 zz, zzp, zzs; 
static float x0, v0;

void traveltime_init(func1 z1 /*z(x)*/, 
					 func1 zp1 /*z'(x)*/,
					 func1 zp2 /*z''(x)*/,
					 float x /*receiver location on the surface*/,
					 float v /*velocity*/)
/*<initialize geometry>*/
{
	zz = z1;
	zzp = zp1;
	zzs = zp2;
	x0 = x;
	v0 = v;
	
}

float traveltime(float x)
/*< traveltime >*/
{
	float t;
	
	t = hypotf(x-x0,zz(x))/v0;
	sf_warning("Location of s&r=%g, Depth(f(x))=%g and t=%g",x0,zz(x),t);

	return t;
}


float dtdx(float x)
/*< The first derivative of traveltime >*/
{
  float num, den, diff;
      
	num =  (x-x0)+zz(x)*zzp(x);
	den = hypotf(x-x0, zz(x))*v0;
	diff = num/den;
	
	return diff;
}


float d2tdx2(float x)
/*< The second derivative of traveltime >*/
{
	float diff2;
	
	diff2 = (zzs(x)*zz(x)+pow(zzp(x),2)+1)/(v0*hypotf(x-x0,zz(x))) - pow(x-x0+zz(x)*zzp(x),2)/(v0*pow(hypotf(x-x0,zz(x)),3));
	
	return diff2;
}
