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

#include <stdio.h>
#include <math.h>

#include "newton.h"
#include <rsf.h>
/*^*/

static func1 zz, zzp, zzs; 
static float v0, xs, xr;

void traveltime_init(func1 z1 /*z(x)*/, 
				   func1 zp1 /*z'(x)*/,
				   func1 zp2 /*z''(x)*/,
				   float s /*the position of source on the surface in km*/,
				   float r /*the position of receiver on the surface in km*/,
				   float v /*velocity */)
/*<initialize geometry>*/
{
	zz = z1;
	zzp = zp1;
	zzs = zp2;
	xs = s;
	xr = r;
	v0 = v;
}


float traveltime(float x)
/*<traveltime>*/
{
	float t;
	
	t = hypotf(x-xs,zz(x))/v0 + hypotf(x-xr,zz(x))/v0;
	sf_warning("x=%g,xs=%g,xr=%g and t=%g",x,xs,xr,t);
	
	return t;
}


float dtdx(float x)
/*<The first derivative of traveltime>*/
{
	float diff;
	
	diff = ((x-xs)+zz(x)*zzp(x))/(v0*hypotf(x-xs, zz(x))) + ((x-xr)+zz(x)*zzp(x))/(v0*hypotf(x-xr, zz(x)));
	
	return diff;
}


float d2tdx2(float x)
/*<The second derivative of traveltime>*/
{
	float diff2;
	
	diff2 = (pow(xs-x,2)*zzp(x)*zzp(x)+(xs-x)*zz(x)*((xs-x)*zzs(x)+2*zzp(x))+pow(zz(x),3)*zzs(x)+zz(x)*zz(x))/(v0*pow(hypotf(x-xs, zz(x)),3)) + (pow(xr-x,2)*zzp(x)*zzp(x)+(xr-x)*zz(x)*((xr-x)*zzs(x)+2*zzp(x))+pow(zz(x),3)*zzs(x)+zz(x)*zz(x))/(v0*pow(hypotf(x-xr, zz(x)),3));
	
	return diff2;
}
		
		
