/*
 *  PStraveltime.c
 *  Pre-Stack
 *
 *  Created by Yanadet Sripanich on 2/28/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include <math.h>
#include "ps_traveltime.h"
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
	
	diff = ((x-xs)+zz(x)*zzp(x))*pow(hypotf(x-xs, zz(x)),-1)/v0 + ((x-xr)+zz(x)*zzp(x))*pow(hypotf(x-xr, zz(x)),-1)/v0;
	
	return diff;
}


float d2tdx2(float x)
/*<The second derivative of traveltime>*/
{
	float diff2;
	
	diff2 = (zz(x)*pow(hypotf(x-xs, zz(x)),2)*zzs(x)+pow((xs-x)*zzp(x)+zz(x),2))/(v0*pow(hypotf(x-xs, zz(x)),3))+(zz(x)*pow(hypotf(x-xr, zz(x)),2)*zzs(x)+pow((xr-x)*zzp(x)+zz(x),2))/(v0*pow(hypotf(x-xr, zz(x)),3));
	
	return diff2;
}
		
		
