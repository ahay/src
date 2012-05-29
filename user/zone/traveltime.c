/*
 *  traveltime.c
 *  Program_1
 *
 *  Created by Yanadet Sripanich on 1/31/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */
#include <rsf.h>

#include "traveltime.h"

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
	sf_warning("%g",t);

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
	
	diff2 = (zz(x)*pow(hypotf(x-x0, zz(x)),2)*zzs(x)+pow((x0-x)*zzp(x)+zz(x),2))/(v0*pow(hypotf(x-x0, zz(x)),3));
	
	return diff2;
}
