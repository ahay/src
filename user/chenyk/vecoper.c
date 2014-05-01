/* Vector operator */
/*
  Copyright (C) 2013 University of Texas at Austin
  
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
#include "vecoper.h"

#ifndef EPS
#define EPS 0.00000000000000000000001
#endif
static int i;

float dotmultsum( float *xx, float *yy, int n)
/*< dot multiplication of vectors and sum up >*/
{	
	float zz=0;
	for(i=0;i<n;i++)
		zz+=xx[i]*yy[i];
	return zz;
}

void scale( float *xx, float *yy, int n, float s)
/*< implement a scaled indentity operator >*/
{
	for(i=0;i<n;i++)
		yy[i]=s*xx[i];
}

void vecabs( float *xx, float *yy, int n)
/*< Absolute value for vector >*/
{
	for(i=0;i<n;i++)
		yy[i]=fabsf(xx[i]);
}

void scalesum(float *xx, float *yy, float *zz, int n, float sx, float sy)
/*< summation between two scaled vector >*/
{
	for(i=0;i<n;i++)
		zz[i]=sx*xx[i]+sy*yy[i];
}

void scalesumreal(float *xx, float yy, float *zz, int n, float sx)
/*< summation between one scaled vector and one scaler  >*/
{
	for(i=0;i<n;i++)
		zz[i]=sx*xx[i]+yy;
}

void vecdiv(float *xx, float *yy, float *zz, int n)
/*< division between two vector >*/
{
	for(i=0;i<n;i++)
		zz[i]=xx[i]/(yy[i]+EPS);
}

void vecmul(float *xx, float *yy, float *zz, int n)
/*< multiplication between two vector >*/
{
	for(i=0;i<n;i++)
		zz[i]=xx[i]*yy[i];
}

float vecmax(float *xx, int n)
/*< maximum value in a vector >*/
{	
	float t=-SF_HUGE;
	for(i=0;i<n;i++)
		if(t<xx[i]) t=xx[i];
	return t;
}

void consvec(float v, int n, float *x)
/*< Create constant-value vector >*/
{	
	for(i=0;i<n;i++)
		x[i]=v;
}



