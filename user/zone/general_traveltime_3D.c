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
#include <stdlib.h>
#include <math.h>
#include "general_traveltime_3D.h"

#ifndef _general_traveltime_3D_h

typedef struct threed {
	float x; /* x-coordinate*/
	float y; /* y-coordinate*/
	float z; /* z-coordinate*/
	float dx1; /* First derivative in x-direction*/
	float dx2; /* Second derivative in x-direction*/
	float dy1; /* First derivative in y-direction*/
	float dy2; /* Second derivative in y-direction*/
	float dxy2; /* Cross-derivative of x and y*/
	float v1; /* Velocity at the reflector from above*/
	float v2; /* Velocity at the reflector from below*/
	float gx1;/* x-direction velocity gradient from above*/
	float gx2;/* x-direction velocity gradient from below*/
	float gy1;/* y-direction velocity gradient from above*/
	float gy2;/* y-direction velocity gradient from below*/
	float gz1;/* z-direction velocity gradient from above*/
	float gz2;/* z-direction velocity gradient from below*/
} threed;
/* Structure pointer */
/*^*/

typedef float (*func1)(int, float, float); 
/* Function pointer for int, float and float (layer,x,y) -> float */
/*^*/

typedef float (*func4)(int, float, float, int); 
/* Function pointer for int, float and float (layer,x,y,m) -> float */
/*^*/

typedef double (*func2)(threed,threed);
/* Function pointer for twod,twod -> double */
/*^*/

typedef struct func3 { 
	func2 T_k;
	func2 T_k_k_1;
	func2 T_k_k_2;
	func2 T_k_k1_1;
	func2 T_k_k1_2;
	func2 T_k_k_k_1;
	func2 T_k_k_k_2;
	func2 T_k_k_k_12;
	func2 T_k_k1_k1_1;
	func2 T_k_k1_k1_2;
	func2 T_k_k1_k1_12;
	func2 T_k_k_k1_1;
	func2 T_k_k_k1_2;
	func2 T_k_k_k1_12;
	func2 T_k_k_k1_21;
	func2 T_k_zk;
	func2 T_k_zk1;
	func2 T_k_zk_zk;
	func2 T_k_zk1_zk1;
	func2 T_k_zk_zk1;
	func2 T_k_k_zk_1;
	func2 T_k_k_zk_2;
	func2 T_k_k1_zk1_1;
	func2 T_k_k1_zk1_2;
	func2 T_k_k_zk1_1;
	func2 T_k_k_zk1_2;
	func2 T_k_k1_zk_1;
	func2 T_k_k1_zk_2;
} func3;
/* Structure pointer */
/*^*/

#endif

static threed y_1k,y_k,y_k1;

/*NOTE T(x_k,x_k+1,z_k,z_k+1) = T_hat(x_k,x_k+1)------------------------------------------------------------------------*/

/* Initialize two segments at a time--------------------------------------------------------------*/

void initialize(int i /* Indicator of layer*/,
				int nr3 /* Number of reflection*/,
				float **xx /* Position*/,
				float *v /* Velocity*/, 
				float *xref /* Reference position*/,
				float *yref /* Reference position*/,
				float *zref /* Reference position*/, 
				float *gx /* x-gradient*/,
				float *gy /* y-gradient*/,
				float *gz /* z-gradient*/,
				func1 z /* z(x,y)*/,
				func4 zder /* z'(x,y)*/,
				func4 zder2_1 /* d2z/dx2 and d2z/dxdy*/,
				func4 zder2_2 /* d2z/dx2 and d2z/dxdy*/) 
/*<Initialize geometry>*/
{
	/* y_1k (y_k-1 th)------------------------------------------------------------------------------*/
	
	y_1k.x = xx[i-1][0];
	y_1k.y = xx[i-1][1];
	y_1k.z = z(i-1,y_1k.x,y_1k.y);
	y_1k.dx1 = zder(i-1,y_1k.x,y_1k.y,0); 
	y_1k.dx2 = zder2_1(i-1,y_1k.x,y_1k.y,0);
	y_1k.dy1 = zder(i-1,y_1k.x,y_1k.y,1); 
	y_1k.dy2 = zder2_2(i-1,y_1k.x,y_1k.y,1);
	y_1k.dxy2 = zder2_1(i-1,y_1k.x,y_1k.y,1); /*or zder2_2(i-1,y_1k.x,y_1k.y,0) (Cross derivative)*/
	
	if (i!=1) {
		y_1k.gx1 = gx[i-2];
		y_1k.gy1 = gy[i-2];
		y_1k.gz1 = gz[i-2];
		y_1k.v1 = v[i-2]+y_1k.gx1*(y_1k.x-xref[i-2])+y_1k.gy1*(y_1k.y-yref[i-2])+y_1k.gz1*(y_1k.z-zref[i-2]);
		y_1k.gx2 = gx[i-1];
		y_1k.gy2 = gy[i-1];
		y_1k.gz2 = gz[i-1];
		y_1k.v2 = v[i-1]+y_1k.gx2*(y_1k.x-xref[i-1])+y_1k.gy2*(y_1k.y-yref[i-1])+y_1k.gz2*(y_1k.z-zref[i-1]);
		
	} else if (i==1) { /* For the air above at the first reflection*/
		
		y_1k.gx1 = 0;
		y_1k.gy1 = 0;
		y_1k.gz1 = 0;
		y_1k.v1 = 0;
		y_1k.gx2 = gx[i-1];
		y_1k.gy2 = gy[i-1];
		y_1k.gz2 = gz[i-1];
		y_1k.v2 = v[i-1]+y_1k.gx2*(y_1k.x-xref[i-1])+y_1k.gy2*(y_1k.y-yref[i-1])+y_1k.gz2*(y_1k.z-zref[i-1]);
	}
	
	/* y_k----------------------------------------------------------------------------------------*/
	
	y_k.x = xx[i][0];
	y_k.y = xx[i][1];
	y_k.z = z(i,y_k.x,y_k.y);
	y_k.dx1 = zder(i,y_k.x,y_k.y,0);
	y_k.dx2 = zder2_1(i,y_k.x,y_k.y,0);
	y_k.dy1 = zder(i,y_k.x,y_k.y,1);
	y_k.dy2 = zder2_2(i,y_k.x,y_k.y,1);
	y_k.dxy2 = zder2_1(i,y_k.x,y_k.y,1);
	
	y_k.gx1 = gx[i-1];
	y_k.gy1 = gy[i-1];
	y_k.gz1 = gz[i-1];
	y_k.v1 = v[i-1]+y_k.gx1*(y_k.x-xref[i-1])+y_k.gy1*(y_k.y-yref[i-1])+y_k.gz1*(y_k.z-zref[i-1]); /*Of the layer from above*/
	y_k.gx2 = gx[i];
	y_k.gy2 = gy[i];
	y_k.gz2 = gz[i];
	y_k.v2 = v[i]+y_k.gx2*(y_k.x-xref[i])+y_k.gy2*(y_k.y-yref[i])+y_k.gz2*(y_k.z-zref[i]); /*Of the layer from below*/
	
	/* y_k1 (y_k+1 th)----------------------------------------------------------------------------*/

	y_k1.x = xx[i+1][0];
	y_k1.y = xx[i+1][1];
	y_k1.z = z(i+1,y_k1.x,y_k1.y);
	y_k1.dx1 = zder(i+1,y_k1.x,y_k1.y,0);
	y_k1.dx2 = zder2_1(i+1,y_k1.x,y_k1.y,0);
	y_k1.dy1 = zder(i+1,y_k1.x,y_k1.y,1);
	y_k1.dy2 = zder2_2(i+1,y_k1.x,y_k1.y,1);
	y_k1.dxy2 = zder2_1(i+1,y_k1.x,y_k1.y,1);
	
	if (i!=nr3) {
		
		y_k1.gx1 = gx[i];
		y_k1.gy1 = gy[i];
		y_k1.gz1 = gz[i];
		y_k1.v1 = v[i]+y_k1.gx1*(y_k1.x-xref[i])+y_k1.gy1*(y_k1.y-yref[i])+y_k1.gz1*(y_k1.z-zref[i]);
		y_k1.gx2 = gx[i+1];
		y_k1.gy2 = gy[i+1];
		y_k1.gz2 = gz[i+1];
		y_k1.v2 = v[i+1]+y_k1.gx2*(y_k1.x-xref[i+1])+y_k1.gy2*(y_k1.y-yref[i+1])+y_k1.gz2*(y_k1.z-zref[i+1]);
		
	} else if (i==nr3) { /* For the air above at the last reflection*/
		
		y_k1.gx1 = gx[i];
		y_k1.gy1 = gy[i];
		y_k1.gz1 = gz[i];
		y_k1.v1 = v[i]+y_k1.gx1*(y_k1.x-xref[i])+y_k1.gy1*(y_k1.y-yref[i])+y_k1.gz1*(y_k1.z-zref[i]);	
		y_k1.gx2 = 0;
		y_k1.gz2 = 0;
		y_k1.v2 = 0;
	}	
}

/* Initialize one segment at a time for computing traveltime-------------------------------------*/

void half_initialize(int i /*Indicator of layer*/,
					 int nr3 /*number of reflection*/,
					 float **xx /*position*/, 
					 float *v /*velocity*/, 
					 float *xref /*reference position*/,
					 float *yref /*reference position*/,
					 float *zref /*reference position*/, 
					 float *gx /*x-gradient*/,
					 float *gy /*y-gradient*/,
					 float *gz /*z-gradient*/,
					 func1 z /*z(x,y)*/,
					 func4 zder /*z'(x,y)*/,
					 func4 zder2_1 /*d2z/dx2 and d2z/dxdy*/,
					 func4 zder2_2 /*d2z/dx2 and d2z/dxdy*/) 
/*<Half Initialize geometry>*/
{
	/* y_k----------------------------------------------------------------------------------------*/
	
	y_k.x = xx[i][0];
	y_k.y = xx[i][1];
	y_k.z = z(i,y_k.x,y_k.y);
	y_k.dx1 = zder(i,y_k.x,y_k.y,0);
	y_k.dx2 = zder2_1(i,y_k.x,y_k.y,0);
	y_k.dy1 = zder(i,y_k.x,y_k.y,1);
	y_k.dy2 = zder2_2(i,y_k.x,y_k.y,1);
	y_k.dxy2 = zder2_1(i,y_k.x,y_k.y,1);
	y_k.gx1 = 0;
	y_k.gy1 = 0;
	y_k.gz1 = 0;
	y_k.v1 = 0;
	y_k.gx2 = gx[i];
		y_k.gy2 = gy[i];
	y_k.gz2 = gz[i];
	y_k.v2 = v[i]+y_k.gx2*(y_k.x-xref[i])+y_k.gy2*(y_k.y-yref[i])+y_k.gz2*(y_k.z-zref[i]);
	
	/* y_k1 (y_k+1 th)-----------------------------------------------------------------------------*/
	
	y_k1.x = xx[i+1][0];
	y_k1.y = xx[i+1][1];
	y_k1.z = z(i+1,y_k1.x,y_k1.y);
	y_k1.dx1 = zder(i+1,y_k1.x,y_k1.y,0);
	y_k1.dx2 = zder2_1(i+1,y_k1.x,y_k1.y,0);
	y_k1.dy1 = zder(i+1,y_k1.x,y_k1.y,1);
	y_k1.dy2 = zder2_2(i+1,y_k1.x,y_k1.y,1);
	y_k1.dxy2 = zder2_1(i+1,y_k1.x,y_k1.y,1);
	y_k1.gx1 = gx[i];
	y_k1.gy1 = gy[i];
	y_k1.gz1 = gz[i];
	y_k1.v1 = v[i]+y_k1.gx1*(y_k1.x-xref[i])+y_k1.gy1*(y_k1.y-yref[i])+y_k1.gz1*(y_k1.z-zref[i]);
	y_k1.gx2 = 0;
	y_k1.gy2 = 0;
	y_k1.gz2 = 0;
	y_k1.v2 = 0;
	
}

/* T_hat functions------------------------------------------------------------------------------------------------------*/

double T_hat_k(func2 T_k)
/*<Traveltime>*/
{
	
	double t_k;
	
	t_k = T_k(y_k,y_k1);
	
	return t_k;
}

/* First Derivative-----------------------------------------------------------------------------------------------------*/

double T_hat_k_k_1(func2 T_k_k_1,func2 T_k_zk)
/*<Derivative of T_hat with respect to x_k>*/
{
	
	double t_k_k_1;
	
	t_k_k_1 = T_k_k_1(y_k,y_k1)+T_k_zk(y_k,y_k1)*y_k.dx1;
	
	return t_k_k_1;
}

double T_hat_k_k_2(func2 T_k_k_2,func2 T_k_zk)
/*<Derivative of T_hat with respect to y_k>*/
{
	
	double t_k_k_2;
	
	t_k_k_2 = T_k_k_2(y_k,y_k1)+T_k_zk(y_k,y_k1)*y_k.dy1;
	
	return t_k_k_2;
}

double T_hat_k_k1_1(func2 T_k_k1_1,func2 T_k_zk1)
/*<Derivative of T_hat with respect to x_k+1>*/
{
	
	double t_k_k1_1;
	
	t_k_k1_1 = T_k_k1_1(y_k,y_k1)+T_k_zk1(y_k,y_k1)*y_k1.dx1;
	
	return t_k_k1_1;
}

double T_hat_k_k1_2(func2 T_k_k1_2,func2 T_k_zk1)
/*<Derivative of T_hat with respect to y_k+1>*/
{
	
	double t_k_k1_2;
	
	t_k_k1_2 = T_k_k1_2(y_k,y_k1)+T_k_zk1(y_k,y_k1)*y_k1.dy1;
	
	return t_k_k1_2;
}

double T_hat_1k_k_1(func2 T_k_k1_1,func2 T_k_zk1)
/*<Derivative of T_hat_k-1th with respect to x_k>*/
{
	
	double t_1k_k_1;
	
	t_1k_k_1 = T_k_k1_1(y_1k,y_k)+T_k_zk1(y_1k,y_k)*y_k.dx1;
	
	return t_1k_k_1;
}

double T_hat_1k_k_2(func2 T_k_k1_2,func2 T_k_zk1)
/*<Derivative of T_hat_k-1th with respect to y_k>*/
{
	
	double t_1k_k_2;
	
	t_1k_k_2 = T_k_k1_2(y_1k,y_k)+T_k_zk1(y_1k,y_k)*y_k.dy1;
	
	return t_1k_k_2;
}

/* Second Derivative----------------------------------------------------------------------------------------------------*/

/* k_k_k Family--------------------------------------*/

double T_hat_k_k_k_1(func2 T_k_k_k_1,func2 T_k_k_zk_1,func2 T_k_zk,func2 T_k_zk_zk)
/*<Second derivative of T_hat with respect to x_k>*/
{
	
	double t_k_k_k_1;
	
	t_k_k_k_1 = T_k_k_k_1(y_k,y_k1)+2*T_k_k_zk_1(y_k,y_k1)*y_k.dx1+T_k_zk_zk(y_k,y_k1)*pow(y_k.dx1,2)+T_k_zk(y_k,y_k1)*y_k.dx2;
	
	return t_k_k_k_1;
}

double T_hat_k_k_k_2(func2 T_k_k_k_2,func2 T_k_k_zk_2,func2 T_k_zk,func2 T_k_zk_zk)
/*<Second derivative of T_hat with respect to y_k>*/
{
	
	double t_k_k_k_2;
	
	t_k_k_k_2 = T_k_k_k_2(y_k,y_k1)+2*T_k_k_zk_2(y_k,y_k1)*y_k.dy1+T_k_zk_zk(y_k,y_k1)*pow(y_k.dy1,2)+T_k_zk(y_k,y_k1)*y_k.dy2;
	
	return t_k_k_k_2;
}

double T_hat_k_k_k_12(func2 T_k_k_k_12,func2 T_k_k_zk_1,func2 T_k_k_zk_2,func2 T_k_zk,func2 T_k_zk_zk)
/*<Second derivative of T_hat with respect to x_k and y_k>*/
{
	
	double t_k_k_k_12;
	
	t_k_k_k_12 = T_k_k_k_12(y_k,y_k1)+T_k_k_zk_1(y_k,y_k1)*y_k.dy1+T_k_k_zk_2(y_k,y_k1)*y_k.dx1+T_k_zk_zk(y_k,y_k1)*y_k.dx1*y_k.dy1+T_k_zk(y_k,y_k1)*y_k.dxy2;
	
	return t_k_k_k_12;
}

/* k_k1_k1 Family--------------------------------------*/

double T_hat_k_k1_k1_1(func2 T_k_k1_k1_1,func2 T_k_k1_zk1_1,func2 T_k_zk1,func2 T_k_zk1_zk1)
/*<Second derivative of T_hat with respect to x_k+1>*/
{
	
	double t_k_k1_k1_1;
	
	t_k_k1_k1_1 = T_k_k1_k1_1(y_k,y_k1)+2*T_k_k1_zk1_1(y_k,y_k1)*y_k1.dx1+T_k_zk1_zk1(y_k,y_k1)*pow(y_k1.dx1,2)+T_k_zk1(y_k,y_k1)*y_k1.dx2;
	
	return t_k_k1_k1_1;
}

double T_hat_k_k1_k1_2(func2 T_k_k1_k1_2,func2 T_k_k1_zk1_2,func2 T_k_zk1,func2 T_k_zk1_zk1)
/*<Second derivative of T_hat with respect to y_k+1>*/
{
	
	double t_k_k1_k1_2;
	
	t_k_k1_k1_2 = T_k_k1_k1_2(y_k,y_k1)+2*T_k_k1_zk1_2(y_k,y_k1)*y_k1.dy1+T_k_zk1_zk1(y_k,y_k1)*pow(y_k1.dy1,2)+T_k_zk1(y_k,y_k1)*y_k1.dy2;
	
	return t_k_k1_k1_2;
}

double T_hat_k_k1_k1_12(func2 T_k_k1_k1_12,func2 T_k_k1_zk1_1,func2 T_k_k1_zk1_2,func2 T_k_zk1,func2 T_k_zk1_zk1)
/*<Second derivative of T_hat with respect to x_k+1 and y_k+1>*/
{
	
	double t_k_k1_k1_12;
	
	t_k_k1_k1_12 = T_k_k1_k1_12(y_k,y_k1)+T_k_k1_zk1_1(y_k,y_k1)*y_k1.dy1+T_k_k1_zk1_2(y_k,y_k1)*y_k1.dx1+T_k_zk1_zk1(y_k,y_k1)*y_k1.dx1*y_k1.dy1+T_k_zk1(y_k,y_k1)*y_k1.dxy2;
	
	return t_k_k1_k1_12;
}

/* 1k_k_k Family--------------------------------------*/

double T_hat_1k_k_k_1(func2 T_k_k1_k1_1,func2 T_k_k1_zk1_1,func2 T_k_zk1,func2 T_k_zk1_zk1)
/*<Second derivative of T_hat_k-1th with respect to x_k>*/
{
	
	double t_1k_k_k_1;
	
	t_1k_k_k_1 = T_k_k1_k1_1(y_1k,y_k)+2*T_k_k1_zk1_1(y_1k,y_k)*y_k.dx1+T_k_zk1_zk1(y_1k,y_k)*pow(y_k.dx1,2)+T_k_zk1(y_1k,y_k)*y_k.dx2;
	
	return t_1k_k_k_1;
}

double T_hat_1k_k_k_2(func2 T_k_k1_k1_2,func2 T_k_k1_zk1_2,func2 T_k_zk1,func2 T_k_zk1_zk1)
/*<Second derivative of T_hat_k-1th with respect to y_k>*/
{
	
	double t_1k_k_k_2;
	
	t_1k_k_k_2 = T_k_k1_k1_2(y_1k,y_k)+2*T_k_k1_zk1_2(y_1k,y_k)*y_k.dy1+T_k_zk1_zk1(y_1k,y_k)*pow(y_k.dy1,2)+T_k_zk1(y_1k,y_k)*y_k.dy2;
	
	return t_1k_k_k_2;
}

double T_hat_1k_k_k_12(func2 T_k_k1_k1_12,func2 T_k_k1_zk1_1,func2 T_k_k1_zk1_2, func2 T_k_zk1,func2 T_k_zk1_zk1)
/*<Second derivative of T_hat_k-1th with respect to x_k and y_k>*/
{
	
	double t_1k_k_k_12;
	
	t_1k_k_k_12 = T_k_k1_k1_12(y_1k,y_k)+T_k_k1_zk1_1(y_1k,y_k)*y_k.dy1+T_k_k1_zk1_2(y_1k,y_k)*y_k.dx1+T_k_zk1_zk1(y_1k,y_k)*y_k.dx1*y_k.dy1+T_k_zk1(y_1k,y_k)*y_k.dxy2;
	
	return t_1k_k_k_12;
}

/* k_k_k1 Family--------------------------------------*/

double T_hat_k_k_k1_1(func2 T_k_k_k1_1,func2 T_k_k1_zk_1,func2 T_k_k_zk1_1, func2 T_k_zk_zk1)
/*<Second derivative of T_hat with respect to x_k and x_k+1>*/
{
	
	double t_k_k_k1_1;
	
	t_k_k_k1_1 = T_k_k_k1_1(y_k,y_k1)+T_k_k_zk1_1(y_k,y_k1)*y_k1.dx1+T_k_k1_zk_1(y_k,y_k1)*y_k.dx1+T_k_zk_zk1(y_k,y_k1)*y_k1.dx1*y_k.dx1;
	
	return t_k_k_k1_1;
}

double T_hat_k_k_k1_2(func2 T_k_k_k1_2,func2 T_k_k1_zk_2,func2 T_k_k_zk1_2, func2 T_k_zk_zk1)
/*<Second derivative of T_hat with respect to y_k and y_k+1>*/
{
	
	double t_k_k_k1_2;
	
	t_k_k_k1_2 = T_k_k_k1_2(y_k,y_k1)+T_k_k_zk1_2(y_k,y_k1)*y_k1.dy1+T_k_k1_zk_2(y_k,y_k1)*y_k.dy1+T_k_zk_zk1(y_k,y_k1)*y_k1.dy1*y_k.dy1;
	
	return t_k_k_k1_2;
}

double T_hat_k_k_k1_12(func2 T_k_k_k1_12,func2 T_k_k1_zk_2,func2 T_k_k_zk1_1, func2 T_k_zk_zk1)
/*<Second derivative of T_hat with respect to x_k and y_k+1>*/
{
	
	double t_k_k_k1_12;
	
	t_k_k_k1_12 = T_k_k_k1_12(y_k,y_k1)+T_k_k_zk1_1(y_k,y_k1)*y_k1.dy1+T_k_k1_zk_2(y_k,y_k1)*y_k.dx1+T_k_zk_zk1(y_k,y_k1)*y_k1.dy1*y_k.dx1;
	
	return t_k_k_k1_12;
}

double T_hat_k_k_k1_21(func2 T_k_k_k1_21,func2 T_k_k1_zk_1,func2 T_k_k_zk1_2, func2 T_k_zk_zk1)
/*<Second derivative of T_hat with respect to y_k and x_k+1>*/
{
	
	double t_k_k_k1_21;
	
	t_k_k_k1_21 = T_k_k_k1_21(y_k,y_k1)+T_k_k_zk1_2(y_k,y_k1)*y_k1.dx1+T_k_k1_zk_1(y_k,y_k1)*y_k.dy1+T_k_zk_zk1(y_k,y_k1)*y_k1.dx1*y_k.dy1;
	
	return t_k_k_k1_21;
}

/* 1k_1k_k Family--------------------------------------*/

double T_hat_1k_1k_k_1(func2 T_k_k_k1_1,func2 T_k_k1_zk_1,func2 T_k_k_zk1_1, func2 T_k_zk_zk1)
/*<Second derivative of T_hat with respect to x_1k and x_k>*/
{
	
	double t_1k_1k_k_1;
	
	t_1k_1k_k_1 = T_k_k_k1_1(y_1k,y_k)+T_k_k_zk1_1(y_1k,y_k)*y_k.dx1+T_k_k1_zk_1(y_1k,y_k)*y_1k.dx1+T_k_zk_zk1(y_1k,y_k)*y_k.dx1*y_1k.dx1;
	
	return t_1k_1k_k_1;
}

double T_hat_1k_1k_k_2(func2 T_k_k_k1_2,func2 T_k_k1_zk_2,func2 T_k_k_zk1_2, func2 T_k_zk_zk1)
/*<Second derivative of T_hat with respect to y_k and y_k+1>*/
{
	
	double t_1k_1k_k_2;
	
	t_1k_1k_k_2 = T_k_k_k1_2(y_1k,y_k)+T_k_k_zk1_2(y_1k,y_k)*y_k.dy1+T_k_k1_zk_2(y_1k,y_k)*y_1k.dy1+T_k_zk_zk1(y_1k,y_k)*y_k.dy1*y_1k.dy1;
	
	return t_1k_1k_k_2;
}

double T_hat_1k_1k_k_12(func2 T_k_k_k1_12,func2 T_k_k1_zk_2,func2 T_k_k_zk1_1, func2 T_k_zk_zk1)
/*<Second derivative of T_hat with respect to x_k and y_k+1>*/
{
	
	double t_1k_1k_k_12;
	
	t_1k_1k_k_12 = T_k_k_k1_12(y_1k,y_k)+T_k_k_zk1_1(y_1k,y_k)*y_k.dy1+T_k_k1_zk_2(y_1k,y_k)*y_1k.dx1+T_k_zk_zk1(y_1k,y_k)*y_k.dy1*y_1k.dx1;
	
	return t_1k_1k_k_12;
}

double T_hat_1k_1k_k_21(func2 T_k_k_k1_21,func2 T_k_k1_zk_1,func2 T_k_k_zk1_2, func2 T_k_zk_zk1)
/*<Second derivative of T_hat with respect to y_k and x_k+1>*/
{
	
	double t_1k_1k_k_21;
	
	t_1k_1k_k_21 = T_k_k_k1_21(y_1k,y_k)+T_k_k_zk1_2(y_1k,y_k)*y_k.dx1+T_k_k1_zk_1(y_1k,y_k)*y_1k.dy1+T_k_zk_zk1(y_1k,y_k)*y_k.dx1*y_1k.dy1;
	
	return t_1k_1k_k_21;
}

