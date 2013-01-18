/*
 *  general_traveltime.c
 *  Multi-Layered
 *
 *  Created by Yanadet Sripanich on 10/11/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "general_traveltime.h"

#ifndef _general_traveltime_h

typedef struct twod {
	float x; /*x-coordinate*/
	float z; /*z-coordinate*/
	float d1; /*First derivative*/
	float d2; /*Second derivative*/
	float v; /*velocity at the reflector*/
	float gx;/*x-direction velocity gradient*/
	float gz;/*z-direction velocity gradient*/
} twod;
/* structure pointer */
/*^*/

typedef float (*func1)(int, float);
/* function pointer for int and float -> float */
/*^*/

typedef float (*func2)(twod,twod);
/* function pointer for twod,twod -> float */
/*^*/

typedef struct func3 {
	func2 T_k;
	func2 T_k_k;
	func2 T_k_k1;
	func2 T_k_k_k;
	func2 T_k_k1_k1;
	func2 T_k_k_k1;
	func2 T_k_zk;
	func2 T_k_zk1;
	func2 T_k_zk_zk;
	func2 T_k_zk1_zk1;
	func2 T_k_zk_zk1;
	func2 T_k_k_zk;
	func2 T_k_k1_zk1;
	func2 T_k_k_zk1;
	func2 T_k_k1_zk;
} func3;
/* structure pointer */
/*^*/

#endif

static twod y_1k,y_k,y_k1;
/*static func3 f;*/

/*T(x_k,x_k+1,z_k,z_k+1) = T_hat(x_k,x_k+1)*/

/*Initialize two segments at a time*/
void initialize(int i /*Indicator of layer*/,
				float x[] /*position*/, 
				float v[] /*velocity*/, 
				float xref[] /*reference position*/, 
				float zref[] /*reference position*/, 
				float gx[] /*x-gradient*/,
				float gz[] /*z-gradient*/,
				func1 z /*z(x)*/,
				func1 zder /*z'(x)*/,
				func1 zder2 /*z''(x)*/)
/*<Initialize geometry>*/
{
	y_1k.x = x[i-1];
	y_1k.z = z(i-1,y_1k.x);
	y_1k.d1 = zder(i-1,y_1k.x);
	y_1k.d2 = zder2(i-1,y_1k.x);
	y_1k.gx = gx[i-1];
	y_1k.gz = gz[i-1];
	y_1k.v = v[i-1]+y_1k.gx*(y_1k.x-xref[i-1])+y_1k.gz*(y_1k.z-zref[i-1]);
	
	y_k.x = x[i];
	y_k.z = z(i,y_k.x);
	y_k.d1 = zder(i,y_k.x);
	y_k.d2 = zder2(i,y_k.x);
	y_k.gx = gx[i];
	y_k.gz = gz[i];
	y_k.v = v[i]+y_k.gx*(y_k.x-xref[i])+y_k.gz*(y_k.z-zref[i]);
	
	y_k1.x = x[i+1];
	y_k1.z = z(i+1,y_k1.x);
	y_k1.d1 = zder(i+1,y_k1.x);
	y_k1.d2 = zder2(i+1,y_k1.x);
	y_k1.gx = gx[i+1];
	y_k1.gz = gz[i+1];
	y_k1.v = v[i+1]+y_k1.gx*(y_k1.x-xref[i+1])+y_k1.gz*(y_k1.z-zref[i+1]);

}

/* Initialize one segment at a time for computing traveltime*/
void half_initialize(int i /*Indicator of layer*/,
				float x[] /*position*/, 
				float v[] /*velocity*/, 
				float xref[] /*reference position*/, 
				float zref[] /*reference position*/, 
				float gx[] /*x-gradient*/,
				float gz[] /*z-gradient*/,
				func1 z /*z(x)*/,
				func1 zder /*z'(x)*/,
				func1 zder2 /*z''(x)*/)
/*<Half Initialize geometry>*/
{
	
	y_k.x = x[i];
	y_k.z = z(i,y_k.x);
	y_k.d1 = zder(i,y_k.x);
	y_k.d2 = zder2(i,y_k.x);
	y_k.gx = gx[i];
	y_k.gz = gz[i];
	y_k.v = v[i]+y_k.gx*(y_k.x-xref[i])+y_k.gz*(y_k.z-zref[i]);
	
	y_k1.x = x[i+1];
	y_k1.z = z(i+1,y_k1.x);
	y_k1.d1 = zder(i+1,y_k1.x);
	y_k1.d2 = zder2(i+1,y_k1.x);
	y_k1.gx = gx[i+1];
	y_k1.gz = gz[i+1];
	y_k1.v = v[i+1]+y_k1.gx*(y_k1.x-xref[i+1])+y_k1.gz*(y_k1.z-zref[i+1]);
	
}


float T_hat_k(func2 T_k)
/*<Traveltime>*/
{
	
	float t_k;
	
	t_k = T_k(y_k,y_k1);
	
	return t_k;
}




float T_hat_k_k(func2 T_k_k,func2 T_k_zk)
/*<Derivative of T_hat with respect to x_k>*/
{
	
	float t_k_k;
	
	t_k_k = T_k_k(y_k,y_k1)+T_k_zk(y_k,y_k1)*y_k.d1;
	
	return t_k_k;
}


float T_hat_k_k1(func2 T_k_k1,func2 T_k_zk1)
/*<Derivative of T_hat with respect to x_k+1>*/
{
	
	float t_k_k1;
	
	t_k_k1 = T_k_k1(y_k,y_k1)+T_k_zk1(y_k,y_k1)*y_k1.d1;
	
	return t_k_k1;
}

float T_hat_1k_k(func2 T_k_k1,func2 T_k_zk1)
/*<Derivative of T_hat_k-1th with respect to x_k>*/
{
	
	float t_1k_k;
	
	t_1k_k = T_k_k1(y_1k,y_k)+T_k_zk1(y_1k,y_k)*y_k.d1;
	
	return t_1k_k;
}


float T_hat_k_k_k(func2 T_k_k_k,func2 T_k_k_zk,func2 T_k_zk,func2 T_k_zk_zk)
/*<Second derivative of T_hat with respect to x_k>*/
{
	
	float t_k_k_k;
	
	t_k_k_k = T_k_k_k(y_k,y_k1)+2*T_k_k_zk(y_k,y_k1)*y_k.d1+T_k_zk_zk(y_k,y_k1)*pow(y_k.d1,2)+T_k_zk(y_k,y_k1)*y_k.d2;
	
	return t_k_k_k;
}


float T_hat_k_k1_k1(func2 T_k_k1_k1,func2 T_k_k1_zk1,func2 T_k_zk1,func2 T_k_zk1_zk1)
/*<Second derivative of T_hat with respect to x_k+1>*/
{
	
	float t_k_k1_k1;
	
	t_k_k1_k1 = T_k_k1_k1(y_k,y_k1)+2*T_k_k1_zk1(y_k,y_k1)*y_k1.d1+T_k_zk1_zk1(y_k,y_k1)*pow(y_k1.d1,2)+T_k_zk1(y_k,y_k1)*y_k1.d2;
	
	return t_k_k1_k1;
}

float T_hat_1k_k_k(func2 T_k_k1_k1,func2 T_k_k1_zk1,func2 T_k_zk1,func2 T_k_zk1_zk1)
/*<Second derivative of T_hat_k-1th with respect to x_k>*/
{
	
	float t_1k_k_k;
	
	t_1k_k_k = T_k_k1_k1(y_1k,y_k)+2*T_k_k1_zk1(y_1k,y_k)*y_k.d1+T_k_zk1_zk1(y_1k,y_k)*pow(y_k.d1,2)+T_k_zk1(y_1k,y_k)*y_k.d2;
	
	return t_1k_k_k;
}

float T_hat_k_k_k1(func2 T_k_k_k1,func2 T_k_k1_zk,func2 T_k_k_zk1, func2 T_k_zk_zk1)
/*<Second derivative of T_hat with respect to x_k and x_k+1>*/
{
	
	float t_k_k_k1;
	
	t_k_k_k1 = T_k_k_k1(y_k,y_k1)+T_k_k_zk1(y_k,y_k1)*y_k1.d1+T_k_k1_zk(y_k,y_k1)*y_k.d1+T_k_zk_zk1(y_k,y_k1)*y_k1.d1*y_k.d1;
	
	return t_k_k_k1;
}

float T_hat_1k_1k_k(func2 T_k_k_k1,func2 T_k_k1_zk,func2 T_k_k_zk1, func2 T_k_zk_zk1)
/*<Second derivative of T_hat with respect to x_1k and x_k>*/
{
	
	float t_1k_1k_k;
	
	t_1k_1k_k = T_k_k_k1(y_1k,y_k)+T_k_k_zk1(y_1k,y_k)*y_k.d1+T_k_k1_zk(y_1k,y_k)*y_1k.d1+T_k_zk_zk1(y_1k,y_k)*y_k.d1*y_1k.d1;
	
	return t_1k_1k_k;
}


/*	int i
 
 func2 functions[] =
 {
 f.T_k,f.T_k_k,f.T_k_k1,f.T_k_k_k,f.T_k_k1_k1,f.T_k_k_k1,
 f.T_k_zk,f.T_k_zk1,f.T_k_zk_zk,f.T_k_zk1_zk1,f.T_k_zk_zk1,
 f.T_k_k_zk,f.T_k_k1_zk1,f.T_k_k_zk1,f.T_k_k1_zk
 };
 
 func2 functions0[] =
 {
 T0_k,T0_k_k,T0_k_k1,T0_k_k_k,T0_k_k1_k1,T0_k_k_k1,
 T0_k_zk,T0_k_zk1,T0_k_zk_zk,T0_k_zk1_zk1,T0_k_zk_zk1,
 T0_k_k_zk,T0_k_k1_zk1,T0_k_k_zk1,T0_k_k1_zk
 };
 
 func2 functions1[] =
 {
 T1_k,T1_k_k,T1_k_k1,T1_k_k_k,T1_k_k1_k1,T1_k_k_k1,
 T1_k_zk,T1_k_zk1,T1_k_zk_zk,T1_k_zk1_zk1,T1_k_zk_zk1,
 T1_k_k_zk,T1_k_k1_zk1,T1_k_k_zk1,T1_k_k1_zk
 };*/


