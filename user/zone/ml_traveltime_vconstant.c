/*
 *  ml_traveltime_vconstant.c
 *  Test_struct
 *
 *  Created by Yanadet Sripanich on 11/10/12.
 * 
 *
 */
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

#include "ml_traveltime_vconstant.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "general_traveltime.h"
/*^*/

#ifndef _ml_traveltime_vconstant_h

typedef struct twod {
	float x; /* x-coordinate*/
	float z; /* z-coordinate*/
	float d1; /* First derivative*/
	float d2; /* Second derivative*/
	float v; /* Velocity at the reflector*/
	float gx;/* x-direction velocity gradient*/
	float gz;/* z-diection velocity gradient*/
} twod;
/* Structure pointer */

#endif



/*Traveltime functions for constant velocity---------------------------------------------------------------------------*/
float T0_k(twod y_k,twod y_k1)
/*<Traveltime>*/
{
	float t_k;
	
	t_k = hypotf(y_k.x-y_k1.x,y_k.z-y_k1.z)/y_k.v2;
	
	return t_k;
	
}

float T0_k_k(twod y_k, twod y_k1) 
/*<Derivative of T with respect to x_k>*/
{
	float t_k_k;
	
	t_k_k = (y_k.x-y_k1.x)/(y_k.v2*hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x));
	
	return t_k_k;
	
}

float T0_k_k1(twod y_k, twod y_k1) 
/*<Derivative of T with respect to x_k1>*/
{
	float t_k_k1;
	
	t_k_k1 = (y_k1.x-y_k.x)/(y_k.v2*hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x));
	
	return t_k_k1;
	
}

float T0_k_k_k(twod y_k, twod y_k1)
/*<Second derivative of T with respect to x_k>*/
{
	float t_k_k_k;
	
	t_k_k_k = 1/(y_k.v2*hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x)) - (pow(y_k.x-y_k1.x,2))/(y_k.v2*pow(hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x),3));
	
	return t_k_k_k;
	
}

float T0_k_k1_k1(twod y_k, twod y_k1)  
/*<Second derivative of T with respect to x_k1>*/
{
	float t_k_k1_k1;
	
	t_k_k1_k1 = 1/(y_k.v2*hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x)) - (pow(y_k.x-y_k1.x,2))/(y_k.v2*pow(hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x),3));
				 
	return t_k_k1_k1;
	
}

float T0_k_k_k1(twod y_k, twod y_k1) 
/*<Second derivative of T with respect to x_k and x_k1>*/
{
	float t_k_k_k1;
	
	t_k_k_k1 = (-1)*(1/(y_k.v2*hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x)) - (pow(y_k.x-y_k1.x,2))/(y_k.v2*pow(hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x),3)));
					 
	return t_k_k_k1;
	
}

float T0_k_zk(twod y_k, twod y_k1)  
/*<Derivative of T with respect to z_k>*/
{
	float t_k_zk;
	
	t_k_zk = (y_k.z-y_k1.z)/(y_k.v2*hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x));
	
	return t_k_zk;
	
}

float T0_k_zk1(twod y_k, twod y_k1) 
/*<Derivative of T with respect to z_k1>*/
{
	float t_k_zk1;
	
	t_k_zk1 = (y_k1.z-y_k.z)/(y_k.v2*hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x));
	
	return t_k_zk1;
	
}

float T0_k_zk_zk(twod y_k, twod y_k1)  
/*<Second Derivative of T with respect to z_k>*/
{
	float t_k_zk_zk;
	
	t_k_zk_zk = 1/(y_k.v2*hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x)) - (pow(y_k.z-y_k1.z,2))/(y_k.v2*pow(hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x),3));
	
	return t_k_zk_zk;
	
}

float T0_k_zk1_zk1(twod y_k, twod y_k1)  
/*<Second Derivative of T with respect to z_k1>*/
{
	float t_k_zk1_zk1;
	
	t_k_zk1_zk1 = 1/(y_k.v2*hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x)) - (pow(y_k.z-y_k1.z,2))/(y_k.v2*pow(hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x),3));
	
	return t_k_zk1_zk1;
	
}

float T0_k_zk_zk1(twod y_k, twod y_k1)  
/*<Second Derivative of T with respect to z_k and z_k1>*/
{
	float t_k_zk_zk1;
	
	t_k_zk_zk1 = (pow(y_k.z-y_k1.z,2))/(y_k.v2*pow(hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x),3)) -  1/(y_k.v2*hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x));
	
	return t_k_zk_zk1;
	
}

float T0_k_k_zk(twod y_k, twod y_k1)  
/*<Second derivative of T with respect to x_k and z_k>*/
{
	float t_k_k_zk;
	
	t_k_k_zk =  ((y_k1.z-y_k.z)*(y_k.x-y_k1.x))/(y_k.v2*pow(hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x),3));	
	
	return t_k_k_zk;
	
}

float T0_k_k1_zk1(twod y_k, twod y_k1) 
/*<Second derivative of T with respect to x_k1 and z_k1>*/
{
	float t_k_k1_zk1;
	
	t_k_k1_zk1 = ((y_k1.z-y_k.z)*(y_k.x-y_k1.x))/(y_k.v2*pow(hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x),3));	
	
	return t_k_k1_zk1;
	
}

float T0_k_k_zk1(twod y_k, twod y_k1)  
/*<Second derivative of T with respect to x_k and z_k1>*/
{
	float t_k_k_zk1;
	
	t_k_k_zk1 = (-1)*(((y_k1.z-y_k.z)*(y_k.x-y_k1.x))/(y_k.v2*pow(hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x),3)));	
	
	return t_k_k_zk1;
	
}

float T0_k_k1_zk(twod y_k, twod y_k1)  
/*<Second derivative of T with respect to x_k1 and z_k>*/
{
	float t_k_k1_zk;
	
	t_k_k1_zk = (-1)*(((y_k1.z-y_k.z)*(y_k.x-y_k1.x))/(y_k.v2*pow(hypotf(y_k.z-y_k1.z,y_k.x-y_k1.x),3)));	
	
	return t_k_k1_zk;
	
}
