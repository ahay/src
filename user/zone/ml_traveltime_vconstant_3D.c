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

#include "ml_traveltime_vconstant_3D.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "general_traveltime_3D.h"
/*^*/

#ifndef _ml_traveltime_vconstant_3D_h

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

#endif

static float eps = 0.0001; /*small constant to avoid division by zero*/

/*Traveltime functions for constant velocity---------------------------------------------------------------------------*/
double T0_k(threed y_k,threed y_k1)
/*<Traveltime>*/
{
	double t_k;
	
	t_k = pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),0.5)/y_k.v2;
	
	return t_k;
	
}

/* First Derivative (6)---------------------------------------------------------------------------------------*/

double T0_k_k_1(threed y_k, threed y_k1) 
/*<Derivative of T with respect to x_k>*/
{
	double t_k_k_1;
	
	t_k_k_1 = (y_k.x-y_k1.x)/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),0.5));
	
	return t_k_k_1;
	
}

double T0_k_k_2(threed y_k, threed y_k1) 
/*<Derivative of T with respect to y_k>*/
{
	double t_k_k_2;
	
	t_k_k_2 = (y_k.y-y_k1.y)/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),0.5));
	
	return t_k_k_2;
	
}

double T0_k_k1_1(threed y_k, threed y_k1) 
/*<Derivative of T with respect to x_k1>*/
{
	double t_k_k1_1;
	
	t_k_k1_1 = (y_k1.x-y_k.x)/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),0.5));
	
	return t_k_k1_1;
	
}

double T0_k_k1_2(threed y_k, threed y_k1) 
/*<Derivative of T with respect to y_k1>*/
{
	double t_k_k1_2;
	
	t_k_k1_2 = (y_k1.y-y_k.y)/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),0.5));
	
	return t_k_k1_2;
	
}

double T0_k_zk(threed y_k, threed y_k1)  
/*<Derivative of T with respect to z_k>*/
{
	double t_k_zk;
	
	
	t_k_zk = (y_k.z-y_k1.z)/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),0.5));
	
	return t_k_zk;
	
}

double T0_k_zk1(threed y_k, threed y_k1) 
/*<Derivative of T with respect to z_k1>*/
{
	double t_k_zk1;
	
	
	t_k_zk1 = (y_k1.z-y_k.z)/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),0.5));
	
	return t_k_zk1;
	
}



/* Second Derivative------------------------------------------------------------------------------------------*/

/* k_k_k Family (6)-----------------------------------*/

double T0_k_k_k_1(threed y_k, threed y_k1)
/*<Second derivative of T with respect to x_k>*/
{
	double t_k_k_k_1;
	
	t_k_k_k_1 = 1/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),0.5)) - (pow(y_k.x-y_k1.x,2))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5));
	
	return t_k_k_k_1;
	
}

double T0_k_k_k_2(threed y_k, threed y_k1)
/*<Second derivative of T with respect to y_k>*/
{
	double t_k_k_k_2;
	
	t_k_k_k_2 = 1/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),0.5)) - (pow(y_k.y-y_k1.y,2))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5));
	
	return t_k_k_k_2;
	
}

double T0_k_k_k_12(threed y_k, threed y_k1)
/*<Second derivative of T with respect to x_k and y_k>*/
{
	double t_k_k_k_12;
	
	t_k_k_k_12 = ((y_k1.y-y_k.y)*(y_k.x-y_k1.x))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5));
	
	return t_k_k_k_12;
	
}

double T0_k_k_zk_1(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k and z_k>*/
{
	double t_k_k_zk_1;
	
	t_k_k_zk_1 =  ((y_k1.z-y_k.z)*(y_k.x-y_k1.x))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5));	
	
	return t_k_k_zk_1;
	
}

double T0_k_k_zk_2(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to y_k and y_k>*/
{
	double t_k_k_zk_2;
	
	t_k_k_zk_2 =  ((y_k1.z-y_k.z)*(y_k.y-y_k1.y))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5));	
	
	return t_k_k_zk_2;
	
}

double T0_k_zk_zk(threed y_k, threed y_k1)  
/*<Second Derivative of T with respect to z_k>*/
{
	double t_k_zk_zk;
	
	t_k_zk_zk = 1/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),0.5)) - (pow(y_k.z-y_k1.z,2))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5));
	
	return t_k_zk_zk;
	
}

/* k_k1_k1 Family & 1k_k_k Family (6)----------------------------------*/

double T0_k_k1_k1_1(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k1>*/
{
	double t_k_k1_k1_1;
	
	t_k_k1_k1_1 = 1/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),0.5)) - (pow(y_k.x-y_k1.x,2))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5));
	
	return t_k_k1_k1_1;
	
}

double T0_k_k1_k1_2(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to y_k1>*/
{
	double t_k_k1_k1_2;
	
	t_k_k1_k1_2 = 1/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),0.5)) - (pow(y_k.y-y_k1.y,2))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5));
	
	return t_k_k1_k1_2;
	
}

double T0_k_k1_k1_12(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k1 and y_k1>*/
{
	double t_k_k1_k1_12;
	
	t_k_k1_k1_12 = ((y_k1.y-y_k.y)*(y_k.x-y_k1.x))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5));
	
	return t_k_k1_k1_12;
	
}

double T0_k_k1_zk1_1(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to x_k1 and z_k1>*/
{
	double t_k_k1_zk1_1;
	
	t_k_k1_zk1_1 = ((y_k1.z-y_k.z)*(y_k.x-y_k1.x))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5));	
	
	return t_k_k1_zk1_1;
	
}

double T0_k_k1_zk1_2(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to y_k1 and z_k1>*/
{
	double t_k_k1_zk1_2;
	
	t_k_k1_zk1_2 = ((y_k1.z-y_k.z)*(y_k.y-y_k1.y))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5));	
	
	return t_k_k1_zk1_2;
	
}

double T0_k_zk1_zk1(threed y_k, threed y_k1)  
/*<Second Derivative of T with respect to z_k1>*/
{
	double t_k_zk1_zk1;
	
	t_k_zk1_zk1 = 1/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),0.5)) - (pow(y_k.z-y_k1.z,2))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5));
	
	return t_k_zk1_zk1;
	
}

/* k_k_k1 Family & 1k_1k_k Family (9)----------------------------------*/

double T0_k_k_k1_1(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to x_k and x_k1>*/
{
	double t_k_k_k1_1;
	
	t_k_k_k1_1 = (-1)*(1/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),0.5)) - (pow(y_k.x-y_k1.x,2))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5)));
	
	return t_k_k_k1_1;
	
}

double T0_k_k_k1_2(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to x_k and x_k1>*/
{
	double t_k_k_k1_2;
	
	t_k_k_k1_2 = (-1)*(1/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),0.5)) - (pow(y_k.y-y_k1.y,2))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5)));
	
	return t_k_k_k1_2;
	
}

double T0_k_k_k1_12(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to x_k and x_k1>*/
{
	double t_k_k_k1_12;
	
	t_k_k_k1_12 = ((y_k.x-y_k1.x)*(y_k.y-y_k1.y))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5));	
	
	return t_k_k_k1_12;
	
}

double T0_k_k_k1_21(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to x_k and x_k1>*/
{
	double t_k_k_k1_21;
	
	t_k_k_k1_21 = ((y_k.x-y_k1.x)*(y_k.y-y_k1.y))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5));	
	
	return t_k_k_k1_21;
	
}

double T0_k_k1_zk_1(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k1 and z_k>*/
{
	double t_k_k1_zk_1;
	
	t_k_k1_zk_1 = ((y_k.z-y_k1.z)*(y_k.x-y_k1.x))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5));	
	
	return t_k_k1_zk_1;
	
}

double T0_k_k1_zk_2(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k1 and z_k>*/
{
	double t_k_k1_zk_2;
	
	t_k_k1_zk_2 = ((y_k.z-y_k1.z)*(y_k.y-y_k1.y))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5));	
	
	return t_k_k1_zk_2;
	
}

double T0_k_k_zk1_1(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k and z_k1>*/
{
	double t_k_k_zk1_1;
	
	t_k_k_zk1_1 = ((y_k.z-y_k1.z)*(y_k.x-y_k1.x))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5));	
	
	return t_k_k_zk1_1;
	
}

double T0_k_k_zk1_2(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k and z_k1>*/
{
	double t_k_k_zk1_2;
	
	t_k_k_zk1_2 = ((y_k.z-y_k1.z)*(y_k.y-y_k1.y))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5));	
	
	return t_k_k_zk1_2;
	
}

double T0_k_zk_zk1(threed y_k, threed y_k1)  
/*<Second Derivative of T with respect to z_k and z_k1>*/
{
	double t_k_zk_zk1;
	
	t_k_zk_zk1 = (pow(y_k.z-y_k1.z,2))/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),1.5)) -  1/(eps+y_k.v2*pow(pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2),0.5));
	
	return t_k_zk_zk1;
	
}

