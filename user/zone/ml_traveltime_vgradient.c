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

#include "ml_traveltime_vgradient.h"

#include "general_traveltime.h"
/*^*/

#ifndef _ml_traveltime_vgradient_h

typedef struct twod {
	float x; /* x-coordinate*/
	float z; /* z-coordinate*/
	float d1; /* First derivative*/
	float d2; /* Second derivative*/
	float v1; /* Velocity at the reflector from above*/
	float v2; /* Velocity at the reflector from below*/
	float gx1;/* x-direction velocity gradient from above*/
	float gx2;/* x-direction velocity gradient from below*/
	float gz1;/* z-direction velocity gradient from above*/
	float gz2;/* z-direction velocity gradient from below*/
	float c111;/* c11 from above*/
	float c112;/* c11 from below*/
	float c331;/* c33 from above*/
	float c332;/* c33 from below*/
	float Q11; /* Q1 (anelliptic parameter) from above*/
	float Q12; /* Q1 (anelliptic parameter) from below*/
	float Q31; /* Q3 (anelliptic parameter) from above*/
	float Q32; /* Q3 (anelliptic parameter) from below*/
	float S11; /* S1 from above*/
	float S12; /* S1 from below*/
	float S31; /* S3 from above*/
	float S32; /* S3 from below*/
} twod;
/* Structure pointer */

#endif

static float eps = 0.0001; /*small constant to avoid division by zero*/

/* Traveltime functions for gradient velocity------------------------------------------------------------------------------*/

double v(twod y_k, twod x_ref)
{
	double v;
	
	v = x_ref.v2 + y_k.gx2*(y_k.x-x_ref.x)+y_k.gz2*(y_k.z-x_ref.z);
	y_k.v2 = v;
	
	return y_k.v2;
} 

double T1_k(twod y_k,twod y_k1)
/*<Traveltime>*/
{
	double t_k;
	
	t_k = (1/hypotf(y_k.gx2,y_k.gz2))*log((1+(pow(hypotf(y_k.gx2,y_k.gz2),2)*pow(hypotf(y_k1.x-y_k.x, y_k1.z-y_k.z),2))/(2*y_k.v2*y_k1.v1))+sqrt(pow(1+(pow(hypotf(y_k.gx2,y_k.gz2),2)*pow(hypotf(y_k1.x-y_k.x, y_k1.z-y_k.z),2))/(2*y_k.v2*y_k1.v1),2)-1));
	
	return t_k;
	
}

double T1_k_k(twod y_k, twod y_k1)
/*<Derivative of T with respect to x_k>*/
{
	double t_k_k,g0;
	
	g0 = hypotf(y_k.gx2,y_k.gz2);
	
	t_k_k = (sqrt(2)*g0*(y_k.x-y_k1.x))/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1)));
	
	return t_k_k;
	
}

double T1_k_k1(twod y_k, twod y_k1)
/*<Derivative of T with respect to x_k1>*/
{
	double t_k_k1,g0;
	
	g0 = hypotf(y_k.gx2,y_k.gz2);
	
	t_k_k1 = (sqrt(2)*g0*(y_k1.x-y_k.x))/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1)));
	
	return t_k_k1;
	
}

double T1_k_k_k(twod y_k, twod y_k1)  
/*<Second derivative of T with respect to x_k>*/
{
	double t_k_k_k,g0;
	
	g0 = hypotf(y_k.gx2,y_k.gz2);
	
	t_k_k_k = ((-1)*sqrt(2)*pow(g0,3)*pow(y_k.x-y_k1.x,2))/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*g0)/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1))) - (sqrt(2)*pow(g0,3)*pow(y_k.x-y_k1.x,2))/(eps+2*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1)));	
	
	return t_k_k_k;
	
}

double T1_k_k1_k1(twod y_k, twod y_k1) 
/*<Second derivative of T with respect to x_k1>*/
{
	double t_k_k1_k1,g0;
	
	g0 = hypotf(y_k.gx2,y_k.gz2);
	
	t_k_k1_k1 = ((-1)*sqrt(2)*pow(g0,3)*pow(y_k.x-y_k1.x,2))/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*g0)/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1))) - (sqrt(2)*pow(g0,3)*pow(y_k.x-y_k1.x,2))/(eps+2*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1)));		
	
	return t_k_k1_k1;
	
}

double T1_k_k_k1(twod y_k, twod y_k1) 
/*<Second derivative of T with respect to x_k and x_k1>*/
{
	double t_k_k_k1,g0;
	
	g0 = hypotf(y_k.gx2,y_k.gz2);
	
	t_k_k_k1 = (sqrt(2)*pow(g0,3)*pow(y_k.x-y_k1.x,2))/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1),3/2)) - (sqrt(2)*g0)/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1))) + (sqrt(2)*pow(g0,3)*pow(y_k.x-y_k1.x,2))/(eps+2*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1)));
	
	return t_k_k_k1;
	
}

double T1_k_zk(twod y_k, twod y_k1) 
/*<Derivative of T with respect to z_k>*/
{
	double t_k_zk,g0;
	
	g0 = hypotf(y_k.gx2,y_k.gz2);
	
	t_k_zk = (sqrt(2)*g0*(y_k.z-y_k1.z))/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1)));
	
	return t_k_zk;
	
}

double T1_k_zk1(twod y_k, twod y_k1) 
/*<Derivative of T with respect to z_k1>*/
{
	double t_k_zk1,g0;
	
	g0 = hypotf(y_k.gx2,y_k.gz2);
	
	t_k_zk1 = (sqrt(2)*g0*(y_k1.z-y_k.z))/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1)));
	
	return t_k_zk1;
	
}

double T1_k_zk_zk(twod y_k, twod y_k1)  
/*<Second Derivative of T with respect to z_k>*/
{
	double t_k_zk_zk,g0;
	
	g0 = hypotf(y_k.gx2,y_k.gz2);
	
	t_k_zk_zk = ((-1)*sqrt(2)*pow(g0,3)*pow(y_k.z-y_k1.z,2))/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*g0)/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1))) - (sqrt(2)*pow(g0,3)*pow(y_k.z-y_k1.z,2))/(eps+2*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1)));	
	
	return t_k_zk_zk;
	
}

double T1_k_zk1_zk1(twod y_k, twod y_k1)  
/*<Second Derivative of T with respect to z_k1>*/
{
	double t_k_zk1_zk1,g0;
	
	g0 = hypotf(y_k.gx2,y_k.gz2);
	
	t_k_zk1_zk1 = ((-1)*sqrt(2)*pow(g0,3)*pow(y_k.z-y_k1.z,2))/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*g0)/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1))) - (sqrt(2)*pow(g0,3)*pow(y_k.z-y_k1.z,2))/(eps+2*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1)));	
	
	return t_k_zk1_zk1;
	
}

double T1_k_zk_zk1(twod y_k, twod y_k1)  
/*<Second Derivative of T with respect to z_k and z_k1>*/
{
	double t_k_zk_zk1,g0;
	
	g0 = hypotf(y_k.gx2,y_k.gz2);
	
	t_k_zk_zk1 = (sqrt(2)*pow(g0,3)*pow(y_k.z-y_k1.z,2))/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1),3/2)) - (sqrt(2)*g0)/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1))) + (sqrt(2)*pow(g0,3)*pow(y_k.z-y_k1.z,2))/(eps+2*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1)));	
	
	return t_k_zk_zk1;
	
}

double T1_k_k_zk(twod y_k, twod y_k1) 
/*<Second derivative of T with respect to x_k and z_k>*/
{
	double t_k_k_zk,g0;
	
	g0 = hypotf(y_k.gx2,y_k.gz2);
	
	t_k_k_zk = (sqrt(2)*pow(g0,3)*(y_k1.z-y_k.z)*(y_k.x-y_k1.x))/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*pow(g0,3)*(y_k1.z-y_k.z)*(y_k.x-y_k1.x))/(eps+2*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1)));		
	
	return t_k_k_zk;
	
}

double T1_k_k1_zk1(twod y_k, twod y_k1) 
/*<Second derivative of T with respect to x_k1 and z_k1>*/
{
	double t_k_k1_zk1,g0;
	
	g0 = hypotf(y_k.gx2,y_k.gz2);
	
	t_k_k1_zk1 = (sqrt(2)*pow(g0,3)*(y_k1.z-y_k.z)*(y_k.x-y_k1.x))/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*pow(g0,3)*(y_k1.z-y_k.z)*(y_k.x-y_k1.x))/(eps+2*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1)));		
	
	return t_k_k1_zk1;
	
}

double T1_k_k_zk1(twod y_k, twod y_k1)  
/*<Second derivative of T with respect to x_k and z_k1>*/
{
	double t_k_k_zk1,g0;
	
	g0 = hypotf(y_k.gx2,y_k.gz2);
	
	
	t_k_k_zk1 = ((sqrt(2)*pow(g0,3)*(y_k.z-y_k1.z)*(y_k.x-y_k1.x))/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*pow(g0,3)*(y_k.z-y_k1.z)*(y_k.x-y_k1.x))/(eps+2*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1))));	
	
	return t_k_k_zk1;
	
}

double T1_k_k1_zk(twod y_k, twod y_k1)  
/*<Second derivative of T with respect to x_k1 and z_k>*/
{
	double t_k_k1_zk,g0;
	
	g0 = hypotf(y_k.gx2,y_k.gz2);
	
	t_k_k1_zk = ((sqrt(2)*pow(g0,3)*(y_k.z-y_k1.z)*(y_k.x-y_k1.x))/(eps+sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*pow(g0,3)*(y_k.z-y_k1.z)*(y_k.x-y_k1.x))/(eps+2*pow(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*(pow(y_k.z-y_k1.z,2)+pow(y_k.x-y_k1.x,2))/(y_k.v2*y_k1.v1))));	
	
	return t_k_k1_zk;
	
}

