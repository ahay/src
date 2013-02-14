/*
 *  ml_traveltime_vgradient_3D.c
 *  
 *
 *  Created by Yanadet Sripanich on 2/6/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "ml_traveltime_vgradient_3D.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "general_traveltime_3D.h"
/*^*/

#ifndef _ml_traveltime_vgradient_3D_h

typedef struct threed {
	float x; /*x-coordinate*/
	float y; /*y-coordinate*/
	float z; /*z-coordinate*/
	float dx1; /*First derivative in x-direction*/
	float dx2; /*Second derivative in x-direction*/
	float dy1; /*First derivative in y-direction*/
	float dy2; /*Second derivative in y-direction*/
	float dxy2; /*Cross-derivative of x and y*/
	float v1; /*velocity at the reflector from above*/
	float v2; /*velocity at the reflector from below*/
	float gx1;/*x-direction velocity gradient from above*/
	float gx2;/*x-direction velocity gradient from below*/
	float gy1;/*y-direction velocity gradient from above*/
	float gy2;/*y-direction velocity gradient from below*/
	float gz1;/*z-direction velocity gradient from above*/
	float gz2;/*z-direction velocity gradient from below*/
} threed;
/* structure pointer */

#endif


float T1_k(threed y_k,threed y_k1)
/*<Traveltime>*/
{
	float t_k,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k = (1/pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5))*log((1+((pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2))*q)/(2*y_k.v2*y_k1.v1))+sqrt(pow(1+((pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2))*q)/(2*y_k.v2*y_k1.v1),2)-1));
	
	return t_k;
	
}

/*First Derivative (6)---------------------------------------------------------------------------------------*/

float T1_k_k_1(threed y_k, threed y_k1) 
/*<Derivative of T with respect to x_k>*/
{
	float t_k_k_1,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k_1 = (sqrt(2)*g0*(y_k.x-y_k1.x))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));
	
	return t_k_k_1;
	
}

float T1_k_k_2(threed y_k, threed y_k1) 
/*<Derivative of T with respect to y_k>*/
{
	float t_k_k_2,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k_2 = (sqrt(2)*g0*(y_k.y-y_k1.y))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));
	
	return t_k_k_2;
	
}

float T1_k_k1_1(threed y_k, threed y_k1) 
/*<Derivative of T with respect to x_k1>*/
{
	float t_k_k1_1,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k1_1 = (sqrt(2)*g0*(y_k1.x-y_k.x))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));
	
	return t_k_k1_1;
	
}

float T1_k_k1_2(threed y_k, threed y_k1) 
/*<Derivative of T with respect to y_k1>*/
{
	float t_k_k1_2,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k1_2 = (sqrt(2)*g0*(y_k1.y-y_k.y))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));
	
	return t_k_k1_2;
	
}

float T1_k_zk(threed y_k, threed y_k1)  
/*<Derivative of T with respect to z_k>*/
{
	float t_k_zk,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_zk = (sqrt(2)*g0*(y_k.z-y_k1.z))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));
	
	return t_k_zk;
	
}

float T1_k_zk1(threed y_k, threed y_k1) 
/*<Derivative of T with respect to z_k1>*/
{
	float t_k_zk1,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_zk1 =  (sqrt(2)*g0*(y_k1.z-y_k.z))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));
	
	return t_k_zk1;
	
}



/*Second Derivative------------------------------------------------------------------------------------------*/

/*k_k_k Family (6)-----------------------------------*/

float T1_k_k_k_1(threed y_k, threed y_k1)
/*<Second derivative of T with respect to x_k>*/
{
	float t_k_k_k_1,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k_k_1 = ((-1)*sqrt(2)*pow(g0,3)*pow(y_k.x-y_k1.x,2))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*g0)/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1))) - (sqrt(2)*pow(g0,3)*pow(y_k.x-y_k1.x,2))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));	
	
	return t_k_k_k_1;
	
}

float T1_k_k_k_2(threed y_k, threed y_k1)
/*<Second derivative of T with respect to y_k>*/
{
	float t_k_k_k_2,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k_k_2 = ((-1)*sqrt(2)*pow(g0,3)*pow(y_k.y-y_k1.y,2))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*g0)/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1))) - (sqrt(2)*pow(g0,3)*pow(y_k.y-y_k1.y,2))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));	
	
	return t_k_k_k_2;
	
}

float T1_k_k_k_12(threed y_k, threed y_k1)
/*<Second derivative of T with respect to x_k and y_k>*/
{
	float t_k_k_k_12,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k_k_12 = (sqrt(2)*pow(g0,3)*(y_k1.y-y_k.y)*(y_k.x-y_k1.x))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*pow(g0,3)*(y_k1.y-y_k.y)*(y_k.x-y_k1.x))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));		
	
	return t_k_k_k_12;
	
}

float T1_k_k_zk_1(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k and z_k>*/
{
	float t_k_k_zk_1,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k_zk_1 =  (sqrt(2)*pow(g0,3)*(y_k1.z-y_k.z)*(y_k.x-y_k1.x))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*pow(g0,3)*(y_k1.z-y_k.z)*(y_k.x-y_k1.x))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));		
	
	return t_k_k_zk_1;
	
}

float T1_k_k_zk_2(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to y_k and y_k>*/
{
	float t_k_k_zk_2,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k_zk_2 =  (sqrt(2)*pow(g0,3)*(y_k1.z-y_k.z)*(y_k.y-y_k1.y))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*pow(g0,3)*(y_k1.z-y_k.z)*(y_k.y-y_k1.y))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));		
	
	return t_k_k_zk_2;
	
}

float T1_k_zk_zk(threed y_k, threed y_k1)  
/*<Second Derivative of T with respect to z_k>*/
{
	float t_k_zk_zk,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_zk_zk = ((-1)*sqrt(2)*pow(g0,3)*pow(y_k.z-y_k1.z,2))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*g0)/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1))) - (sqrt(2)*pow(g0,3)*pow(y_k.z-y_k1.z,2))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));	
	
	return t_k_zk_zk;
	
}

/*k_k1_k1 Family & 1k_k_k Family (6)----------------------------------*/

float T1_k_k1_k1_1(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k1>*/
{
	float t_k_k1_k1_1,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k1_k1_1 = ((-1)*sqrt(2)*pow(g0,3)*pow(y_k.x-y_k1.x,2))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*g0)/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1))) - (sqrt(2)*pow(g0,3)*pow(y_k.x-y_k1.x,2))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));	
	
	return t_k_k1_k1_1;
	
}

float T1_k_k1_k1_2(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to y_k1>*/
{
	float t_k_k1_k1_2,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k1_k1_2 = ((-1)*sqrt(2)*pow(g0,3)*pow(y_k.y-y_k1.y,2))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*g0)/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1))) - (sqrt(2)*pow(g0,3)*pow(y_k.y-y_k1.y,2))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));	
	
	return t_k_k1_k1_2;
	
}

float T1_k_k1_k1_12(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k1 and y_k1>*/
{
	float t_k_k1_k1_12,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k1_k1_12 = (sqrt(2)*pow(g0,3)*(y_k1.y-y_k.y)*(y_k.x-y_k1.x))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*pow(g0,3)*(y_k1.y-y_k.y)*(y_k.x-y_k1.x))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));		
	
	return t_k_k1_k1_12;
	
}

float T1_k_k1_zk1_1(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to x_k1 and z_k1>*/
{
	float t_k_k1_zk1_1,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k1_zk1_1 = (sqrt(2)*pow(g0,3)*(y_k1.z-y_k.z)*(y_k.x-y_k1.x))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*pow(g0,3)*(y_k1.z-y_k.z)*(y_k.x-y_k1.x))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));		
	
	return t_k_k1_zk1_1;
	
}

float T1_k_k1_zk1_2(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to y_k1 and z_k1>*/
{
	float t_k_k1_zk1_2,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k1_zk1_2 = (sqrt(2)*pow(g0,3)*(y_k1.z-y_k.z)*(y_k.y-y_k1.y))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*pow(g0,3)*(y_k1.z-y_k.z)*(y_k.y-y_k1.y))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));		
	
	return t_k_k1_zk1_2;
	
}

float T1_k_zk1_zk1(threed y_k, threed y_k1)  
/*<Second Derivative of T with respect to z_k1>*/
{
	float t_k_zk1_zk1,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_zk1_zk1 = ((-1)*sqrt(2)*pow(g0,3)*pow(y_k.z-y_k1.z,2))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*g0)/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1))) - (sqrt(2)*pow(g0,3)*pow(y_k.z-y_k1.z,2))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));	
	
	return t_k_zk1_zk1;
	
}

/*k_k_k1 Family & 1k_1k_k Family (9)----------------------------------*/

float T1_k_k_k1_1(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to x_k and x_k1>*/
{
	float t_k_k_k1_1,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k_k1_1 = (sqrt(2)*pow(g0,3)*pow(y_k.x-y_k1.x,2))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) - (sqrt(2)*g0)/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1))) + (sqrt(2)*pow(g0,3)*pow(y_k.x-y_k1.x,2))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));	
	
	return t_k_k_k1_1;
	
}

float T1_k_k_k1_2(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to x_k and x_k1>*/
{
	float t_k_k_k1_2,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k_k1_2 = (sqrt(2)*pow(g0,3)*pow(y_k.y-y_k1.y,2))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) - (sqrt(2)*g0)/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1))) + (sqrt(2)*pow(g0,3)*pow(y_k.y-y_k1.y,2))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));	
	
	return t_k_k_k1_2;
	
}

float T1_k_k_k1_12(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to x_k and x_k1>*/
{
	float t_k_k_k1_12,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k_k1_12 = (sqrt(2)*pow(g0,3)*(y_k.x-y_k1.x)*(y_k.y-y_k1.y))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*pow(g0,3)*(y_k.x-y_k1.x)*(y_k.y-y_k1.y))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));		
	
	return t_k_k_k1_12;
	
}

float T1_k_k_k1_21(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to x_k and x_k1>*/
{
	float t_k_k_k1_21,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k_k1_21 = (sqrt(2)*pow(g0,3)*(y_k.x-y_k1.x)*(y_k.y-y_k1.y))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*pow(g0,3)*(y_k.x-y_k1.x)*(y_k.y-y_k1.y))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));		
	
	return t_k_k_k1_21;
	
}

float T1_k_k1_zk_1(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k1 and z_k>*/
{
	float t_k_k1_zk_1,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k1_zk_1 = (sqrt(2)*pow(g0,3)*(y_k.z-y_k1.z)*(y_k.x-y_k1.x))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*pow(g0,3)*(y_k.z-y_k1.z)*(y_k.x-y_k1.x))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));		
	
	return t_k_k1_zk_1;
	
}

float T1_k_k1_zk_2(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k1 and z_k>*/
{
	float t_k_k1_zk_2,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k1_zk_2 = (sqrt(2)*pow(g0,3)*(y_k.z-y_k1.z)*(y_k.y-y_k1.y))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*pow(g0,3)*(y_k.z-y_k1.z)*(y_k.y-y_k1.y))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));		
	
	return t_k_k1_zk_2;
	
}

float T1_k_k_zk1_1(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k and z_k1>*/
{
	float t_k_k_zk1_1,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k_zk1_1 = (sqrt(2)*pow(g0,3)*(y_k.z-y_k1.z)*(y_k.x-y_k1.x))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*pow(g0,3)*(y_k.z-y_k1.z)*(y_k.x-y_k1.x))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));		
	
	return t_k_k_zk1_1;
	
}

float T1_k_k_zk1_2(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k and z_k1>*/
{
	float t_k_k_zk1_2,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_k_zk1_2 = (sqrt(2)*pow(g0,3)*(y_k.z-y_k1.z)*(y_k.y-y_k1.y))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) + (sqrt(2)*pow(g0,3)*(y_k.z-y_k1.z)*(y_k.y-y_k1.y))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));		
	
	return t_k_k_zk1_2;
	
}

float T1_k_zk_zk1(threed y_k, threed y_k1)  
/*<Second Derivative of T with respect to z_k and z_k1>*/
{
	float t_k_zk_zk1,g0,q;
	
	q = pow(y_k.x-y_k1.x,2) + pow(y_k.y-y_k1.y,2) + pow(y_k.z-y_k1.z,2);
	
	g0 = pow(pow(y_k.gx2,2)+pow(y_k.gy2,2)+pow(y_k.gz2,2),0.5);
	
	t_k_zk_zk1 = (sqrt(2)*pow(g0,3)*pow(y_k.z-y_k1.z,2))/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*pow(y_k.v2,2)*pow(y_k1.v1,2)*pow(pow(g0,2)*q/(y_k.v2*y_k1.v1),3/2)) - (sqrt(2)*g0)/(sqrt(pow(g0,2)*q/(2*y_k.v2*y_k1.v1))+2*y_k.v2*y_k1.v1*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1))) + (sqrt(2)*pow(g0,3)*pow(y_k.z-y_k1.z,2))/(2*pow(pow(g0,2)*q/(2*y_k.v2*y_k1.v1)+2,3/2)*pow(y_k.v2,2)*pow(y_k1.v1,2)*sqrt(pow(g0,2)*q/(y_k.v2*y_k1.v1)));	
	
	return t_k_zk_zk1;
	
}



