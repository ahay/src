/*
 Copy_k1.yight (C) 2009 University of Texas at Austin
 
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

#include "ml_traveltime_vgradient_3D.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "general_traveltime_3D.h"
/*^*/

#ifndef _ml_traveltime_vgradient_3D_h

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
	float y_k.gx21;/* x-direction velocity gradient from above*/
	float y_k.gx22;/* x-direction velocity gradient from below*/
	float y_k.gy21;/* y-direction velocity gradient from above*/
	float y_k.gy22;/* y-direction velocity gradient from below*/
	float y_k.gz21;/* z-direction velocity gradient from above*/
	float y_k.gz22;/* z-direction velocity gradient from below*/
	float c111;/* c11 from above*/
	float c112;/* c11 from below*/
	float c331;/* c33 from above*/
	float c332;/* c33 from below*/
	float Q11; /* Q1 (horizontal anelliptic parameter) from above*/
	float Q12; /* Q1 (horizontal anelliptic parameter) from below*/
	float Q31; /* Q3 (vertical anelliptic parameter) from above*/
	float Q32; /* Q3 (vertical anelliptic parameter) from below*/
	float S11; /* S1 from above*/
	float S12; /* S1 from below*/
	float S31; /* S3 from above*/
	float S32; /* S3 from below*/
} threed;
/* Structure pointer */

#endif

/* Traveltime functions for gradient velocity------------------------------------------------------------------------------*/
double T1_k(threed y_k,threed y_k1)
/*<Traveltime>*/
{
	double t_k;
	
	
	t_k = log(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
        (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
     pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
            pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
              pow(-y_k1.z + y_k.z,2)))/2.,2),0.5))*
   pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5);
	
	return t_k;
	
}

/* First Derivative (6)---------------------------------------------------------------------------------------*/

double T1_k_k_1(threed y_k, threed y_k1) 
/*<Derivative of T with respect to x_k>*/
{
	double t_k_k_1;
	
	
	t_k_k_1 = ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (2*(-y_k1.x + y_k.x)*pow(y_k.v2,-1) - 
       y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
          pow(-y_k1.z + y_k.z,2)) - (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
          y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
             pow(-y_k1.z + y_k.z,2)))*
        (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
     pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
       pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1))/2.;
	
	return t_k_k_1;
	
}

double T1_k_k_2(threed y_k, threed y_k1) 
/*<Derivative of T with respect to y_k>*/
{
	double t_k_k_2;
	
	
	t_k_k_2 = ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (2*(-y_k1.y + y_k.y)*pow(y_k.v2,-1) - 
       y_k.gy2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
          pow(-y_k1.z + y_k.z,2)) - (2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
          y_k.gy2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
             pow(-y_k1.z + y_k.z,2)))*
        (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
     pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
       pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1))/2.;
	
	return t_k_k_2;
	
}

double T1_k_k1_1(threed y_k, threed y_k1) 
/*<Derivative of T with respect to x_k1>*/
{
	double t_k_k1_1;
	
	t_k_k1_1 = ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
       y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
          pow(-y_k1.z + y_k.z,2)) + (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
          y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
             pow(-y_k1.z + y_k.z,2)))*
        (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
     pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
       pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1))/2.;
	
	return t_k_k1_1;
	
}

double T1_k_k1_2(threed y_k, threed y_k1) 
/*<Derivative of T with respect to y_k1>*/
{
	double t_k_k1_2;
	

	
	t_k_k1_2 = ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
       y_k.gy2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
          pow(-y_k1.z + y_k.z,2)) + (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
          y_k.gy2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
             pow(-y_k1.z + y_k.z,2)))*
        (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
     pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
       pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1))/2.;
	
	return t_k_k1_2;
	
}

double T1_k_zk(threed y_k, threed y_k1)  
/*<Derivative of T with respect to z_k>*/
{
	double t_k_zk;
	
	
	t_k_zk = ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (2*(-y_k1.z + y_k.z)*pow(y_k.v2,-1) - 
       y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
          pow(-y_k1.z + y_k.z,2)) - (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
          y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
             pow(-y_k1.z + y_k.z,2)))*
        (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
     pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
       pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1))/2.;
	
	return t_k_zk;
	
}

double T1_k_zk1(threed y_k, threed y_k1) 
/*<Derivative of T with respect to z_k1>*/
{
	double t_k_zk1;
	
	
	t_k_zk1 =  ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
       y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
          pow(-y_k1.z + y_k.z,2)) + (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
          y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
             pow(-y_k1.z + y_k.z,2)))*
        (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
     pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
       pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1))/2.;
	
	return t_k_zk1;
	
}



/* Second Derivative------------------------------------------------------------------------------------------*/

/* k_k_k Family (6)-----------------------------------*/

double T1_k_k_k_1(threed y_k, threed y_k1)
/*<Second derivative of T with respect to x_k>*/
{
	double t_k_k_k_1;
	
	
	t_k_k_k_1 = pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
   (-(pow(((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
            (2*(-y_k1.x + y_k.x)*pow(y_k.v2,-1) - 
              y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)) - 
              (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
                 y_k.gx2*pow(y_k.v2,-2)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
               (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.)*pow(-1 + 
                 pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                      pow(y_k.v2,-1)*
                      (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                        pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5)))/2.,2)*
        pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
             pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-2)) + 
     (2*y_k.gx2*(y_k1.x - y_k.x)*(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
         pow(y_k.v2,-2) + (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
         pow(y_k.v2,-1) + pow(y_k.gx2,2)*(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
         pow(y_k1.v1,-1)*pow(y_k.v2,-3)*
         (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) - 
        pow(-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
                y_k.gx2*pow(y_k.v2,-2)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))/
           2.,2)*pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))/2.,2)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.,
            2),-1.5) + (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
         (2*y_k.gx2*(y_k1.x - y_k.x)*pow(y_k.v2,-2) + pow(y_k.v2,-1) + 
           pow(y_k.gx2,2)*pow(y_k.v2,-3)*
            (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
         (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))/2.)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.,
            2),-0.5) + pow(-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
              pow(y_k1.v1,-1)*(2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
                y_k.gx2*pow(y_k.v2,-2)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))/
           2.,2)*pow(-1 + pow(1 + 
             ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.,
            2),-0.5))*pow(1 + 
        ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1));	
	
	return t_k_k_k_1;
	
}

double T1_k_k_k_2(threed y_k, threed y_k1)
/*<Second derivative of T with respect to y_k>*/
{
	double t_k_k_k_2;
	

	t_k_k_k_2 = pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
   (-(pow(((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
            (2*(-y_k1.y + y_k.y)*pow(y_k.v2,-1) - 
              y_k.gy2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)) - 
              (2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
                 y_k.gy2*pow(y_k.v2,-2)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
               (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.)*pow(-1 + 
                 pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                      pow(y_k.v2,-1)*
                      (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                        pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5)))/2.,2)*
        pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
             pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-2)) + 
     (2*y_k.gy2*(y_k1.y - y_k.y)*(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
         pow(y_k.v2,-2) + (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
         pow(y_k.v2,-1) + pow(y_k.gy2,2)*(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
         pow(y_k1.v1,-1)*pow(y_k.v2,-3)*
         (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) - 
        pow(-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              (2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
                y_k.gy2*pow(y_k.v2,-2)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))/
           2.,2)*pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))/2.,2)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.,
            2),-1.5) + (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
         (2*y_k.gy2*(y_k1.y - y_k.y)*pow(y_k.v2,-2) + pow(y_k.v2,-1) + 
           pow(y_k.gy2,2)*pow(y_k.v2,-3)*
            (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
         (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))/2.)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.,
            2),-0.5) + pow(-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
              pow(y_k1.v1,-1)*(2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
                y_k.gy2*pow(y_k.v2,-2)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))/
           2.,2)*pow(-1 + pow(1 + 
             ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.,
            2),-0.5))*pow(1 + 
        ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1));	
	
	return t_k_k_k_2;
	
}

double T1_k_k_k_12(threed y_k, threed y_k1)
/*<Second derivative of T with respect to x_k and y_k>*/
{
	double t_k_k_k_12;
	
	
	t_k_k_k_12 = ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
          (2*(-y_k1.x + y_k.x)*pow(y_k.v2,-1) - 
            y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
               y_k.gx2*pow(y_k.v2,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          (2*(-y_k1.y + y_k.y)*pow(y_k.v2,-1) - 
            y_k.gy2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
               y_k.gy2*pow(y_k.v2,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                   pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                   (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                 2.,2),0.5),-2)) + 
       (4*y_k.gy2*(y_k1.x - y_k.x)*pow(y_k.v2,-2) + 4*y_k.gx2*(y_k1.y - y_k.y)*pow(y_k.v2,-2) + 
          4*y_k.gx2*y_k.gy2*pow(y_k.v2,-3)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) - 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
             y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
             y_k.gy2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-1.5) + 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
             y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
             y_k.gy2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5) + 
          4*((y_k.gy2*(y_k1.x - y_k.x) + y_k.gx2*(y_k1.y - y_k.y))*pow(y_k.v2,-2) + 
             y_k.gx2*y_k.gy2*pow(y_k.v2,-3)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
           (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5))*pow(1 + 
          ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-1)))/4.;		
	
	return t_k_k_k_12;
	
}

double T1_k_k_zk_1(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k and z_k>*/
{
	double t_k_k_zk_1;
	
	
	t_k_k_zk_1 =  ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
          (2*(-y_k1.x + y_k.x)*pow(y_k.v2,-1) - 
            y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
               y_k.gx2*pow(y_k.v2,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          (2*(-y_k1.z + y_k.z)*pow(y_k.v2,-1) - 
            y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
               y_k.gz2*pow(y_k.v2,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                   pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                   (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                 2.,2),0.5),-2)) + 
       (4*y_k.gz2*(y_k1.x - y_k.x)*pow(y_k.v2,-2) + 4*y_k.gx2*(y_k1.z - y_k.z)*pow(y_k.v2,-2) + 
          4*y_k.gx2*y_k.gz2*pow(y_k.v2,-3)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) - 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
             y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
             y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-1.5) + 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
             y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
             y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5) + 
          4*((y_k.gz2*(y_k1.x - y_k.x) + y_k.gx2*(y_k1.z - y_k.z))*pow(y_k.v2,-2) + 
             y_k.gx2*y_k.gz2*pow(y_k.v2,-3)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
           (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5))*pow(1 + 
          ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-1)))/4.;		
	
	return t_k_k_zk_1;
	
}

double T1_k_k_zk_2(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to y_k and y_k>*/
{
	double t_k_k_zk_2;
	
	
	t_k_k_zk_2 =  ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
          (2*(-y_k1.y + y_k.y)*pow(y_k.v2,-1) - 
            y_k.gy2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
               y_k.gy2*pow(y_k.v2,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          (2*(-y_k1.z + y_k.z)*pow(y_k.v2,-1) - 
            y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
               y_k.gz2*pow(y_k.v2,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                   pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                   (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                 2.,2),0.5),-2)) + 
       (4*y_k.gz2*(y_k1.y - y_k.y)*pow(y_k.v2,-2) + 4*y_k.gy2*(y_k1.z - y_k.z)*pow(y_k.v2,-2) + 
          4*y_k.gy2*y_k.gz2*pow(y_k.v2,-3)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) - 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
           (2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
             y_k.gy2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
             y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-1.5) + 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
           (2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
             y_k.gy2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
             y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5) + 
          4*((y_k.gz2*(y_k1.y - y_k.y) + y_k.gy2*(y_k1.z - y_k.z))*pow(y_k.v2,-2) + 
             y_k.gy2*y_k.gz2*pow(y_k.v2,-3)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
           (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5))*pow(1 + 
          ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-1)))/4.;		
	
	return t_k_k_zk_2;
	
}

double T1_k_zk_zk(threed y_k, threed y_k1)  
/*<Second Derivative of T with respect to z_k>*/
{
	double t_k_zk_zk;
	
	
	t_k_zk_zk = pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
   (-(pow(((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
            (2*(-y_k1.z + y_k.z)*pow(y_k.v2,-1) - 
              y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)) - 
              (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
                 y_k.gz2*pow(y_k.v2,-2)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
               (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.)*pow(-1 + 
                 pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                      pow(y_k.v2,-1)*
                      (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                        pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5)))/2.,2)*
        pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
             pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-2)) + 
     (2*y_k.gz2*(y_k1.z - y_k.z)*(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
         pow(y_k.v2,-2) + (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
         pow(y_k.v2,-1) + pow(y_k.gz2,2)*(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
         pow(y_k1.v1,-1)*pow(y_k.v2,-3)*
         (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) - 
        pow(-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
                y_k.gz2*pow(y_k.v2,-2)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))/
           2.,2)*pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))/2.,2)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.,
            2),-1.5) + (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
         (2*y_k.gz2*(y_k1.z - y_k.z)*pow(y_k.v2,-2) + pow(y_k.v2,-1) + 
           pow(y_k.gz2,2)*pow(y_k.v2,-3)*
            (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
         (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))/2.)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.,
            2),-0.5) + pow(-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
              pow(y_k1.v1,-1)*(2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
                y_k.gz2*pow(y_k.v2,-2)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))/
           2.,2)*pow(-1 + pow(1 + 
             ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.,
            2),-0.5))*pow(1 + 
        ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1));	
	
	return t_k_zk_zk;
	
}

/* k_k1_k1 Family & 1k_k_k Family (6)----------------------------------*/

double T1_k_k1_k1_1(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k1>*/
{
	double t_k_k1_k1_1;
	
	
	t_k_k1_k1_1 = pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
   (-(pow(((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
            (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
              y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)) + 
              (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
                 y_k.gx2*pow(y_k1.v1,-2)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
               (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.)*pow(-1 + 
                 pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                      pow(y_k.v2,-1)*
                      (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                        pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5)))/2.,2)*
        pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
             pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-2)) + 
     (2*y_k.gx2*(-y_k1.x + y_k.x)*(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-2)*
         pow(y_k.v2,-1) + (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
         pow(y_k.v2,-1) + pow(y_k.gx2,2)*(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
         pow(y_k1.v1,-3)*pow(y_k.v2,-1)*
         (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) - 
        pow(((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
             (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
               y_k.gx2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))/2.
           ,2)*pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))/2.,2)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.,
            2),-1.5) + (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
         (2*y_k.gx2*(-y_k1.x + y_k.x)*pow(y_k1.v1,-2) + pow(y_k1.v1,-1) + 
           pow(y_k.gx2,2)*pow(y_k1.v1,-3)*
            (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
         (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))/2.)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.,
            2),-0.5) + pow(((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
             pow(y_k.v2,-1)*(2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
               y_k.gx2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))/2.
           ,2)*pow(-1 + pow(1 + 
             ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.,
            2),-0.5))*pow(1 + 
        ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1));	
	
	return t_k_k1_k1_1;
	
}

double T1_k_k1_k1_2(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to y_k1>*/
{
	double t_k_k1_k1_2;
	
	
	t_k_k1_k1_2 = pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
   (-(pow(((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
            (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
              y_k.gy2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)) + 
              (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
                 y_k.gy2*pow(y_k1.v1,-2)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
               (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.)*pow(-1 + 
                 pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                      pow(y_k.v2,-1)*
                      (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                        pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5)))/2.,2)*
        pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
             pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-2)) + 
     (2*y_k.gy2*(-y_k1.y + y_k.y)*(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-2)*
         pow(y_k.v2,-1) + (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
         pow(y_k.v2,-1) + pow(y_k.gy2,2)*(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
         pow(y_k1.v1,-3)*pow(y_k.v2,-1)*
         (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) - 
        pow(((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
             (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
               y_k.gy2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))/2.
           ,2)*pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))/2.,2)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.,
            2),-1.5) + (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
         (2*y_k.gy2*(-y_k1.y + y_k.y)*pow(y_k1.v1,-2) + pow(y_k1.v1,-1) + 
           pow(y_k.gy2,2)*pow(y_k1.v1,-3)*
            (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
         (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))/2.)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.,
            2),-0.5) + pow(((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
             pow(y_k.v2,-1)*(2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
               y_k.gy2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))/2.
           ,2)*pow(-1 + pow(1 + 
             ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.,
            2),-0.5))*pow(1 + 
        ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1));	
	
	return t_k_k1_k1_2;
	
}

double T1_k_k1_k1_12(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k1 and y_k1>*/
{
	double t_k_k1_k1_12;
	
	
	t_k_k1_k1_12 = ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
          (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
            y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
               y_k.gx2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
            y_k.gy2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
               y_k.gy2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                   pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                   (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                 2.,2),0.5),-2)) + 
       (4*y_k.gy2*(-y_k1.x + y_k.x)*pow(y_k1.v1,-2) + 4*y_k.gx2*(-y_k1.y + y_k.y)*pow(y_k1.v1,-2) + 
          4*y_k.gx2*y_k.gy2*pow(y_k1.v1,-3)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) - 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
             y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
             y_k.gy2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-1.5) + 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
             y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
             y_k.gy2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5) + 
          4*((y_k.gy2*(-y_k1.x + y_k.x) + y_k.gx2*(-y_k1.y + y_k.y))*pow(y_k1.v1,-2) + 
             y_k.gx2*y_k.gy2*pow(y_k1.v1,-3)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
           (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5))*pow(1 + 
          ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-1)))/4.;		
	
	return t_k_k1_k1_12;
	
}

double T1_k_k1_zk1_1(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to x_k1 and z_k1>*/
{
	double t_k_k1_zk1_1;
	
	
	t_k_k1_zk1_1 = ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
          (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
            y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
               y_k.gx2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
            y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
               y_k.gz2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                   pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                   (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                 2.,2),0.5),-2)) + 
       (4*y_k.gz2*(-y_k1.x + y_k.x)*pow(y_k1.v1,-2) + 4*y_k.gx2*(-y_k1.z + y_k.z)*pow(y_k1.v1,-2) + 
          4*y_k.gx2*y_k.gz2*pow(y_k1.v1,-3)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) - 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
             y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
             y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-1.5) + 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
             y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
             y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5) + 
          4*((y_k.gz2*(-y_k1.x + y_k.x) + y_k.gx2*(-y_k1.z + y_k.z))*pow(y_k1.v1,-2) + 
             y_k.gx2*y_k.gz2*pow(y_k1.v1,-3)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
           (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5))*pow(1 + 
          ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-1)))/4.;		
	
	return t_k_k1_zk1_1;
	
}

double T1_k_k1_zk1_2(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to y_k1 and z_k1>*/
{
	double t_k_k1_zk1_2;
	
	
	t_k_k1_zk1_2 = ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
          (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
            y_k.gy2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
               y_k.gy2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
            y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
               y_k.gz2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                   pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                   (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                 2.,2),0.5),-2)) + 
       (4*y_k.gz2*(-y_k1.y + y_k.y)*pow(y_k1.v1,-2) + 4*y_k.gy2*(-y_k1.z + y_k.z)*pow(y_k1.v1,-2) + 
          4*y_k.gy2*y_k.gz2*pow(y_k1.v1,-3)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) - 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
           (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
             y_k.gy2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
             y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-1.5) + 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
           (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
             y_k.gy2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
             y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5) + 
          4*((y_k.gz2*(-y_k1.y + y_k.y) + y_k.gy2*(-y_k1.z + y_k.z))*pow(y_k1.v1,-2) + 
             y_k.gy2*y_k.gz2*pow(y_k1.v1,-3)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
           (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5))*pow(1 + 
          ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-1)))/4.;		
	
	return t_k_k1_zk1_2;
	
}

double T1_k_zk1_zk1(threed y_k, threed y_k1)  
/*<Second Derivative of T with respect to z_k1>*/
{
	double t_k_zk1_zk1;
	
	
	t_k_zk1_zk1 = pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
   (-(pow(((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
            (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
              y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)) + 
              (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
                 y_k.gz2*pow(y_k1.v1,-2)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
               (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.)*pow(-1 + 
                 pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                      pow(y_k.v2,-1)*
                      (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                        pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5)))/2.,2)*
        pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
             pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-2)) + 
     (2*y_k.gz2*(-y_k1.z + y_k.z)*(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-2)*
         pow(y_k.v2,-1) + (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
         pow(y_k.v2,-1) + pow(y_k.gz2,2)*(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
         pow(y_k1.v1,-3)*pow(y_k.v2,-1)*
         (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) - 
        pow(((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
             (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
               y_k.gz2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))/2.
           ,2)*pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))/2.,2)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.,
            2),-1.5) + (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
         (2*y_k.gz2*(-y_k1.z + y_k.z)*pow(y_k1.v1,-2) + pow(y_k1.v1,-1) + 
           pow(y_k.gz2,2)*pow(y_k1.v1,-3)*
            (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
         (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))/2.)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.,
            2),-0.5) + pow(((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
             pow(y_k.v2,-1)*(2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
               y_k.gz2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))/2.
           ,2)*pow(-1 + pow(1 + 
             ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.,
            2),-0.5))*pow(1 + 
        ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1));	
	
	return t_k_zk1_zk1;
	
}

/* k_k_k1 Family & 1k_1k_k Family (9)----------------------------------*/

double T1_k_k_k1_1(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to x_k and x_k1>*/
{
	double t_k_k_k1_1;
	
	
	t_k_k_k1_1 = ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
            y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
               y_k.gx2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          (2*(-y_k1.x + y_k.x)*pow(y_k.v2,-1) - 
            y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
               y_k.gx2*pow(y_k.v2,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                   pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                   (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                 2.,2),0.5),-2)) + 
       (4*y_k.gx2*(-y_k1.x + y_k.x)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) + 
          4*y_k.gx2*(y_k1.x - y_k.x)*pow(y_k1.v1,-2)*pow(y_k.v2,-1) - 
          4*pow(y_k1.v1,-1)*pow(y_k.v2,-1) + 
          2*pow(y_k.gx2,2)*pow(y_k1.v1,-2)*pow(y_k.v2,-2)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) + 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
             y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
             y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-1.5) - 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
             y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
             y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5) + 
          2*(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.)*
           (-2*pow(y_k1.v1,-1)*(y_k.gx2*(y_k1.x - y_k.x)*pow(y_k.v2,-2) + pow(y_k.v2,-1)) + 
             pow(y_k1.v1,-2)*(2*y_k.gx2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
                pow(y_k.gx2,2)*pow(y_k.v2,-2)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5))*pow(1 + 
          ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-1)))/4.;	
	
	return t_k_k_k1_1;
	
}

double T1_k_k_k1_2(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to y_k and y_k1>*/
{
	double t_k_k_k1_2;
	
	
	t_k_k_k1_2 = ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
            y_k.gy2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
               y_k.gy2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          (2*(-y_k1.y + y_k.y)*pow(y_k.v2,-1) - 
            y_k.gy2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
               y_k.gy2*pow(y_k.v2,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                   pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                   (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                 2.,2),0.5),-2)) + 
       (4*y_k.gy2*(-y_k1.y + y_k.y)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) + 
          4*y_k.gy2*(y_k1.y - y_k.y)*pow(y_k1.v1,-2)*pow(y_k.v2,-1) - 
          4*pow(y_k1.v1,-1)*pow(y_k.v2,-1) + 
          2*pow(y_k.gy2,2)*pow(y_k1.v1,-2)*pow(y_k.v2,-2)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) + 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
             y_k.gy2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
             y_k.gy2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-1.5) - 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
             y_k.gy2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
             y_k.gy2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5) + 
          2*(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.)*
           (-2*pow(y_k1.v1,-1)*(y_k.gy2*(y_k1.y - y_k.y)*pow(y_k.v2,-2) + pow(y_k.v2,-1)) + 
             pow(y_k1.v1,-2)*(2*y_k.gy2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
                pow(y_k.gy2,2)*pow(y_k.v2,-2)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5))*pow(1 + 
          ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-1)))/4.;	
	
	return t_k_k_k1_2;
	
}

double T1_k_k_k1_12(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to x_k and y_k1>*/
{
	double t_k_k_k1_12;
	
	
	t_k_k_k1_12 = ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
            y_k.gy2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
               y_k.gy2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          (2*(-y_k1.x + y_k.x)*pow(y_k.v2,-1) - 
            y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
               y_k.gx2*pow(y_k.v2,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                   pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                   (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                 2.,2),0.5),-2)) + 
       (4*y_k.gx2*(-y_k1.y + y_k.y)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) + 
          4*y_k.gy2*(y_k1.x - y_k.x)*pow(y_k1.v1,-2)*pow(y_k.v2,-1) + 
          2*y_k.gx2*y_k.gy2*pow(y_k1.v1,-2)*pow(y_k.v2,-2)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) + 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
             y_k.gy2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
             y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-1.5) - 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
             y_k.gy2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
             y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5) - 
          2*(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.)*
           (2*y_k.gx2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) - 
             y_k.gy2*pow(y_k1.v1,-2)*(2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
                y_k.gx2*pow(y_k.v2,-2)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5))*pow(1 + 
          ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-1)))/4.;		
	
	return t_k_k_k1_12;
	
}

double T1_k_k_k1_21(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to y_k and x_k1>*/
{
	double t_k_k_k1_21;
	
	
	t_k_k_k1_21 = ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
            y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
               y_k.gx2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          (2*(-y_k1.y + y_k.y)*pow(y_k.v2,-1) - 
            y_k.gy2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
               y_k.gy2*pow(y_k.v2,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                   pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                   (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                 2.,2),0.5),-2)) + 
       (4*y_k.gy2*(-y_k1.x + y_k.x)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) + 
          4*y_k.gx2*(y_k1.y - y_k.y)*pow(y_k1.v1,-2)*pow(y_k.v2,-1) + 
          2*y_k.gx2*y_k.gy2*pow(y_k1.v1,-2)*pow(y_k.v2,-2)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) + 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
             y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
             y_k.gy2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-1.5) - 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
             y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
             y_k.gy2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5) - 
          2*(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.)*
           (2*y_k.gy2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) - 
             y_k.gx2*pow(y_k1.v1,-2)*(2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
                y_k.gy2*pow(y_k.v2,-2)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5))*pow(1 + 
          ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-1)))/4.;		
	
	return t_k_k_k1_21;
	
}

double T1_k_k1_zk_1(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k1 and z_k>*/
{
	double t_k_k1_zk_1;
	
	
	t_k_k1_zk_1 = ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
            y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
               y_k.gx2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          (2*(-y_k1.z + y_k.z)*pow(y_k.v2,-1) - 
            y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
               y_k.gz2*pow(y_k.v2,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                   pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                   (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                 2.,2),0.5),-2)) + 
       (4*y_k.gz2*(-y_k1.x + y_k.x)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) + 
          4*y_k.gx2*(y_k1.z - y_k.z)*pow(y_k1.v1,-2)*pow(y_k.v2,-1) + 
          2*y_k.gx2*y_k.gz2*pow(y_k1.v1,-2)*pow(y_k.v2,-2)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) + 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
             y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
             y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-1.5) - 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
             y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
             y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5) - 
          2*(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.)*
           (2*y_k.gz2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) - 
             y_k.gx2*pow(y_k1.v1,-2)*(2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
                y_k.gz2*pow(y_k.v2,-2)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5))*pow(1 + 
          ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-1)))/4.;		
	
	return t_k_k1_zk_1;
	
}

double T1_k_k1_zk_2(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to y_k1 and z_k>*/
{
	double t_k_k1_zk_2;
	
	
	t_k_k1_zk_2 = ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
            y_k.gy2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
               y_k.gy2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          (2*(-y_k1.z + y_k.z)*pow(y_k.v2,-1) - 
            y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
               y_k.gz2*pow(y_k.v2,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                   pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                   (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                 2.,2),0.5),-2)) + 
       (4*y_k.gz2*(-y_k1.y + y_k.y)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) + 
          4*y_k.gy2*(y_k1.z - y_k.z)*pow(y_k1.v1,-2)*pow(y_k.v2,-1) + 
          2*y_k.gy2*y_k.gz2*pow(y_k1.v1,-2)*pow(y_k.v2,-2)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) + 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
             y_k.gy2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
             y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-1.5) - 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1) - 
             y_k.gy2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
             y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5) - 
          2*(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.)*
           (2*y_k.gz2*(y_k1.y - y_k.y)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) - 
             y_k.gy2*pow(y_k1.v1,-2)*(2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
                y_k.gz2*pow(y_k.v2,-2)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5))*pow(1 + 
          ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-1)))/4.;		
	
	return t_k_k1_zk_2;
	
}

double T1_k_k_zk1_1(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k and z_k1>*/
{
	double t_k_k_zk1_1;
	
	
	t_k_k_zk1_1 = ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
            y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
               y_k.gz2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          (2*(-y_k1.x + y_k.x)*pow(y_k.v2,-1) - 
            y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
               y_k.gx2*pow(y_k.v2,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                   pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                   (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                 2.,2),0.5),-2)) + 
       (4*y_k.gx2*(-y_k1.z + y_k.z)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) + 
          4*y_k.gz2*(y_k1.x - y_k.x)*pow(y_k1.v1,-2)*pow(y_k.v2,-1) + 
          2*y_k.gx2*y_k.gz2*pow(y_k1.v1,-2)*pow(y_k.v2,-2)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) + 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
             y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
             y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-1.5) - 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
             y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
             y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5) - 
          2*(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.)*
           (2*y_k.gx2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) - 
             y_k.gz2*pow(y_k1.v1,-2)*(2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
                y_k.gx2*pow(y_k.v2,-2)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5))*pow(1 + 
          ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-1)))/4.;		
	
	return t_k_k_zk1_1;
	
}

double T1_k_k_zk1_2(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to y_k and z_k1>*/
{
	double t_k_k_zk1_2;
	
	
	t_k_k_zk1_2 = ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
            y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
               y_k.gz2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          (2*(-y_k1.y + y_k.y)*pow(y_k.v2,-1) - 
            y_k.gy2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
               y_k.gy2*pow(y_k.v2,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                   pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                   (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                 2.,2),0.5),-2)) + 
       (4*y_k.gy2*(-y_k1.z + y_k.z)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) + 
          4*y_k.gz2*(y_k1.y - y_k.y)*pow(y_k1.v1,-2)*pow(y_k.v2,-1) + 
          2*y_k.gy2*y_k.gz2*pow(y_k1.v1,-2)*pow(y_k.v2,-2)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) + 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
             y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
             y_k.gy2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-1.5) - 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
             y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
             y_k.gy2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5) - 
          2*(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.)*
           (2*y_k.gy2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) - 
             y_k.gz2*pow(y_k1.v1,-2)*(2*(y_k1.y - y_k.y)*pow(y_k.v2,-1) + 
                y_k.gy2*pow(y_k.v2,-2)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5))*pow(1 + 
          ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-1)))/4.;		
	
	return t_k_k_zk1_2;
	
}

double T1_k_zk_zk1(threed y_k, threed y_k1)  
/*<Second Derivative of T with respect to z_k and z_k1>*/
{
	double t_k_zk_zk1;
	
	t_k_zk_zk1 = ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
     pow(pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
            y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
               y_k.gz2*pow(y_k1.v1,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          (2*(-y_k1.z + y_k.z)*pow(y_k.v2,-1) - 
            y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
               pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
               y_k.gz2*pow(y_k.v2,-2)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.)*pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                    pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))
                   /2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
               pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                 pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                   pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                   (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                 2.,2),0.5),-2)) + 
       (4*y_k.gz2*(-y_k1.z + y_k.z)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) + 
          4*y_k.gz2*(y_k1.z - y_k.z)*pow(y_k1.v1,-2)*pow(y_k.v2,-1) - 
          4*pow(y_k1.v1,-1)*pow(y_k.v2,-1) + 
          2*pow(y_k.gz2,2)*pow(y_k1.v1,-2)*pow(y_k.v2,-2)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)) + 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
             y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
             y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-1.5) - 
          (pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
             y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
             y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5) + 
          2*(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + 
                  pow(-y_k1.z + y_k.z,2)))/2.)*
           (-2*pow(y_k1.v1,-1)*(y_k.gz2*(y_k1.z - y_k.z)*pow(y_k.v2,-2) + pow(y_k.v2,-1)) + 
             pow(y_k1.v1,-2)*(2*y_k.gz2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
                pow(y_k.gz2,2)*pow(y_k.v2,-2)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2))))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                  pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/
                2.,2),-0.5))*pow(1 + 
          ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gy2,2) + pow(y_k.gz2,2))*
                 pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                 (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.y + y_k.y,2) + pow(-y_k1.z + y_k.z,2)))/2.
              ,2),0.5),-1)))/4.;	
	
	return t_k_zk_zk1;
	
}



