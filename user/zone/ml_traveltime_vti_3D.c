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

#include "ml_traveltime_vti_3D.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "general_traveltime_3D.h"
/*^*/

#ifndef _ml_traveltime_vti_3D_h

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

/*Traveltime functions for constant velocity---------------------------------------------------------------------------*/
double T2_k(threed y_k,threed y_k1)
/*<Traveltime>*/
{
	double t_k;
	
	t_k = sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
     (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))) + 
    ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
       sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
         (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
          (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
     (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)));
	
	return t_k;
	
}

/* First Derivative (6)---------------------------------------------------------------------------------------*/

double T2_k_k_1(threed y_k, threed y_k1) 
/*<Derivative of T with respect to x_k>*/
{
	double t_k_k_1;
	
	t_k_k_1 = (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
      ((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
        (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
         pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
     (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
           (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
     ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
        ((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
               (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
          (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
             (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
          (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
          (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
        sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
     (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
     (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
        sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))/
   (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
        (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
        + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));
	
	return t_k_k_1;
	
}

double T2_k_k_2(threed y_k, threed y_k1) 
/*<Derivative of T with respect to y_k>*/
{
	double t_k_k_2;
	
	t_k_k_2 = (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
      ((2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
        (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
         pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
     (2*(-y_k.y + y_k1.y)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
           (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
     ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
        ((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
               (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
          (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
             (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
          (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
          (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
        sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
     (2*y_k.S12*(-y_k.y + y_k1.y)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
     (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
        sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))/
   (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
        (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
        + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));
	
	return t_k_k_2;
	
}

double T2_k_k1_1(threed y_k, threed y_k1) 
/*<Derivative of T with respect to x_k1>*/
{
	double t_k_k1_1;
	
	t_k_k1_1 = (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
      ((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
        (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
         pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
     (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
           (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
     ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
        ((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
               (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
          (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
             (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
          (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
          (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
        sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
     (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
     (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
        sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))/
   (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
        (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
        + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));
	
	return t_k_k1_1;
	
}

double T2_k_k1_2(threed y_k, threed y_k1) 
/*<Derivative of T with respect to y_k1>*/
{
	double t_k_k1_2;
	
	t_k_k1_2 = (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
      ((-2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
        (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
         pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
     (2*(-y_k.y + y_k1.y)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
           (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
     ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
        ((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
               (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
          (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
             (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
          (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
          (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
        sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
     (2*y_k.S12*(-y_k.y + y_k1.y)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
     (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
        sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))/
   (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
        (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
        + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));
	
	return t_k_k1_2;
	
}

double T2_k_zk(threed y_k, threed y_k1)  
/*<Derivative of T with respect to z_k>*/
{
	double t_k_zk;
	
	
	t_k_zk = (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
      ((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
        (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
         pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
     (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
           (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
     ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
        ((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
               (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
          (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
          (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
          (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
        sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
     (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
     (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
        sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))/
   (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
        (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
        + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));
	
	return t_k_zk;
	
}

double T2_k_zk1(threed y_k, threed y_k1) 
/*<Derivative of T with respect to z_k1>*/
{
	double t_k_zk1;
	
	
	t_k_zk1 = (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
      ((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
        (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
         pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
     (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
           (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
     ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
        ((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
               (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
          (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
          (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
          (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
        sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
     (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
     (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
        sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))/
   (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
        (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
        + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));
	
	return t_k_zk1;
	
}



/* Second Derivative------------------------------------------------------------------------------------------*/

/* k_k_k Family (6)-----------------------------------*/

double T2_k_k_k_1(threed y_k, threed y_k1)
/*<Second derivative of T with respect to x_k>*/
{
	double t_k_k_k_1;
	
	t_k_k_k_1 = -pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
        ((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
          (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
           pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
       (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
             (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
       ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          ((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                 (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                  pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
            (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
            (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
            (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
          sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
       (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
       (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2),2)/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((8*y_k.S12*pow(-y_k.x + y_k1.x,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (2*y_k.S12)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (8*pow(-y_k.x + y_k1.x,2)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) + 
         (2*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
      (4*(-y_k.x + y_k1.x)*((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 + 
      (2*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
            (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         pow((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),2))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((8*pow(-y_k.x + y_k1.x,2))/pow(y_k.c112,2) + (4*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-8*y_k.Q12*pow(-y_k.x + y_k1.x,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (2*y_k.Q12)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (8*pow(-y_k.x + y_k1.x,2)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) - 
                (2*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (16*pow(y_k.S12,2)*pow(-y_k.x + y_k1.x,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) - 
           (16*y_k.S12*pow(-y_k.x + y_k1.x,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (16*y_k.S12*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (16*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (2*y_k.S12*(-y_k.x + y_k1.x)*((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (8*y_k.S12*pow(-y_k.x + y_k1.x,2)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (2*y_k.S12*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
      (8*pow(-y_k.x + y_k1.x,2)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) - 
      (2*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_k_k_1;
	
}

double T2_k_k_k_2(threed y_k, threed y_k1)
/*<Second derivative of T with respect to y_k>*/
{
	double t_k_k_k_2;
	
	t_k_k_k_2 = -pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
        ((2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
          (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
           pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
       (2*(-y_k.y + y_k1.y)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
             (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
       ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          ((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                 (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                  pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
            (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
            (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
            (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
          sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
       (2*y_k.S12*(-y_k.y + y_k1.y)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
       (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2),2)/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((8*y_k.S12*pow(-y_k.y + y_k1.y,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (2*y_k.S12)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (8*pow(-y_k.y + y_k1.y,2)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) + 
         (2*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
      (4*(-y_k.y + y_k1.y)*((2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 + 
      (2*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
            (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         pow((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),2))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((8*pow(-y_k.y + y_k1.y,2))/pow(y_k.c112,2) + (4*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-8*y_k.Q12*pow(-y_k.y + y_k1.y,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (2*y_k.Q12)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (8*pow(-y_k.y + y_k1.y,2)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) - 
                (2*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (16*pow(y_k.S12,2)*pow(-y_k.y + y_k1.y,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) - 
           (16*y_k.S12*pow(-y_k.y + y_k1.y,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (16*y_k.S12*pow(-y_k.y + y_k1.y,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (16*pow(-y_k.y + y_k1.y,2)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (2*y_k.S12*(-y_k.y + y_k1.y)*((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (8*y_k.S12*pow(-y_k.y + y_k1.y,2)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (2*y_k.S12*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
      (8*pow(-y_k.y + y_k1.y,2)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) - 
      (2*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_k_k_2;
	
}

double T2_k_k_k_12(threed y_k, threed y_k1)
/*<Second derivative of T with respect to x_k and y_k>*/
{
	double t_k_k_k_12;
	
	t_k_k_k_12 = -((((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))*
       (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.y + y_k1.y)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S12*(-y_k.y + y_k1.y)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((8*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (8*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)) - 
      (2*(-y_k.y + y_k1.y)*((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      (2*(-y_k.x + y_k1.x)*((2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((8*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y))/pow(y_k.c112,2) + (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-8*y_k.Q12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (8*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (16*pow(y_k.S12,2)*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) - 
           (16*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (16*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (16*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))
                ))/(y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S12*(-y_k.y + y_k1.y)*((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S12*(-y_k.x + y_k1.x)*((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (8*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (8*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));
	
	return t_k_k_k_12;
	
}

double T2_k_k_zk_1(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k and z_k>*/
{
	double t_k_k_zk_1;
	
	t_k_k_zk_1 =  -((((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))*
       (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((4*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
         (4*y_k.S32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)) - 
      (2*(-y_k.z + y_k1.z)*((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 - 
      (2*(-y_k.x + y_k1.x)*((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/(y_k.c112*y_k.c332) + (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-4*y_k.Q12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
                (4*y_k.Q32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (16*y_k.S12*y_k.S32*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) - 
           (8*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*y_k.S32*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*y_k.S32*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S12*(-y_k.x + y_k1.x)*((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S32*(-y_k.z + y_k1.z)*((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (4*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
      (4*y_k.S32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));	
	
	return t_k_k_zk_1;
	
}

double T2_k_k_zk_2(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to y_k and y_k>*/
{
	double t_k_k_zk_2;
	
	t_k_k_zk_2 =  -((((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.y + y_k1.y)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S12*(-y_k.y + y_k1.y)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))*
       (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((4*y_k.S12*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
         (4*y_k.S32*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)) - 
      (2*(-y_k.z + y_k1.z)*((2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 - 
      (2*(-y_k.y + y_k1.y)*((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/(y_k.c112*y_k.c332) + (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-4*y_k.Q12*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
                (4*y_k.Q32*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (16*y_k.S12*y_k.S32*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) - 
           (8*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*y_k.S32*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*y_k.S32*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S12*(-y_k.y + y_k1.y)*((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S32*(-y_k.z + y_k1.z)*((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (4*y_k.S12*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
      (4*y_k.S32*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));	
	
	return t_k_k_zk_2;
	
}

double T2_k_zk_zk(threed y_k, threed y_k1)  
/*<Second Derivative of T with respect to z_k>*/
{
	double t_k_zk_zk;
	
	t_k_zk_zk = -pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
        ((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
          (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
           pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
       (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
             (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
       ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          ((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                 (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                  pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
            (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
            (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
            (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
          sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
       (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
       (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2),2)/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((8*y_k.S32*pow(-y_k.z + y_k1.z,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (2*y_k.S32)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (8*pow(-y_k.z + y_k1.z,2)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) + 
         (2*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
      (4*(-y_k.z + y_k1.z)*((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 + 
      (2*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
            (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         pow((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),2))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((8*pow(-y_k.z + y_k1.z,2))/pow(y_k.c332,2) + (4*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-8*y_k.Q32*pow(-y_k.z + y_k1.z,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (2*y_k.Q32)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (8*pow(-y_k.z + y_k1.z,2)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) - 
                (2*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (16*pow(y_k.S32,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,4)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) - 
           (16*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,4)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (20*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (20*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (2*y_k.S32*(-y_k.z + y_k1.z)*((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (8*y_k.S32*pow(-y_k.z + y_k1.z,2)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (2*y_k.S32*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
      (8*pow(-y_k.z + y_k1.z,2)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) - 
      (2*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_zk_zk;
	
}

/* k_k1_k1 Family & 1k_k_k Family (6)----------------------------------*/

double T2_k_k1_k1_1(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k1>*/
{
	double t_k_k1_k1_1;
	
	t_k_k1_k1_1 = -pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
        ((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
          (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
           pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
       (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
             (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
       ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          ((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                 (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                  pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
            (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
            (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
            (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
          sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
       (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
       (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2),2)/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((8*y_k.S12*pow(-y_k.x + y_k1.x,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (2*y_k.S12)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (8*pow(-y_k.x + y_k1.x,2)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) + 
         (2*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
      (4*(-y_k.x + y_k1.x)*((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 + 
      (2*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
            (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         pow((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),2))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((8*pow(-y_k.x + y_k1.x,2))/pow(y_k.c112,2) + (4*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-8*y_k.Q12*pow(-y_k.x + y_k1.x,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (2*y_k.Q12)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (8*pow(-y_k.x + y_k1.x,2)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) - 
                (2*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (16*pow(y_k.S12,2)*pow(-y_k.x + y_k1.x,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) - 
           (16*y_k.S12*pow(-y_k.x + y_k1.x,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (16*y_k.S12*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (16*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (2*y_k.S12*(-y_k.x + y_k1.x)*((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (8*y_k.S12*pow(-y_k.x + y_k1.x,2)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (2*y_k.S12*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
      (8*pow(-y_k.x + y_k1.x,2)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) - 
      (2*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_k1_k1_1;
	
}

double T2_k_k1_k1_2(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to y_k1>*/
{
	double t_k_k1_k1_2;
	
	t_k_k1_k1_2 = -pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
        ((-2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
          (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
           pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
       (2*(-y_k.y + y_k1.y)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
             (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
       ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          ((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                 (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                  pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
            (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
            (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
            (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
          sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
       (2*y_k.S12*(-y_k.y + y_k1.y)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
       (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2),2)/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((8*y_k.S12*pow(-y_k.y + y_k1.y,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (2*y_k.S12)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (8*pow(-y_k.y + y_k1.y,2)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) + 
         (2*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
      (4*(-y_k.y + y_k1.y)*((-2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 + 
      (2*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
            (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         pow((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),2))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((8*pow(-y_k.y + y_k1.y,2))/pow(y_k.c112,2) + (4*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-8*y_k.Q12*pow(-y_k.y + y_k1.y,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (2*y_k.Q12)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (8*pow(-y_k.y + y_k1.y,2)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) - 
                (2*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (16*pow(y_k.S12,2)*pow(-y_k.y + y_k1.y,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) - 
           (16*y_k.S12*pow(-y_k.y + y_k1.y,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (16*y_k.S12*pow(-y_k.y + y_k1.y,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (16*pow(-y_k.y + y_k1.y,2)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (2*y_k.S12*(-y_k.y + y_k1.y)*((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (8*y_k.S12*pow(-y_k.y + y_k1.y,2)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (2*y_k.S12*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
      (8*pow(-y_k.y + y_k1.y,2)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) - 
      (2*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_k1_k1_2;
	
}

double T2_k_k1_k1_12(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k1 and y_k1>*/
{
	double t_k_k1_k1_12;
	
	t_k_k1_k1_12 = -((((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))*
       (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((-2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.y + y_k1.y)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S12*(-y_k.y + y_k1.y)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((8*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (8*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)) + 
      (2*(-y_k.y + y_k1.y)*((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 + 
      (2*(-y_k.x + y_k1.x)*((-2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((8*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y))/pow(y_k.c112,2) + (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-8*y_k.Q12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (8*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (16*pow(y_k.S12,2)*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) - 
           (16*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (16*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (16*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))
                ))/(y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S12*(-y_k.y + y_k1.y)*((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S12*(-y_k.x + y_k1.x)*((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (8*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (8*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));
	
	return t_k_k1_k1_12;
	
}

double T2_k_k1_zk1_1(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to x_k1 and z_k1>*/
{
	double t_k_k1_zk1_1;
	
	t_k_k1_zk1_1 = -((((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))*
       (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((4*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
         (4*y_k.S32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)) + 
      (2*(-y_k.z + y_k1.z)*((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 + 
      (2*(-y_k.x + y_k1.x)*((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/(y_k.c112*y_k.c332) + (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-4*y_k.Q12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
                (4*y_k.Q32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (16*y_k.S12*y_k.S32*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) - 
           (8*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*y_k.S32*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*y_k.S32*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S12*(-y_k.x + y_k1.x)*((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S32*(-y_k.z + y_k1.z)*((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (4*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
      (4*y_k.S32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_k1_zk1_1;
	
}

double T2_k_k1_zk1_2(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to y_k1 and z_k1>*/
{
	double t_k_k1_zk1_2;
	
	t_k_k1_zk1_2 = -((((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((-2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.y + y_k1.y)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S12*(-y_k.y + y_k1.y)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))*
       (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((4*y_k.S12*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
         (4*y_k.S32*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)) + 
      (2*(-y_k.z + y_k1.z)*((-2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 + 
      (2*(-y_k.y + y_k1.y)*((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/(y_k.c112*y_k.c332) + (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-4*y_k.Q12*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
                (4*y_k.Q32*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (16*y_k.S12*y_k.S32*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) - 
           (8*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*y_k.S32*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*y_k.S32*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S12*(-y_k.y + y_k1.y)*((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S32*(-y_k.z + y_k1.z)*((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (4*y_k.S12*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
      (4*y_k.S32*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));	
	
	return t_k_k1_zk1_2;
	
}

double T2_k_zk1_zk1(threed y_k, threed y_k1)  
/*<Second Derivative of T with respect to z_k1>*/
{
	double t_k_zk1_zk1;
	
	t_k_zk1_zk1 = -pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
        ((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
          (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
           pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
       (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
             (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
       ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          ((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                 (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                  pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
            (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
            (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
            (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
          sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
       (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
       (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                  (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2),2)/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((8*y_k.S32*pow(-y_k.z + y_k1.z,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (2*y_k.S32)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (8*pow(-y_k.z + y_k1.z,2)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) + 
         (2*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
      (4*(-y_k.z + y_k1.z)*((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 + 
      (2*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
            (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         pow((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),2))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((8*pow(-y_k.z + y_k1.z,2))/pow(y_k.c332,2) + (4*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-8*y_k.Q32*pow(-y_k.z + y_k1.z,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (2*y_k.Q32)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (8*pow(-y_k.z + y_k1.z,2)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) - 
                (2*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (16*pow(y_k.S32,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,4)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) - 
           (16*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,4)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (20*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (20*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (2*y_k.S32*(-y_k.z + y_k1.z)*((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (8*y_k.S32*pow(-y_k.z + y_k1.z,2)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (2*y_k.S32*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
      (8*pow(-y_k.z + y_k1.z,2)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) - 
      (2*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_zk1_zk1;
	
}

/* k_k_k1 Family & 1k_1k_k Family (9)----------------------------------*/

double T2_k_k_k1_1(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to x_k and x_k1>*/
{
	double t_k_k_k1_1;
	
	t_k_k_k1_1 = -((((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))*
       (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((-8*y_k.S12*pow(-y_k.x + y_k1.x,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
         (2*y_k.S12)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
         (8*pow(-y_k.x + y_k1.x,2)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) - 
         (2*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
      (2*(-y_k.x + y_k1.x)*((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      (2*(-y_k.x + y_k1.x)*((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      (2*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
            (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-8*pow(-y_k.x + y_k1.x,2))/pow(y_k.c112,2) - (4*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((8*y_k.Q12*pow(-y_k.x + y_k1.x,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
                (2*y_k.Q12)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (8*pow(-y_k.x + y_k1.x,2)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) + 
                (2*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (16*pow(y_k.S12,2)*pow(-y_k.x + y_k1.x,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) + 
           (16*y_k.S12*pow(-y_k.x + y_k1.x,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (16*y_k.S12*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (16*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S12*(-y_k.x + y_k1.x)*((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S12*(-y_k.x + y_k1.x)*((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (8*y_k.S12*pow(-y_k.x + y_k1.x,2)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
      (2*y_k.S12*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
      (8*pow(-y_k.x + y_k1.x,2)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) + 
      (2*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_k_k1_1;
	
}

double T2_k_k_k1_2(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to x_k and x_k1>*/
{
	double t_k_k_k1_2;
	
	t_k_k_k1_2 = -((((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((-2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.y + y_k1.y)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S12*(-y_k.y + y_k1.y)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))*
       (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.y + y_k1.y)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S12*(-y_k.y + y_k1.y)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((-8*y_k.S12*pow(-y_k.y + y_k1.y,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
         (2*y_k.S12)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
         (8*pow(-y_k.y + y_k1.y,2)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) - 
         (2*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
      (2*(-y_k.y + y_k1.y)*((2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      (2*(-y_k.y + y_k1.y)*((-2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      (2*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
            (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-8*pow(-y_k.y + y_k1.y,2))/pow(y_k.c112,2) - (4*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((8*y_k.Q12*pow(-y_k.y + y_k1.y,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
                (2*y_k.Q12)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (8*pow(-y_k.y + y_k1.y,2)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) + 
                (2*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (16*pow(y_k.S12,2)*pow(-y_k.y + y_k1.y,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) + 
           (16*y_k.S12*pow(-y_k.y + y_k1.y,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (16*y_k.S12*pow(-y_k.y + y_k1.y,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (16*pow(-y_k.y + y_k1.y,2)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S12*(-y_k.y + y_k1.y)*((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S12*(-y_k.y + y_k1.y)*((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (8*y_k.S12*pow(-y_k.y + y_k1.y,2)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
      (2*y_k.S12*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
      (8*pow(-y_k.y + y_k1.y,2)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) + 
      (2*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_k_k1_2;
	
}

double T2_k_k_k1_12(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to x_k and x_k1>*/
{
	double t_k_k_k1_12;
	
	t_k_k_k1_12 = -((((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))*
       (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((-2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.y + y_k1.y)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S12*(-y_k.y + y_k1.y)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((-8*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
         (8*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)) + 
      (2*(-y_k.y + y_k1.y)*((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      (2*(-y_k.x + y_k1.x)*((-2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-8*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y))/pow(y_k.c112,2) + (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((8*y_k.Q12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
                (8*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (16*pow(y_k.S12,2)*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) + 
           (16*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (16*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (16*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))
                ))/(y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S12*(-y_k.y + y_k1.y)*((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S12*(-y_k.x + y_k1.x)*((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (8*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
      (8*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));	
	
	return t_k_k_k1_12;
	
}

double T2_k_k_k1_21(threed y_k, threed y_k1) 
/*<Second derivative of T with respect to x_k and x_k1>*/
{
	double t_k_k_k1_21;
	
	t_k_k_k1_21 = -((((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))*
       (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.y + y_k1.y)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S12*(-y_k.y + y_k1.y)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((-8*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
         (8*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)) - 
      (2*(-y_k.y + y_k1.y)*((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 + 
      (2*(-y_k.x + y_k1.x)*((2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-8*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y))/pow(y_k.c112,2) + (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((8*y_k.Q12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
                (8*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (16*pow(y_k.S12,2)*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) + 
           (16*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (16*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (16*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))
                ))/(y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S12*(-y_k.y + y_k1.y)*((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S12*(-y_k.x + y_k1.x)*((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (8*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
      (8*(-y_k.x + y_k1.x)*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));	
	
	return t_k_k_k1_21;
	
}

double T2_k_k1_zk_1(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k1 and z_k>*/
{
	double t_k_k1_zk_1;
	
	t_k_k1_zk_1 = -((((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))*
       (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((-4*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (4*y_k.S32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
         (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)) - 
      (2*(-y_k.z + y_k1.z)*((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 + 
      (2*(-y_k.x + y_k1.x)*((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/(y_k.c112*y_k.c332) + (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((4*y_k.Q12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (4*y_k.Q32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
                (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (16*y_k.S12*y_k.S32*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) + 
           (8*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*y_k.S32*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*y_k.S32*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S12*(-y_k.x + y_k1.x)*((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S32*(-y_k.z + y_k1.z)*((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (4*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (4*y_k.S32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
      (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));	
	
	return t_k_k1_zk_1;
	
}

double T2_k_k1_zk_2(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k1 and z_k>*/
{
	double t_k_k1_zk_2;
	
	t_k_k1_zk_2 = -((((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((-2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.y + y_k1.y)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S12*(-y_k.y + y_k1.y)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))*
       (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((-4*y_k.S12*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (4*y_k.S32*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
         (8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)) - 
      (2*(-y_k.z + y_k1.z)*((-2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 + 
      (2*(-y_k.y + y_k1.y)*((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/(y_k.c112*y_k.c332) + (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((4*y_k.Q12*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (4*y_k.Q32*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
                (8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (16*y_k.S12*y_k.S32*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) + 
           (8*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*y_k.S32*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*y_k.S32*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S12*(-y_k.y + y_k1.y)*((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S32*(-y_k.z + y_k1.z)*((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (4*y_k.S12*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (4*y_k.S32*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
      (8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));	
	
	return t_k_k1_zk_2;
	
}

double T2_k_k_zk1_1(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k and z_k1>*/
{
	double t_k_k_zk1_1;
	
	t_k_k_zk1_1 = -((((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))*
       (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((-4*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (4*y_k.S32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
         (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)) + 
      (2*(-y_k.z + y_k1.z)*((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 - 
      (2*(-y_k.x + y_k1.x)*((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/(y_k.c112*y_k.c332) + (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((4*y_k.Q12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (4*y_k.Q32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
                (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (16*y_k.S12*y_k.S32*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) + 
           (8*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*y_k.S32*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*y_k.S32*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S12*(-y_k.x + y_k1.x)*((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.x + y_k1.x)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S32*(-y_k.z + y_k1.z)*((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.x + y_k1.x)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (4*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (4*y_k.S32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
      (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));	
	
	return t_k_k_zk1_1;
	
}

double T2_k_k_zk1_2(threed y_k, threed y_k1)  
/*<Second derivative of T with respect to x_k and z_k1>*/
{
	double t_k_k_zk1_2;
	
	t_k_k_zk1_2 = -((((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.y + y_k1.y)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S12*(-y_k.y + y_k1.y)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))*
       (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((-4*y_k.S12*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (4*y_k.S32*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
         (8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)) + 
      (2*(-y_k.z + y_k1.z)*((2*y_k.S12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 - 
      (2*(-y_k.y + y_k1.y)*((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/(y_k.c112*y_k.c332) + (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((4*y_k.Q12*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (4*y_k.Q32*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
                (8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (16*y_k.S12*y_k.S32*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) + 
           (8*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*y_k.S32*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*y_k.S32*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S12*(-y_k.y + y_k1.y)*((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.y + y_k1.y)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S32*(-y_k.z + y_k1.z)*((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.y + y_k1.y)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.y + y_k1.y))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.y + y_k1.y)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(-y_k.y + y_k1.y)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.y + y_k1.y)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (4*y_k.S12*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (4*y_k.S32*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
      (8*(-y_k.y + y_k1.y)*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));	
	
	return t_k_k_zk1_2;
	
}

double T2_k_zk_zk1(threed y_k, threed y_k1)  
/*<Second Derivative of T with respect to z_k and z_k1>*/
{
	double t_k_zk_zk1;
	
	t_k_zk_zk1 = -((((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))*
       (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
          ((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
             pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
               (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
         ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            ((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                    pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                    (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   (((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
       ((-8*y_k.S32*pow(-y_k.z + y_k1.z,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) + 
         (2*y_k.S32)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
         (8*pow(-y_k.z + y_k1.z,2)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) - 
         (2*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
          pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)) + 
      (2*(-y_k.z + y_k1.z)*((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 - 
      (2*(-y_k.z + y_k1.z)*((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/
            pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 - 
      (2*(1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/
            (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 - 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         pow(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-8*pow(-y_k.z + y_k1.z,2))/pow(y_k.c332,2) - (4*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((8*y_k.Q32*pow(-y_k.z + y_k1.z,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
                (2*y_k.Q32)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (8*pow(-y_k.z + y_k1.z,2)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) + 
                (2*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (16*pow(y_k.S32,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,4)*
              (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) + 
           (16*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,4)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (20*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (20*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S32*(-y_k.z + y_k1.z)*((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S32*(-y_k.z + y_k1.z)*((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.z + y_k1.z)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.z + y_k1.z)*((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/
                 pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,3)*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2)*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (8*y_k.S32*pow(-y_k.z + y_k1.z,2)*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2) - 
      (2*y_k.S32*sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)) - 
      (8*pow(-y_k.z + y_k1.z,2)*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),3) + 
      (2*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2),2))/
    (2.*sqrt(((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2)))\
         + ((y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
           sqrt(pow((pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2))*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                   (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2)) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
         (pow(-y_k.x + y_k1.x,2) + pow(-y_k.y + y_k1.y,2) + pow(-y_k.z + y_k1.z,2))));
	
	return t_k_zk_zk1;
	
}

