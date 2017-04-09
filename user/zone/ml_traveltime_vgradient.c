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

#include "general_traveltime.h"
/*^*/

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
	
	
	t_k = log(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
        (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
     pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
            (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),0.5))*
   pow(pow(y_k.gx2,2) + pow(y_k.gz2,2),-0.5);
	
	return t_k;
	
}

double T1_k_k(twod y_k, twod y_k1)
/*<Derivative of T with respect to x_k>*/
{
	double t_k_k;
	
	t_k_k = ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(pow(y_k.gx2,2) + pow(y_k.gz2,2),-0.5)*
     (2*(-y_k1.x + y_k.x)*pow(y_k.v2,-1) - 
       y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) - 
       (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
          y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
        (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
               (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
     pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
       pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1))/2.;
	
	return t_k_k;
	
}

double T1_k_k1(twod y_k, twod y_k1)
/*<Derivative of T with respect to x_k1>*/
{
	double t_k_k1;
	
	t_k_k1 = ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*pow(pow(y_k.gx2,2) + pow(y_k.gz2,2),-0.5)*
     (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
       y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) + 
       (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
          y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
        (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
               (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
     pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
       pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1))/2.;
	
	return t_k_k1;
	
}


double T1_k_zk(twod y_k, twod y_k1) 
/*<Derivative of T with respect to z_k>*/
{
	double t_k_zk;
	
	
	t_k_zk = ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(pow(y_k.gx2,2) + pow(y_k.gz2,2),-0.5)*
     (2*(-y_k1.z + y_k.z)*pow(y_k.v2,-1) - 
       y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) - 
       (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
          y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
        (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
               (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
     pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
       pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1))/2.;
	
	return t_k_zk;
	
}

double T1_k_zk1(twod y_k, twod y_k1) 
/*<Derivative of T with respect to z_k1>*/
{
	double t_k_zk1;
	
	t_k_zk1 = ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*pow(pow(y_k.gx2,2) + pow(y_k.gz2,2),-0.5)*
     (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
       y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) + 
       (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
          y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
        (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
               (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
     pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
       pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1))/2.;
	
	return t_k_zk1;
	
}


/*Second derivatives*/

double T1_k_k_k(twod y_k, twod y_k1)  
/*<Second derivative of T with respect to x_k>*/
{
	double t_k_k_k;
	
	t_k_k_k = pow(pow(y_k.gx2,2) + pow(y_k.gz2,2),-0.5)*
   (-(pow(((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
            (2*(-y_k1.x + y_k.x)*pow(y_k.v2,-1) - 
              y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) - 
              (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
                 y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
               (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
               pow(-1 + pow(1 + 
                   ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                      (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5)))/
          2.,2)*pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                 pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
           0.5),-2)) + (2*y_k.gx2*(y_k1.x - y_k.x)*(pow(y_k.gx2,2) + pow(y_k.gz2,2))*
         pow(y_k1.v1,-1)*pow(y_k.v2,-2) + 
        (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1) + 
        pow(y_k.gx2,2)*(pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-3)*
         (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) - 
        pow(-((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
                y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2))))/2.,2)
          *pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
          -1.5) + (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
         (2*y_k.gx2*(y_k1.x - y_k.x)*pow(y_k.v2,-2) + pow(y_k.v2,-1) + 
           pow(y_k.gx2,2)*pow(y_k.v2,-3)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
         (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
          -0.5) + pow(-((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
                y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2))))/2.,2)
          *pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
          -0.5))*pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
               (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1));	
	
	return t_k_k_k;
	
}

double T1_k_k1_k1(twod y_k, twod y_k1) 
/*<Second derivative of T with respect to x_k1>*/
{
	double t_k_k1_k1;
	
	t_k_k1_k1 = pow(pow(y_k.gx2,2) + pow(y_k.gz2,2),-0.5)*
   (-(pow(((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
            (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
              y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) + 
              (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
                 y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
               (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
               pow(-1 + pow(1 + 
                   ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                      (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5)))/
          2.,2)*pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                 pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
           0.5),-2)) + (2*y_k.gx2*(-y_k1.x + y_k.x)*(pow(y_k.gx2,2) + pow(y_k.gz2,2))*
         pow(y_k1.v1,-2)*pow(y_k.v2,-1) + 
        (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1) + 
        pow(y_k.gx2,2)*(pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-3)*pow(y_k.v2,-1)*
         (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) - 
        pow(((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
             (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
               y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2))))/2.,2)*
         pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
          -1.5) + (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
         (2*y_k.gx2*(-y_k1.x + y_k.x)*pow(y_k1.v1,-2) + pow(y_k1.v1,-1) + 
           pow(y_k.gx2,2)*pow(y_k1.v1,-3)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
         (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
          -0.5) + pow(((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
             (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
               y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2))))/2.,2)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
          -0.5))*pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
               (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1));		
	
	return t_k_k1_k1;
	
}

double T1_k_k_k1(twod y_k, twod y_k1) 
/*<Second derivative of T with respect to x_k and x_k1>*/
{
	double t_k_k_k1;
	
	t_k_k_k1 = ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(pow(y_k.gx2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
            y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
               y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
             pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
          (2*(-y_k1.x + y_k.x)*pow(y_k.v2,-1) - 
            y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
               y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
             pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
               (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                   pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
             0.5),-2)) + (4*y_k.gx2*(-y_k1.x + y_k.x)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) + 
          4*y_k.gx2*(y_k1.x - y_k.x)*pow(y_k1.v1,-2)*pow(y_k.v2,-1) - 
          4*pow(y_k1.v1,-1)*pow(y_k.v2,-1) + 
          2*pow(y_k.gx2,2)*pow(y_k1.v1,-2)*pow(y_k.v2,-2)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) + 
          (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
             y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
             y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
            -1.5) - (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
             y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
             y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
            -0.5) + (2 + pow(y_k.gx2,2)*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) + 
             pow(y_k.gz2,2)*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           (-2*pow(y_k1.v1,-1)*(y_k.gx2*(y_k1.x - y_k.x)*pow(y_k.v2,-2) + pow(y_k.v2,-1)) + 
             pow(y_k1.v1,-2)*(2*y_k.gx2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
                pow(y_k.gx2,2)*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))
             )*pow(-1 + pow(1 + 
               ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
        pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                 pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
           0.5),-1)))/4.;
	
	return t_k_k_k1;
	
}


double T1_k_zk_zk(twod y_k, twod y_k1)  
/*<Second Derivative of T with respect to z_k>*/
{
	double t_k_zk_zk;
	
	t_k_zk_zk =pow(pow(y_k.gx2,2) + pow(y_k.gz2,2),-0.5)*
   (-(pow(((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
            (2*(-y_k1.z + y_k.z)*pow(y_k.v2,-1) - 
              y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) - 
              (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
                 y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
               (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
               pow(-1 + pow(1 + 
                   ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                      (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5)))/
          2.,2)*pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                 pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
           0.5),-2)) + (2*y_k.gz2*(y_k1.z - y_k.z)*(pow(y_k.gx2,2) + pow(y_k.gz2,2))*
         pow(y_k1.v1,-1)*pow(y_k.v2,-2) + 
        (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1) + 
        pow(y_k.gz2,2)*(pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-3)*
         (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) - 
        pow(-((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
                y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2))))/2.,2)
          *pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
          -1.5) + (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
         (2*y_k.gz2*(y_k1.z - y_k.z)*pow(y_k.v2,-2) + pow(y_k.v2,-1) + 
           pow(y_k.gz2,2)*pow(y_k.v2,-3)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
         (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
          -0.5) + pow(-((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
              (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
                y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2))))/2.,2)
          *pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
          -0.5))*pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
               (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1));	
	
	return t_k_zk_zk;
	
}

double T1_k_zk1_zk1(twod y_k, twod y_k1)  
/*<Second Derivative of T with respect to z_k1>*/
{
	double t_k_zk1_zk1;
	
	t_k_zk1_zk1 = pow(pow(y_k.gx2,2) + pow(y_k.gz2,2),-0.5)*
   (-(pow(((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
            (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
              y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) + 
              (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
                 y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
               (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
               pow(-1 + pow(1 + 
                   ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                      (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5)))/
          2.,2)*pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                 pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
           0.5),-2)) + (2*y_k.gz2*(-y_k1.z + y_k.z)*(pow(y_k.gx2,2) + pow(y_k.gz2,2))*
         pow(y_k1.v1,-2)*pow(y_k.v2,-1) + 
        (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1) + 
        pow(y_k.gz2,2)*(pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-3)*pow(y_k.v2,-1)*
         (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) - 
        pow(((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
             (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
               y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2))))/2.,2)*
         pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
          -1.5) + (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
         (2*y_k.gz2*(-y_k1.z + y_k.z)*pow(y_k1.v1,-2) + pow(y_k1.v1,-1) + 
           pow(y_k.gz2,2)*pow(y_k1.v1,-3)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
         (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
          -0.5) + pow(((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
             (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
               y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2))))/2.,2)*
         pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
          -0.5))*pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
        pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
               (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),0.5),-1));	
	
	return t_k_zk1_zk1;
	
}

double T1_k_zk_zk1(twod y_k, twod y_k1)  
/*<Second Derivative of T with respect to z_k and z_k1>*/
{
	double t_k_zk_zk1;
	
	
	t_k_zk_zk1 = ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(pow(y_k.gx2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
            y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
               y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
             pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
          (2*(-y_k1.z + y_k.z)*pow(y_k.v2,-1) - 
            y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
               y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
             pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
               (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                   pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
             0.5),-2)) + (4*y_k.gz2*(-y_k1.z + y_k.z)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) + 
          4*y_k.gz2*(y_k1.z - y_k.z)*pow(y_k1.v1,-2)*pow(y_k.v2,-1) - 
          4*pow(y_k1.v1,-1)*pow(y_k.v2,-1) + 
          2*pow(y_k.gz2,2)*pow(y_k1.v1,-2)*pow(y_k.v2,-2)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) + 
          (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
             y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
             y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
            -1.5) - (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
             y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
             y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
            -0.5) + (2 + pow(y_k.gx2,2)*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) + 
             pow(y_k.gz2,2)*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           (-2*pow(y_k1.v1,-1)*(y_k.gz2*(y_k1.z - y_k.z)*pow(y_k.v2,-2) + pow(y_k.v2,-1)) + 
             pow(y_k1.v1,-2)*(2*y_k.gz2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
                pow(y_k.gz2,2)*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))
             )*pow(-1 + pow(1 + 
               ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
        pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
             (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                 pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
           0.5),-1)))/4.;	
	
	return t_k_zk_zk1;
	
}

double T1_k_k_zk(twod y_k, twod y_k1) 
/*<Second derivative of T with respect to x_k and z_k>*/
{
	double t_k_k_zk;
	
	t_k_k_zk = ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(pow(y_k.gx2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
          (2*(-y_k1.x + y_k.x)*pow(y_k.v2,-1) - 
            y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
               y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
             pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
          (2*(-y_k1.z + y_k.z)*pow(y_k.v2,-1) - 
            y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
               y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
             pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
               (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                   pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
             0.5),-2)) + (4*y_k.gz2*(y_k1.x - y_k.x)*pow(y_k.v2,-2) + 
          4*y_k.gx2*(y_k1.z - y_k.z)*pow(y_k.v2,-2) + 
          4*y_k.gx2*y_k.gz2*pow(y_k.v2,-3)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) - 
          (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
             y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
             y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
            -1.5) + (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
             y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
             y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
            -0.5) + 2*((y_k.gz2*(y_k1.x - y_k.x) + y_k.gx2*(y_k1.z - y_k.z))*pow(y_k.v2,-2) + 
             y_k.gx2*y_k.gz2*pow(y_k.v2,-3)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           (2 + pow(y_k.gx2,2)*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) + 
             pow(y_k.gz2,2)*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
              (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
            -0.5))*pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
             pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                 pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
           0.5),-1)))/4.;		
	
	return t_k_k_zk;
	
}

double T1_k_k1_zk1(twod y_k, twod y_k1) 
/*<Second derivative of T with respect to x_k1 and z_k1>*/
{
	double t_k_k1_zk1;
	
	t_k_k1_zk1 = ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*pow(pow(y_k.gx2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
          (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
            y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
               y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
             pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
          (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
            y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
               y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
             pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
               (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                   pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
             0.5),-2)) + (4*y_k.gz2*(-y_k1.x + y_k.x)*pow(y_k1.v1,-2) + 
          4*y_k.gx2*(-y_k1.z + y_k.z)*pow(y_k1.v1,-2) + 
          4*y_k.gx2*y_k.gz2*pow(y_k1.v1,-3)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) - 
          (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
             y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
             y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
            -1.5) + (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k.v2,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
             y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
             y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
            -0.5) + 4*((y_k.gz2*(-y_k1.x + y_k.x) + y_k.gx2*(-y_k1.z + y_k.z))*pow(y_k1.v1,-2) + 
             y_k.gx2*y_k.gz2*pow(y_k1.v1,-3)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
            -0.5))*pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
             pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                 pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
           0.5),-1)))/4.;		
	
	return t_k_k1_zk1;
	
}

double T1_k_k_zk1(twod y_k, twod y_k1)  
/*<Second derivative of T with respect to x_k and z_k1>*/
{
	double t_k_k_zk1;
	
	t_k_k_zk1 = ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(pow(y_k.gx2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
            y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
               y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
             pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
          (2*(-y_k1.x + y_k.x)*pow(y_k.v2,-1) - 
            y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
               y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
             pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
               (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                   pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
             0.5),-2)) + (4*y_k.gx2*(-y_k1.z + y_k.z)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) + 
          4*y_k.gz2*(y_k1.x - y_k.x)*pow(y_k1.v1,-2)*pow(y_k.v2,-1) + 
          2*y_k.gx2*y_k.gz2*pow(y_k1.v1,-2)*pow(y_k.v2,-2)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) + 
          (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
             y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
             y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
            -1.5) - (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1) - 
             y_k.gz2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
             y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
            -0.5) - 2*(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
           (2*y_k.gx2*(y_k1.z - y_k.z)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) - 
             y_k.gz2*pow(y_k1.v1,-2)*(2*(y_k1.x - y_k.x)*pow(y_k.v2,-1) + 
                y_k.gx2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2))))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
            -0.5))*pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
             pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                 pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
           0.5),-1)))/4.;	
	
	return t_k_k_zk1;
	
}

double T1_k_k1_zk(twod y_k, twod y_k1)  
/*<Second derivative of T with respect to x_k1 and z_k>*/
{
	double t_k_k1_zk;
	
	t_k_k1_zk = ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(pow(y_k.gx2,2) + pow(y_k.gz2,2),-0.5)*
     (-((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
          (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
            y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) + 
            (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
               y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
             pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
          (2*(-y_k1.z + y_k.z)*pow(y_k.v2,-1) - 
            y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) - 
            (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
               y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
             (1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                  (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
             pow(-1 + pow(1 + 
                 ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                    (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),-0.5))*
          pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
               (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
            pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                   pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
             0.5),-2)) + (4*y_k.gz2*(-y_k1.x + y_k.x)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) + 
          4*y_k.gx2*(y_k1.z - y_k.z)*pow(y_k1.v1,-2)*pow(y_k.v2,-1) + 
          2*y_k.gx2*y_k.gz2*pow(y_k1.v1,-2)*pow(y_k.v2,-2)*
           (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)) + 
          (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
             y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
             y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
                (pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2)*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
            -1.5) - (pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*pow(y_k.v2,-1)*
           (2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1) - 
             y_k.gx2*pow(y_k1.v1,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           (2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
             y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
            -0.5) - 2*(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.)*
           (2*y_k.gz2*(y_k1.x - y_k.x)*pow(y_k1.v1,-1)*pow(y_k.v2,-2) - 
             y_k.gx2*pow(y_k1.v1,-2)*(2*(y_k1.z - y_k.z)*pow(y_k.v2,-1) + 
                y_k.gz2*pow(y_k.v2,-2)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2))))*
           pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                  pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
            -0.5))*pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
             pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2. + 
          pow(-1 + pow(1 + ((pow(y_k.gx2,2) + pow(y_k.gz2,2))*pow(y_k1.v1,-1)*
                 pow(y_k.v2,-1)*(pow(-y_k1.x + y_k.x,2) + pow(-y_k1.z + y_k.z,2)))/2.,2),
           0.5),-1)))/4.;	
	
	return t_k_k1_zk;
	
}

