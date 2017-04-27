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

#include "general_traveltime.h"
/*^*/


/* Traveltime functions for velocity in TI media------------------------------------------------------------------------------*/

double T2_k(twod y_k,twod y_k1)
/*<Traveltime>*/
{
	double t_k;
	
	t_k = sqrt((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*(1 - 
       (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
    ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
         (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
            (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
          (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)));
	return t_k;
	
}

double T2_k_k(twod y_k, twod y_k1)
/*<Derivative of T with respect to x_k>*/
{
	double t_k_k;
	
	t_k_k = ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
        (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) - 
     (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
     ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((-4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
          (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
               (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
          (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
          (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
          (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
     (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
     (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
        sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2))/
   (2.*sqrt((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
        (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
       ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))));
	
	
	return t_k_k;
	
}

double T2_k_k1(twod y_k, twod y_k1)
/*<Derivative of T with respect to x_k1>*/
{
	double t_k_k1;
	
	
	t_k_k1 = ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
        (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) + 
     (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
     ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
          (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
               (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
          (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
          (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
          (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
     (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
     (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
        sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2))/
   (2.*sqrt((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
        (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
       ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_k1;
	
}

double T2_k_k_k(twod y_k, twod y_k1)  
/*<Second derivative of T with respect to x_k>*/
{
	double t_k_k_k;
	
	t_k_k_k = -pow((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
          (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) - 
       (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
       ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((-4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                 (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
            (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
            (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                 (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
            (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
       (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
       (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2),2)/
    (4.*pow((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((8*y_k.S12*pow(-y_k.x + y_k1.x,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (2*y_k.S12)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
         (8*pow(-y_k.x + y_k1.x,2)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3) + 
         (2*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) - 
      (4*(-y_k.x + y_k1.x)*((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 + 
      (2*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 - 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*pow((-4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))),2))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*pow(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((8*pow(-y_k.x + y_k1.x,2))/pow(y_k.c112,2) + 
           (4*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-8*y_k.Q12*pow(-y_k.x + y_k1.x,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + (2*y_k.Q12)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (8*pow(-y_k.x + y_k1.x,2)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3) - 
                (2*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (16*pow(y_k.S12,2)*pow(-y_k.x + y_k1.x,4)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) - 
           (16*y_k.S12*pow(-y_k.x + y_k1.x,4)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (20*y_k.S12*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (20*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (2*y_k.S12*(-y_k.x + y_k1.x)*((-4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (8*y_k.S12*pow(-y_k.x + y_k1.x,2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (2*y_k.S12*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
      (8*pow(-y_k.x + y_k1.x,2)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3) - 
      (2*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2))/
    (2.*sqrt((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_k_k;
	
}

double T2_k_k1_k1(twod y_k, twod y_k1) 
/*<Second derivative of T with respect to x_k1>*/
{
	double t_k_k1_k1;
	
	t_k_k1_k1 = -pow((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
          (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) + 
       (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
       ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                 (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
            (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
            (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                 (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
            (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
       (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
       (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2),2)/
    (4.*pow((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((8*y_k.S12*pow(-y_k.x + y_k1.x,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (2*y_k.S12)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
         (8*pow(-y_k.x + y_k1.x,2)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3) + 
         (2*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) + 
      (4*(-y_k.x + y_k1.x)*((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 + 
      (2*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 - 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*pow((4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))),2))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*pow(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((8*pow(-y_k.x + y_k1.x,2))/pow(y_k.c112,2) + 
           (4*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-8*y_k.Q12*pow(-y_k.x + y_k1.x,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + (2*y_k.Q12)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (8*pow(-y_k.x + y_k1.x,2)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3) - 
                (2*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (16*pow(y_k.S12,2)*pow(-y_k.x + y_k1.x,4)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) - 
           (16*y_k.S12*pow(-y_k.x + y_k1.x,4)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (20*y_k.S12*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (20*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (2*y_k.S12*(-y_k.x + y_k1.x)*((4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (8*y_k.S12*pow(-y_k.x + y_k1.x,2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (2*y_k.S12*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
      (8*pow(-y_k.x + y_k1.x,2)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3) - 
      (2*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2))/
    (2.*sqrt((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_k1_k1;
	
}

double T2_k_k_k1(twod y_k, twod y_k1) 
/*<Second derivative of T with respect to x_k and x_k1>*/
{
    double t_k_k_k1;
	
/*	g0 = hypotf(y_k.gx2,y_k.gz2); */
	
	t_k_k_k1 = -(((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                   (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2))*
       ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((-4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                   (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((-8*y_k.S12*pow(-y_k.x + y_k1.x,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + 
         (2*y_k.S12)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
         (8*pow(-y_k.x + y_k1.x,2)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3) - 
         (2*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) + 
      (2*(-y_k.x + y_k1.x)*((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      (2*(-y_k.x + y_k1.x)*((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      (2*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 - 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((-4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*pow(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((-8*pow(-y_k.x + y_k1.x,2))/pow(y_k.c112,2) - 
           (4*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((8*y_k.Q12*pow(-y_k.x + y_k1.x,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - (2*y_k.Q12)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (8*pow(-y_k.x + y_k1.x,2)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3) + 
                (2*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (16*pow(y_k.S12,2)*pow(-y_k.x + y_k1.x,4)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) + 
           (16*y_k.S12*pow(-y_k.x + y_k1.x,4)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (20*y_k.S12*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (20*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S12*(-y_k.x + y_k1.x)*((-4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S12*(-y_k.x + y_k1.x)*((4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (8*y_k.S12*pow(-y_k.x + y_k1.x,2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - 
      (2*y_k.S12*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
      (8*pow(-y_k.x + y_k1.x,2)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3) + 
      (2*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2))/
    (2.*sqrt((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_k_k1;
	
}

double T2_k_zk(twod y_k, twod y_k1) 
/*<Derivative of T with respect to z_k>*/
{
	double t_k_zk;
	
	t_k_zk = ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
        (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) - 
     (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
     ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((-4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
          (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
               (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
          (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
          (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
          (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
     (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
     (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
        sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2))/
   (2.*sqrt((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
        (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
       ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_zk;
	
}

double T2_k_zk1(twod y_k, twod y_k1) 
/*<Derivative of T with respect to z_k1>*/
{
	double t_k_zk1;
	
	t_k_zk1 = ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
        (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) + 
     (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
     ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
          (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
               (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
          (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
          (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
          (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
      (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
     (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
     (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
        sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
          (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
             (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
           (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2))/
   (2.*sqrt((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
        (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
       ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_zk1;
	
}

double T2_k_zk_zk(twod y_k, twod y_k1)  
/*<Second Derivative of T with respect to z_k>*/
{
	double t_k_zk_zk;
	
	t_k_zk_zk = -pow((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
          (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) - 
       (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
       ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((-4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                 (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
            (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
            (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + 
                 (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
            (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
       (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
       (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2),2)/
    (4.*pow((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((8*y_k.S32*pow(-y_k.z + y_k1.z,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (2*y_k.S32)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
         (8*pow(-y_k.z + y_k1.z,2)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3) + 
         (2*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) - 
      (4*(-y_k.z + y_k1.z)*((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 + 
      (2*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 - 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*pow((-4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))),2))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*pow(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((8*pow(-y_k.z + y_k1.z,2))/pow(y_k.c332,2) + 
           (4*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-8*y_k.Q32*pow(-y_k.z + y_k1.z,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + (2*y_k.Q32)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (8*pow(-y_k.z + y_k1.z,2)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3) - 
                (2*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (16*pow(y_k.S32,2)*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,4)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) - 
           (16*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,4)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (20*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (20*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.x + y_k1.x,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (2*y_k.S32*(-y_k.z + y_k1.z)*((-4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (8*y_k.S32*pow(-y_k.z + y_k1.z,2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (2*y_k.S32*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
      (8*pow(-y_k.z + y_k1.z,2)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3) - 
      (2*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2))/
    (2.*sqrt((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_zk_zk;
	
}

double T2_k_zk1_zk1(twod y_k, twod y_k1)  
/*<Second Derivative of T with respect to z_k1>*/
{
	double t_k_zk1_zk1;
	
	t_k_zk1_zk1 = -pow((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
          (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) + 
       (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
       ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                 (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
            (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
            (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + 
                 (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
            (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
        (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
       (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
       (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
          sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
            (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
               (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
             (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2),2)/
    (4.*pow((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((8*y_k.S32*pow(-y_k.z + y_k1.z,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (2*y_k.S32)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
         (8*pow(-y_k.z + y_k1.z,2)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3) + 
         (2*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) + 
      (4*(-y_k.z + y_k1.z)*((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 + 
      (2*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 - 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*pow((4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))),2))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*pow(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((8*pow(-y_k.z + y_k1.z,2))/pow(y_k.c332,2) + 
           (4*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-8*y_k.Q32*pow(-y_k.z + y_k1.z,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + (2*y_k.Q32)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (8*pow(-y_k.z + y_k1.z,2)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3) - 
                (2*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (16*pow(y_k.S32,2)*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,4)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) - 
           (16*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,4)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (20*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (20*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.x + y_k1.x,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (2*y_k.S32*(-y_k.z + y_k1.z)*((4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (8*y_k.S32*pow(-y_k.z + y_k1.z,2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (2*y_k.S32*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
      (8*pow(-y_k.z + y_k1.z,2)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3) - 
      (2*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2))/
    (2.*sqrt((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_zk1_zk1;
	
}

double T2_k_zk_zk1(twod y_k, twod y_k1)  
/*<Second Derivative of T with respect to z_k and z_k1>*/
{
	double t_k_zk_zk1;
	
	
	t_k_zk_zk1 = -(((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
         ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + 
                   (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2))*
       ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
         ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((-4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + 
                   (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((-8*y_k.S32*pow(-y_k.z + y_k1.z,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + 
         (2*y_k.S32)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
         (8*pow(-y_k.z + y_k1.z,2)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3) - 
         (2*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) + 
      (2*(-y_k.z + y_k1.z)*((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 - 
      (2*(-y_k.z + y_k1.z)*((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 - 
      (2*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 - 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((-4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*pow(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((-8*pow(-y_k.z + y_k1.z,2))/pow(y_k.c332,2) - 
           (4*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((8*y_k.Q32*pow(-y_k.z + y_k1.z,2))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - (2*y_k.Q32)/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (8*pow(-y_k.z + y_k1.z,2)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3) + 
                (2*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (16*pow(y_k.S32,2)*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,4)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) + 
           (16*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,4)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (20*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (20*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*pow(-y_k.x + y_k1.x,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S32*(-y_k.z + y_k1.z)*((-4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S32*(-y_k.z + y_k1.z)*((4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (8*y_k.S32*pow(-y_k.z + y_k1.z,2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - 
      (2*y_k.S32*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
      (8*pow(-y_k.z + y_k1.z,2)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3) + 
      (2*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2))/
    (2.*sqrt((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_zk_zk1;
	
}

double T2_k_k_zk(twod y_k, twod y_k1) 
/*<Second derivative of T with respect to x_k and z_k>*/
{
	double t_k_k_zk;
	
	t_k_k_zk = -(((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((-4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                   (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2))*
       ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
         ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((-4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + 
                   (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((4*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + 
         (4*y_k.S32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3)) - 
      (2*(-y_k.z + y_k1.z)*((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 - 
      (2*(-y_k.x + y_k1.x)*((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((-4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((-4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*pow(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/(y_k.c112*y_k.c332) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-4*y_k.Q12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - 
                (4*y_k.Q32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (16*y_k.S12*y_k.S32*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) - 
           (8*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,3)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*y_k.S32*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,3)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*y_k.S12*pow(-y_k.x + y_k1.x,3)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*y_k.S32*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*pow(-y_k.x + y_k1.x,3)*(-y_k.z + y_k1.z)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S12*(-y_k.x + y_k1.x)*((-4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S32*(-y_k.z + y_k1.z)*((-4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (4*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - 
      (4*y_k.S32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3))/
    (2.*sqrt((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_k_zk;
	
}

double T2_k_k1_zk1(twod y_k, twod y_k1) 
/*<Second derivative of T with respect to x_k1 and z_k1>*/
{
	double t_k_k1_zk1;
	
	t_k_k1_zk1 = -(((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                   (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2))*
       ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
         ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + 
                   (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((4*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + 
         (4*y_k.S32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3)) + 
      (2*(-y_k.z + y_k1.z)*((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 + 
      (2*(-y_k.x + y_k1.x)*((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*pow(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/(y_k.c112*y_k.c332) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-4*y_k.Q12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - 
                (4*y_k.Q32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (16*y_k.S12*y_k.S32*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) - 
           (8*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,3)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*y_k.S32*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,3)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*y_k.S12*pow(-y_k.x + y_k1.x,3)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*y_k.S32*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*pow(-y_k.x + y_k1.x,3)*(-y_k.z + y_k1.z)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S12*(-y_k.x + y_k1.x)*((4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S32*(-y_k.z + y_k1.z)*((4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (4*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - 
      (4*y_k.S32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3))/
    (2.*sqrt((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_k1_zk1;
	
}

double T2_k_k_zk1(twod y_k, twod y_k1)  
/*<Second derivative of T with respect to x_k and z_k1>*/
{
	double t_k_k_zk1;
	
	
	t_k_k_zk1 = -(((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((-4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                   (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2))*
       ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
         ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + 
                   (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((-4*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (4*y_k.S32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + 
         (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3)) + 
      (2*(-y_k.z + y_k1.z)*((2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 - 
      (2*(-y_k.x + y_k1.x)*((-2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((-4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*pow(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((-8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/(y_k.c112*y_k.c332) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((4*y_k.Q12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (4*y_k.Q32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - 
                (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (16*y_k.S12*y_k.S32*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) + 
           (8*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,3)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*y_k.S32*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,3)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*y_k.S12*pow(-y_k.x + y_k1.x,3)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*y_k.S32*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*pow(-y_k.x + y_k1.x,3)*(-y_k.z + y_k1.z)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S12*(-y_k.x + y_k1.x)*((4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S32*(-y_k.z + y_k1.z)*((-4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (4*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (4*y_k.S32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - 
      (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3))/
    (2.*sqrt((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_k_zk1;
	
}

double T2_k_k1_zk(twod y_k, twod y_k1)  
/*<Second derivative of T with respect to x_k1 and z_k>*/
{
	double t_k_k1_zk;
	
	t_k_k1_zk = -(((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
            (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) + 
         (2*(-y_k.x + y_k1.x)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c112 + 
         ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                   (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
              (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + 
                   (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
         (2*y_k.S12*(-y_k.x + y_k1.x)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
         (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2))*
       ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
            (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)) - 
         (2*(-y_k.z + y_k1.z)*(1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/y_k.c332 + 
         ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((-4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                   (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
              (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
              (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + 
                   (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
              (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
          (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
         (2*y_k.S32*(-y_k.z + y_k1.z)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
         (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
            sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
              (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                 (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
               (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
    (4.*pow((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)),1.5)) + 
   ((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*((-4*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - 
         (4*y_k.S32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + 
         (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3)) - 
      (2*(-y_k.z + y_k1.z)*((-2*y_k.S12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
           (2*(-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c332 + 
      (2*(-y_k.x + y_k1.x)*((2*y_k.S32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
           (2*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/y_k.c112 - 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((-4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))*
         ((4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (4.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*pow(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))),1.5)) + 
      ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*((-8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/(y_k.c112*y_k.c332) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((4*y_k.Q12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + 
                (4*y_k.Q32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - 
                (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (16*y_k.S12*y_k.S32*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),3)) + 
           (8*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,3)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*y_k.S32*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,3)*(-1 + 
                (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*y_k.S12*pow(-y_k.x + y_k1.x,3)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (8*y_k.S32*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (8*pow(-y_k.x + y_k1.x,3)*(-y_k.z + y_k1.z)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (2.*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (y_k.S12*(-y_k.x + y_k1.x)*((-4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      ((-y_k.x + y_k1.x)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((-4*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c332 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((-2*y_k.Q32*(-y_k.z + y_k1.z))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) + 
                (2*(-y_k.z + y_k1.z)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*y_k.S32*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) - 
           (4*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,3)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*pow(-y_k.x + y_k1.x,2)*(-y_k.z + y_k1.z)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) - 
      (y_k.S32*(-y_k.z + y_k1.z)*((4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       ((pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      ((-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         ((4*(-y_k.x + y_k1.x)*(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332))/y_k.c112 + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              ((2*y_k.Q12*(-y_k.x + y_k1.x))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2)) - 
                (2*(-y_k.x + y_k1.x)*(y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2)))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) - 
           (4*y_k.S12*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*pow(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2),2)) + 
           (4*pow(-y_k.x + y_k1.x,3)*pow(-y_k.z + y_k1.z,2)*(-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/
                 (pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/(y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))) + 
           (4*(-y_k.x + y_k1.x)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/
       (pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))))) + 
      (4*y_k.S12*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) + 
      (4*y_k.S32*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),2) - 
      (8*(-y_k.x + y_k1.x)*(-y_k.z + y_k1.z)*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*
         sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
           (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
              (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
            (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/pow(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2),3))/
    (2.*sqrt((pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332)*
         (1 - (y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))) + 
        ((y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2))*sqrt(pow(pow(-y_k.x + y_k1.x,2)/y_k.c112 + pow(-y_k.z + y_k1.z,2)/y_k.c332,2) + 
             (2*pow(-y_k.x + y_k1.x,2)*pow(-y_k.z + y_k1.z,2)*(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))*
                (-1 + (y_k.Q12*pow(-y_k.x + y_k1.x,2) + y_k.Q32*pow(-y_k.z + y_k1.z,2))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))))/
              (y_k.c112*y_k.c332*(y_k.S12*pow(-y_k.x + y_k1.x,2) + y_k.S32*pow(-y_k.z + y_k1.z,2)))))/(pow(-y_k.x + y_k1.x,2) + pow(-y_k.z + y_k1.z,2))));
	return t_k_k1_zk;
	
}

