/* Generalized Dix inversion and quartic parameter estimation*/
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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrixops_2D.h"


void gendix(float **inverted, float **eff, float *twtime, int nlayer) {
/*<Generalized Dix inversion componentwise>*/

int i;
float **A, **Bold, **B;

A = sf_floatalloc2(2,2); Bold = sf_floatalloc2(2,2); B = sf_floatalloc2(2,2);

for (i = 0; i< nlayer; i++) {

    A[0][0] = eff[i][0]; A[0][1] = eff[i][1]; A[1][0] = eff[i][1]; A[1][1] = eff[i][2];
    mat_inverse(A);
    
    if (A[1][0] != A[0][1]) sf_warning("Gendix: B is not symmetric %f and %f",A[1][0], A[0][1]);

    if (i != 0 ){
        B[0][0] = (twtime[i]*A[0][0] - twtime[i-1]*Bold[0][0])/((twtime[i]-twtime[i-1]));
        B[0][1] = (twtime[i]*A[1][0] - twtime[i-1]*Bold[1][0])/((twtime[i]-twtime[i-1])); B[1][0] = B[0][1];
        B[1][1] = (twtime[i]*A[1][1] - twtime[i-1]*Bold[1][1])/((twtime[i]-twtime[i-1]));
        
        mat_inverse(B);
        
        inverted[i][0] = B[0][0];
        inverted[i][1] = B[0][1];
        inverted[i][2] = B[1][1];
        
        Bold[0][0] = A[0][0]; Bold[0][1] = A[0][1]; Bold[1][0] = A[1][0]; Bold[1][1] = A[1][1];
    } else { /* First layer */
        Bold[0][0] = A[0][0]; Bold[0][1] = A[0][1]; Bold[1][0] = A[1][0]; Bold[1][1] = A[1][1];
        mat_inverse(A);
        inverted[i][0] = A[0][0];
        inverted[i][1] = A[1][0];
        inverted[i][2] = A[1][1];
    }
}

}


void quartic(float **inverted, float **eff, float *twtime, int nlayer) {
/*<Quartic coefficient inversion componentwise>*/

int i;
float a11,a12,a22;
float a1111,a1112,a1121,a1122,a1211,a1212,a1221,a1222,a2111,a2112,a2121,a2122,a2211,a2212,a2221,a2222;
float h1111,h1112,h1121,h1122,h1211,h1212,h1221,h1222,h2111,h2112,h2121,h2122,h2211,h2212,h2221,h2222;
float h1111old,h1112old,h1121old,h1122old,h1211old,h1212old,h1221old,h1222old,h2111old,h2112old,h2121old,h2122old,h2211old,h2212old,h2221old,h2222old;
float h1111int,h1112int,h1121int,h1122int,h1211int,h1212int,h1221int,h1222int,h2111int,h2112int,h2121int,h2122int,h2211int,h2212int,h2221int,h2222int;


/* Loop over layers */
for (i = 0; i< nlayer; i++) {

    /* Due to the symmetry of time tensor */
    a11 = eff[i][0]; a12 = eff[i][1]; a22 = eff[i][2]; 
    a1111 = eff[i][3]; 
    a1112 = eff[i][4]; a1121 = a1112; a1211 = a1112; a2111 = a1112; 
    a1122 = eff[i][5]; a1221 = a1122; a2211 = a1122; a1212 = a1122; a2121 = a1122; a2112 = a1122;
    a1222 = eff[i][6]; a2122 = a1222; a2212 = a1222; a2221 = a1222; 
    a2222 = eff[i][7];
    
        /*fourth-ranked h tensor componentwise (only five independent parameters)*/
         h1111 = (3*pow(pow(a12,2) - a11*a22,2)*(pow(a12,2) + pow(a22,2))*twtime[i] + 12*(pow(a12,3)*(-(a1122*a12) + a11*a1222) + pow(a12,2)*(-2*a11*a1122 + a12*(2*a1112 + a1222))*a22 + a12*(a11*a1112 - (a1111 + 3*a1122)*a12)*pow(a22,2) + 3*a1112*a12*pow(a22,3) - a1111*pow(a22,4))*pow(twtime[i],3))/(2.*pow(pow(a12,2) - a11*a22,4));
                
        h1112 = (-6*(pow(a12,4)*a1222 + a1112*pow(a22,4) - a1122*a12*pow(a22,2)*(a11 + 3*a22) + pow(a12,2)*a22*(2*a11*a1222 + (a1112 + 3*a1222)*a22) - pow(a12,3)*(2*a1122*a22 + (a11 + a22)*a2222))*pow(twtime[i],3))/pow(pow(a12,2) - a11*a22,4);
                
        h1121 = (-3*a12*(a11 + a22)*pow(pow(a12,2) - a11*a22,2)*twtime[i] + 12*(pow(a12,2)*(a11*a1122*a12 - (pow(a11,2) + pow(a12,2))*a1222) + a12*(2*pow(a11,2)*a1122 - 2*a11*a1112*a12 + 3*a1122*pow(a12,2))*a22 + (-(pow(a11,2)*a1112) + a11*a1111*a12 - 3*a1112*pow(a12,2))*pow(a22,2) + a1111*a12*pow(a22,3))*pow(twtime[i],3))/(2.*pow(pow(a12,2) - a11*a22,4));
                
        h1122 = (-6*(-(a12*a22*(3*pow(a12,2)*a1222 - 3*a1122*a12*a22 + a1112*pow(a22,2))) - a11*a12*(pow(a12,2)*a1222 - 2*a1122*a12*a22 + a1112*pow(a22,2)) + pow(a12,4)*a2222 + pow(a11,2)*(-2*a12*a1222*a22 + a1122*pow(a22,2) + pow(a12,2)*a2222))*pow(twtime[i],3))/pow(pow(a12,2) - a11*a22,4);
               
        h1211 = (6*(pow(a11,2)*a12*(-(a12*a1222) + a1122*a22) + a11*(-(pow(a12,2)*(2*a1112 + a1222)*a22) - a1112*pow(a22,3) + 2*a1122*a12*(pow(a12,2) + pow(a22,2))) + a12*((a1111 + a1122)*pow(a12,2)*a22 + a1111*pow(a22,3) - a1112*(pow(a12,3) + 2*a12*pow(a22,2))))*pow(twtime[i],3))/pow(pow(a12,2) - a11*a22,4);
               
        h1212 = (3*pow(pow(a12,2) - a11*a22,2)*(pow(a12,2) + pow(a22,2))*twtime[i] - 12*(a1122*(pow(a12,4) + a11*pow(a22,3) + 2*pow(a12,2)*a22*(a11 + a22)) - a12*(pow(a11,2)*a1222*a22 + pow(a12,2)*(a1112 + a1222)*a22 + a1112*pow(a22,3) + 2*a11*a1222*(pow(a12,2) + pow(a22,2))) + a11*pow(a12,2)*(a11 + a22)*a2222)*pow(twtime[i],3))/(2.*pow(pow(a12,2) - a11*a22,4));
                
        h1221 = (-6*(2*pow(a11,2)*a12*(a1122*a12 - a1112*a22) + pow(a11,3)*(-(a12*a1222) + a1122*a22) + a11*a12*(-(pow(a12,2)*(a1112 + a1222)) + (a1111 + 2*a1122)*a12*a22 - a1112*pow(a22,2)) + pow(a12,2)*(a1122*pow(a12,2) + a22*(-2*a1112*a12 + a1111*a22)))*pow(twtime[i],3))/pow(pow(a12,2) - a11*a22,4);
             
        h1222 = (-3*a12*(a11 + a22)*pow(pow(a12,2) - a11*a22,2)*twtime[i] + 12*(2*pow(a11,2)*a12*(-(a12*a1222) + a1122*a22) - pow(a12,2)*(pow(a12,2)*a1222 - 2*a1122*a12*a22 + a1112*pow(a22,2)) + pow(a11,3)*(-(a1222*a22) + a12*a2222) + a11*a12*(a1122*(pow(a12,2) + pow(a22,2)) + a12*(-((a1112 + 2*a1222)*a22) + a12*a2222)))*pow(twtime[i],3))/(2.*pow(pow(a12,2) - a11*a22,4));
                
        h2111 = (-3*a12*(a11 + a22)*pow(pow(a12,2) - a11*a22,2)*twtime[i] + 12*(pow(a11,2)*a12*(-(a12*a1222) + a1122*a22) + a11*(-(pow(a12,2)*(2*a1112 + a1222)*a22) - a1112*pow(a22,3) + 2*a1122*a12*(pow(a12,2) + pow(a22,2))) + a12*((a1111 + a1122)*pow(a12,2)*a22 + a1111*pow(a22,3) - a1112*(pow(a12,3) + 2*a12*pow(a22,2))))*pow(twtime[i],3))/(2.*pow(pow(a12,2) - a11*a22,4));
              
        h2112 = (-6*(a1122*(pow(a12,4) + a11*pow(a22,3) + 2*pow(a12,2)*a22*(a11 + a22)) - a12*(pow(a11,2)*a1222*a22 + pow(a12,2)*(a1112 + a1222)*a22 + a1112*pow(a22,3) + 2*a11*a1222*(pow(a12,2) + pow(a22,2))) + a11*pow(a12,2)*(a11 + a22)*a2222)*pow(twtime[i],3))/pow(pow(a12,2) - a11*a22,4);
               
        h2121 = (3*(pow(a11,2) + pow(a12,2))*pow(pow(a12,2) - a11*a22,2)*twtime[i] + 12*(2*pow(a11,2)*a12*(-(a1122*a12) + a1112*a22) + pow(a11,3)*(a12*a1222 - a1122*a22) + a11*a12*(pow(a12,2)*(a1112 + a1222) - (a1111 + 2*a1122)*a12*a22 + a1112*pow(a22,2)) - pow(a12,2)*(a1122*pow(a12,2) + a22*(-2*a1112*a12 + a1111*a22)))*pow(twtime[i],3))/(2.*pow(pow(a12,2) - a11*a22,4));
              
        h2122 = (6*(2*pow(a11,2)*a12*(-(a12*a1222) + a1122*a22) - pow(a12,2)*(pow(a12,2)*a1222 - 2*a1122*a12*a22 + a1112*pow(a22,2)) + pow(a11,3)*(-(a1222*a22) + a12*a2222) + a11*a12*(a1122*(pow(a12,2) + pow(a22,2)) + a12*(-((a1112 + 2*a1222)*a22) + a12*a2222)))*pow(twtime[i],3))/pow(pow(a12,2) - a11*a22,4);
               
        h2211 = (6*(pow(a11,3)*a12*a1222 + a11*a12*(3*a1112*pow(a12,2) - 2*a1122*a12*a22 + 2*a1112*pow(a22,2)) - pow(a12,2)*(-(a1112*a12*a22) + a1111*(pow(a12,2) + pow(a22,2))) + pow(a11,2)*(a12*a1222*a22 - a1122*(3*pow(a12,2) + pow(a22,2))))*pow(twtime[i],3))/pow(pow(a12,2) - a11*a22,4);
           
        h2212 = (-3*a12*(a11 + a22)*pow(pow(a12,2) - a11*a22,2)*twtime[i] + 12*(a11*a12*(3*a1122*pow(a12,2) - 2*a12*a1222*a22 + 2*a1122*pow(a22,2)) - pow(a12,2)*(-(a1122*a12*a22) + a1112*(pow(a12,2) + pow(a22,2))) + pow(a11,3)*a12*a2222 + pow(a11,2)*(-(a1222*(3*pow(a12,2) + pow(a22,2))) + a12*a22*a2222))*pow(twtime[i],3))/(2.*pow(pow(a12,2) - a11*a22,4));
                
        h2221 = (-6*(-3*pow(a11,3)*a1122*a12 + pow(a11,4)*a1222 + pow(a12,3)*(a1112*a12 - a1111*a22) - a11*pow(a12,2)*(a1111*a12 + 2*a1122*a12 - 2*a1112*a22) + pow(a11,2)*a12*(3*a1112*a12 + a12*a1222 - a1122*a22))*pow(twtime[i],3))/pow(pow(a12,2) - a11*a22,4);
               
        h2222 = (3*(pow(a11,2) + pow(a12,2))*pow(pow(a12,2) - a11*a22,2)*twtime[i] - 12*(-3*pow(a11,3)*a12*a1222 + pow(a12,3)*(a1122*a12 - a1112*a22) - a11*pow(a12,2)*(a1112*a12 + 2*a12*a1222 - 2*a1122*a22) + pow(a11,4)*a2222 + pow(a11,2)*a12*(3*a1122*a12 - a1222*a22 + a12*a2222))*pow(twtime[i],3))/(2.*pow(pow(a12,2) - a11*a22,4));
        
    
    if (i != 0 ){

        
        /* Subtract to get interval parameters */
        h1111int = h1111 - h1111old;
        h1112int = h1112 - h1112old;
        h1121int = h1121 - h1121old;
        h1122int = h1122 - h1122old;
        h1211int = h1211 - h1211old;
        h1212int = h1212 - h1212old;
        h1221int = h1221 - h1221old;
        h1222int = h1222 - h1222old;
        h2111int = h2111 - h2111old;
        h2112int = h2112 - h2112old;
        h2121int = h2121 - h2121old;
        h2122int = h2122 - h2122old;
        h2211int = h2211 - h2211old;
        h2212int = h2212 - h2212old;
        h2221int = h2221 - h2221old;
        h2222int = h2222 - h2222old;
        
        
        /* Computation of interval Aijkl */
        a11 = inverted[i][0]; a12 = inverted[i][1]; a22 = inverted[i][2]; // From gendix
        a1111 = -(2*pow(a11,4)*h1111int + 8*pow(a11,3)*a12*h1112int + 8*a11*pow(a12,3)*h1222int + 2*pow(a12,4)*h2222int - 3*pow(a11,2)*(-4*pow(a12,2)*h1122int + (twtime[i]-twtime[i-1])))/(12.*pow((twtime[i]-twtime[i-1]),3));

        a1112 = -(2*pow(a11,3)*(a12*h1111int + a22*h1112int) + 6*pow(a11,2)*a12*(a12*h1112int + a22*h1122int) + 2*pow(a12,3)*(a12*h1222int + a22*h2222int) - 3*a11*a12*(-2*a12*(a12*h1122int + a22*h1222int) + (twtime[i]-twtime[i-1])))/(12.*pow((twtime[i]-twtime[i-1]),3));

        a1121 = -(2*pow(a11,3)*(a12*h1111int + a22*h1112int) + 6*pow(a11,2)*a12*(a12*h1112int + a22*h1122int) + 2*pow(a12,3)*(a12*h1222int + a22*h2222int) - 3*a11*a12*(-2*a12*(a12*h1122int + a22*h1222int) + (twtime[i]-twtime[i-1])))/(12.*pow((twtime[i]-twtime[i-1]),3));

        a1122 = (-2*(pow(a11,2)*(pow(a12,2)*h1111int + 2*a12*a22*h1112int + pow(a22,2)*h1122int) + 2*a11*a12*(pow(a12,2)*h1112int + 2*a12*a22*h1122int + pow(a22,2)*h1222int) + pow(a12,2)*(pow(a12,2)*h1122int + 2*a12*a22*h1222int + pow(a22,2)*h2222int)) + 3*a11*a22*(twtime[i]-twtime[i-1]))/(12.*pow((twtime[i]-twtime[i-1]),3));

        a1211 = -(2*pow(a11,3)*(a12*h1111int + a22*h1112int) + 6*pow(a11,2)*a12*(a12*h1112int + a22*h1122int) + 2*pow(a12,3)*(a12*h1222int + a22*h2222int) - 3*a11*a12*(-2*a12*(a12*h1122int + a22*h1222int) + (twtime[i]-twtime[i-1])))/(12.*pow((twtime[i]-twtime[i-1]),3));

        a1212 = -(2*pow(a11,2)*(pow(a12,2)*h1111int + 2*a12*a22*h1112int + pow(a22,2)*h1122int) + 4*a11*a12*(pow(a12,2)*h1112int + 2*a12*a22*h1122int + pow(a22,2)*h1222int) + pow(a12,2)*(2*pow(a12,2)*h1122int + 4*a12*a22*h1222int + 2*pow(a22,2)*h2222int - 3*(twtime[i]-twtime[i-1])))/(12.*pow((twtime[i]-twtime[i-1]),3));

        a1221 = -(2*pow(a11,2)*(pow(a12,2)*h1111int + 2*a12*a22*h1112int + pow(a22,2)*h1122int) + 4*a11*a12*(pow(a12,2)*h1112int + 2*a12*a22*h1122int + pow(a22,2)*h1222int) + pow(a12,2)*(2*pow(a12,2)*h1122int + 4*a12*a22*h1222int + 2*pow(a22,2)*h2222int - 3*(twtime[i]-twtime[i-1])))/(12.*pow((twtime[i]-twtime[i-1]),3));

        a1222 = -(2*a11*(pow(a12,3)*h1111int + 3*pow(a12,2)*a22*h1112int + 3*a12*pow(a22,2)*h1122int + pow(a22,3)*h1222int) + a12*(2*pow(a12,3)*h1112int + 6*pow(a12,2)*a22*h1122int + 6*a12*pow(a22,2)*h1222int + 2*pow(a22,3)*h2222int - 3*a22*(twtime[i]-twtime[i-1])))/(12.*pow((twtime[i]-twtime[i-1]),3));

        a2111 = -(2*pow(a11,3)*(a12*h1111int + a22*h1112int) + 6*pow(a11,2)*a12*(a12*h1112int + a22*h1122int) + 2*pow(a12,3)*(a12*h1222int + a22*h2222int) - 3*a11*a12*(-2*a12*(a12*h1122int + a22*h1222int) + (twtime[i]-twtime[i-1])))/(12.*pow((twtime[i]-twtime[i-1]),3));

        a2112 = -(2*pow(a11,2)*(pow(a12,2)*h1111int + 2*a12*a22*h1112int + pow(a22,2)*h1122int) + 4*a11*a12*(pow(a12,2)*h1112int + 2*a12*a22*h1122int + pow(a22,2)*h1222int) + pow(a12,2)*(2*pow(a12,2)*h1122int + 4*a12*a22*h1222int + 2*pow(a22,2)*h2222int - 3*(twtime[i]-twtime[i-1])))/(12.*pow((twtime[i]-twtime[i-1]),3));

        a2121 = -(2*pow(a11,2)*(pow(a12,2)*h1111int + 2*a12*a22*h1112int + pow(a22,2)*h1122int) + 4*a11*a12*(pow(a12,2)*h1112int + 2*a12*a22*h1122int + pow(a22,2)*h1222int) + pow(a12,2)*(2*pow(a12,2)*h1122int + 4*a12*a22*h1222int + 2*pow(a22,2)*h2222int - 3*(twtime[i]-twtime[i-1])))/(12.*pow((twtime[i]-twtime[i-1]),3));

        a2122 = -(2*a11*(pow(a12,3)*h1111int + 3*pow(a12,2)*a22*h1112int + 3*a12*pow(a22,2)*h1122int + pow(a22,3)*h1222int) + a12*(2*pow(a12,3)*h1112int + 6*pow(a12,2)*a22*h1122int + 6*a12*pow(a22,2)*h1222int + 2*pow(a22,3)*h2222int - 3*a22*(twtime[i]-twtime[i-1])))/(12.*pow((twtime[i]-twtime[i-1]),3));

        a2211 = (-2*(pow(a11,2)*(pow(a12,2)*h1111int + 2*a12*a22*h1112int + pow(a22,2)*h1122int) + 2*a11*a12*(pow(a12,2)*h1112int + 2*a12*a22*h1122int + pow(a22,2)*h1222int) + pow(a12,2)*(pow(a12,2)*h1122int + 2*a12*a22*h1222int + pow(a22,2)*h2222int)) + 3*a11*a22*(twtime[i]-twtime[i-1]))/(12.*pow((twtime[i]-twtime[i-1]),3));

        a2212 = -(2*a11*(pow(a12,3)*h1111int + 3*pow(a12,2)*a22*h1112int + 3*a12*pow(a22,2)*h1122int + pow(a22,3)*h1222int) + a12*(2*pow(a12,3)*h1112int + 6*pow(a12,2)*a22*h1122int + 6*a12*pow(a22,2)*h1222int + 2*pow(a22,3)*h2222int - 3*a22*(twtime[i]-twtime[i-1])))/(12.*pow((twtime[i]-twtime[i-1]),3));

        a2221 = -(2*a11*(pow(a12,3)*h1111int + 3*pow(a12,2)*a22*h1112int + 3*a12*pow(a22,2)*h1122int + pow(a22,3)*h1222int) + a12*(2*pow(a12,3)*h1112int + 6*pow(a12,2)*a22*h1122int + 6*a12*pow(a22,2)*h1222int + 2*pow(a22,3)*h2222int - 3*a22*(twtime[i]-twtime[i-1])))/(12.*pow((twtime[i]-twtime[i-1]),3));

        a2222 = (-2*(pow(a12,4)*h1111int + 4*pow(a12,3)*a22*h1112int + 6*pow(a12,2)*pow(a22,2)*h1122int + 4*a12*pow(a22,3)*h1222int + pow(a22,4)*h2222int) + 3*pow(a22,2)*(twtime[i]-twtime[i-1]))/(12.*pow((twtime[i]-twtime[i-1]),3));
        
        
        if (abs((a1112+a1121+a1211+a2111)/4-a1112) > 1e-5 ) sf_warning("Quartic: interval A is not symmetric at x^3 y : %f, %f, %f, %f",a1112, a1121,a1211,a2111);
        if (abs((a1122+a1221+a2211+a2112+a1212+a2121)/6-a1122) > 1e-5 ) sf_warning("Quartic: interval A is not symmetric at x^2 y^2 : %f, %f, %f, %f, %f, %f",a1122, a1221,a2211,a2112,a1212,a2121);
        if (abs((a2221+a2212+a2122+a1222)/4-a1222) > 1e-5 ) sf_warning("Quartic: interval A is not symmetric at x y^3 : %f, %f, %f, %f",a2221, a2212,a2122,a1222);
        
        inverted[i][3] = a1111;
        inverted[i][4] = a1112;
        inverted[i][5] = a1122;
        inverted[i][6] = a1222;
        inverted[i][7] = a2222;
        
        /* Store eff values upto i*/
        h1111old = h1111;
        h1112old = h1112;
        h1121old = h1121;
        h1122old = h1122;
        h1211old = h1211;
        h1212old = h1212;
        h1221old = h1221;
        h1222old = h1222;
        h2111old = h2111;
        h2112old = h2112;
        h2121old = h2121;
        h2122old = h2122;
        h2211old = h2211;
        h2212old = h2212;
        h2221old = h2221;
        h2222old = h2222;
        
    } else { /* First layer */

        h1111old = h1111;
        h1112old = h1112;
        h1121old = h1121;
        h1122old = h1122;
        h1211old = h1211;
        h1212old = h1212;
        h1221old = h1221;
        h1222old = h1222;
        h2111old = h2111;
        h2112old = h2112;
        h2121old = h2121;
        h2122old = h2122;
        h2211old = h2211;
        h2212old = h2212;
        h2221old = h2221;
        h2222old = h2222;

        inverted[i][3] = a1111;
        inverted[i][4] = a1112;
        inverted[i][5] = a1122;
        inverted[i][6] = a1222;
        inverted[i][7] = a2222;

    }

}

}

