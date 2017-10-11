/* Compute bond transformation of the Cij. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include <math.h>
#include <rsf.h>
#include "bond.h"

/*Follow Lapilli and Fowler (2013) of the first two rotations*/
void bond ( float *phi  /* Rz around z (1st)*/, 
		    float *the  /* Ry'around y' (2nd)*/,
		    float *c11o, 
		    float *c12o, 
		    float *c13o, 
		    float *c14o, 
		    float *c15o, 
		    float *c16o, 
		    float *c22o, 
		    float *c23o, 
		    float *c24o, 
		    float *c25o, 
		    float *c26o, 
		    float *c33o, 
		    float *c34o, 
		    float *c35o, 
		    float *c36o, 
		    float *c44o, 
		    float *c45o, 
		    float *c46o, 
		    float *c55o, 
		    float *c56o, 
		    float *c66o)
/*< Bond rotations >*/
{
    float c11r,c12r,c13r,c14r,c15r,c16r,c22r,c23r,c24r,c25r,c26r,c33r,c34r,c35r,c36r,c44r,c45r,c46r,c55r,c56r,c66r;
    float c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66;
    float s1 = sin(*phi); float c1 = cos(*phi);
    float s2 = sin(*the); float c2 = cos(*the);
    
    c11=*c11o; c12=*c12o; c13=*c13o;
    c14=*c14o; c15=*c15o; c16=*c16o;
    c22=*c22o; c23=*c23o; c24=*c24o;
    c25=*c25o; c26=*c26o; c33=*c33o;
    c34=*c34o; c35=*c35o; c36=*c36o;
    c44=*c44o; c45=*c45o; c46=*c46o;
    c55=*c55o; c56=*c56o; c66=*c66o;
    
    c11r = c22*pow(s1,4) + 4*c1*pow(s1,3)*(c2*c26 - c24*s2) + 2*pow(c1,2)*pow(s1,2)*(c12*pow(c2,2) + 2*pow(c2,2)*c66 - 2*c2*(c25 + 2*c46)*s2 + (c23 + 2*c44)*pow(s2,2)) + 4*pow(c1,3)*s1*(c16*pow(c2,3) - s2*(c14*pow(c2,2) + 2*pow(c2,2)*c56 - c2*(c36 + 2*c45)*s2 + c34*pow(s2,2))) + pow(c1,4)*(c11*pow(c2,4) + s2*(-4*c15*pow(c2,3) + s2*(2*c13*pow(c2,2) + 4*pow(c2,2)*c55 - 4*c2*c35*s2 + c33*pow(s2,2))));
   
    c12r = pow(c1,2)*pow(s1,2)*(c11*pow(c2,4) + c22 - 4*pow(c2,2)*c66 - 4*c15*pow(c2,3)*s2 + 8*c2*c46*s2 + 2*c13*pow(c2,2)*pow(s2,2) - 4*c44*pow(s2,2) + 4*pow(c2,2)*c55*pow(s2,2) - 4*c2*c35*pow(s2,3) + c33*pow(s2,4)) + pow(c1,4)*(c12*pow(c2,2) + s2*(-2*c2*c25 + c23*s2)) + pow(s1,4)*(c12*pow(c2,2) + s2*(-2*c2*c25 + c23*s2)) - 2*pow(c1,3)*s1*(c16*pow(c2,3) - pow(c2,2)*(c14 + 2*c56)*s2 + s2*(c24 - c34*pow(s2,2)) + c2*(-c26 + (c36 + 2*c45)*pow(s2,2))) + 2*c1*pow(s1,3)*(c16*pow(c2,3) - pow(c2,2)*(c14 + 2*c56)*s2 + s2*(c24 - c34*pow(s2,2)) + c2*(-c26 + (c36 + 2*c45)*pow(s2,2)));

    c13r = pow(s1,2)*(pow(c2,2)*c23 + 2*c2*c25*s2 + c12*pow(s2,2)) + 2*c1*s1*(pow(c2,3)*c36 - pow(c2,2)*(c34 - 2*c56)*s2 + c2*(c16 - 2*c45)*pow(s2,2) - c14*pow(s2,3)) + pow(c1,2)*(c13*(pow(c2,4) + pow(s2,4)) + c2*s2*(-2*pow(c2,2)*c35 + c2*(c11 + c33 - 4*c55)*s2 + 2*c35*pow(s2,2) + 2*c15*(pow(c2,2) - pow(s2,2))));
    
    
    c14r = pow(s1,3)*(-(pow(c2,2)*c25) + c2*(-c12 + c23)*s2 + c25*pow(s2,2)) + pow(c1,3)*(c14*pow(c2,3) + s2*(c16*pow(c2,2) - 2*pow(c2,2)*c45 + c2*(c34 - 2*c56)*s2 + c36*pow(s2,2))) + c1*pow(s1,2)*(-2*pow(c2,3)*c56 + 2*pow(c2,2)*(-c16 + c36 + c45)*s2 + s2*(c26 - 2*c45*pow(s2,2)) + c2*(c24 + 2*(c14 - c34 + c56)*pow(s2,2))) + pow(c1,2)*s1*(pow(c2,3)*(-c11 + c13 + 2*c55)*s2 - 2*c46*pow(s2,2) + c35*pow(s2,4) - c15*(pow(c2,4) - 3*pow(c2,2)*pow(s2,2)) + pow(c2,2)*(2*c46 - 3*c35*pow(s2,2)) + c2*(-2*c44*s2 + 2*c66*s2 + (-c13 + c33 - 2*c55)*pow(s2,3)));
    
    c15r = pow(s1,3)*(c2*c24 + c26*s2) + c1*pow(s1,2)*(pow(c2,2)*(c25 + 2*c46) + c2*(c12 - c23 - 2*c44 + 2*c66)*s2 - (c25 + 2*c46)*pow(s2,2)) + pow(c1,2)*s1*(2*pow(c2,3)*c56 + pow(c2,2)*(3*c16 - 2*c36 - 4*c45)*s2 + c2*(3*c34 - 4*c56)*pow(s2,2) + (c36 + 2*c45)*pow(s2,3) + c14*(pow(c2,3) - 2*c2*pow(s2,2))) + pow(c1,3)*(c15*(pow(c2,4) - 3*pow(c2,2)*pow(s2,2)) + s2*(c11*pow(c2,3) - c13*pow(c2,3) - 2*pow(c2,3)*c55 + 3*pow(c2,2)*c35*s2 + c13*c2*pow(s2,2) - c2*c33*pow(s2,2) + 2*c2*c55*pow(s2,2) - c35*pow(s2,3)));
    
    c16r = c1*(pow(c1,2)*s1*(c12*pow(c2,2) - c11*pow(c2,4) + 2*pow(c2,2)*c66 + 4*c15*pow(c2,3)*s2 - 2*c2*c25*s2 - 4*c2*c46*s2 - 2*c13*pow(c2,2)*pow(s2,2) + c23*pow(s2,2) + 2*c44*pow(s2,2) - 4*pow(c2,2)*c55*pow(s2,2) + 4*c2*c35*pow(s2,3) - c33*pow(s2,4)) + pow(s1,3)*(-(c12*pow(c2,2)) + c22 + s2*(2*c2*c25 - c23*s2)) + pow(c1,3)*(c16*pow(c2,3) - s2*(c14*pow(c2,2) + 2*pow(c2,2)*c56 - c2*(c36 + 2*c45)*s2 + c34*pow(s2,2))) + c1*pow(s1,2)*(-2*c16*pow(c2,3) - 3*c24*s2 + 2*pow(c2,2)*(c14 + 2*c56)*s2 + 2*c34*pow(s2,3) + c2*(3*c26 - 2*(c36 + 2*c45)*pow(s2,2))));
    
    c22r = pow(c1,4)*c22 + pow(c1,3)*(-4*c2*c26*s1 + 4*c24*s1*s2) + 2*pow(c1,2)*pow(s1,2)*(c12*pow(c2,2) + 2*pow(c2,2)*c66 - 2*c2*(c25 + 2*c46)*s2 + (c23 + 2*c44)*pow(s2,2)) - 4*c1*pow(s1,3)*(c16*pow(c2,3) - s2*(c14*pow(c2,2) + 2*pow(c2,2)*c56 - c2*(c36 + 2*c45)*s2 + c34*pow(s2,2))) + pow(s1,4)*(c11*pow(c2,4) + s2*(-4*c15*pow(c2,3) + s2*(2*c13*pow(c2,2) + 4*pow(c2,2)*c55 - 4*c2*c35*s2 + c33*pow(s2,2))));
    
    c23r = pow(c1,2)*(pow(c2,2)*c23 + 2*c2*c25*s2 + c12*pow(s2,2)) - 2*c1*s1*(pow(c2,3)*c36 - pow(c2,2)*(c34 - 2*c56)*s2 + c2*(c16 - 2*c45)*pow(s2,2) - c14*pow(s2,3)) + pow(s1,2)*(c13*(pow(c2,4) + pow(s2,4)) + c2*s2*(-2*pow(c2,2)*c35 + c2*(c11 + c33 - 4*c55)*s2 + 2*c35*pow(s2,2) + 2*c15*(pow(c2,2) - pow(s2,2))));
    
    c24r = pow(c1,3)*(c2*c24 + c26*s2) + pow(c1,2)*s1*(-(pow(c2,2)*(c25 + 2*c46)) + c2*(-c12 + c23 + 2*c44 - 2*c66)*s2 + (c25 + 2*c46)*pow(s2,2)) + c1*pow(s1,2)*(2*pow(c2,3)*c56 + pow(c2,2)*(3*c16 - 2*c36 - 4*c45)*s2 + c2*(3*c34 - 4*c56)*pow(s2,2) + (c36 + 2*c45)*pow(s2,3) + c14*(pow(c2,3) - 2*c2*pow(s2,2))) + pow(s1,3)*(-(c15*(pow(c2,4) - 3*pow(c2,2)*pow(s2,2))) + 
s2*(-(c11*pow(c2,3)) + c13*pow(c2,3) + 2*pow(c2,3)*c55 - 3*pow(c2,2)*c35*s2 - c13*c2*pow(s2,2) + c2*c33*pow(s2,2) - 2*c2*c55*pow(s2,2) + c35*pow(s2,3)));
    
    c25r = pow(c1,3)*(pow(c2,2)*c25 + c2*(c12 - c23)*s2 - c25*pow(s2,2)) + pow(s1,3)*(c14*pow(c2,3) + s2*(c16*pow(c2,2) - 2*pow(c2,2)*c45 + c2*(c34 - 2*c56)*s2 + c36*pow(s2,2))) + c1*pow(s1,2)*(pow(c2,3)*(c11 - c13 - 2*c55)*s2 + 2*c46*pow(s2,2) - c35*pow(s2,4) + c15*(pow(c2,4) - 3*pow(c2,2)*pow(s2,2)) + pow(c2,2)*(-2*c46 + 3*c35*pow(s2,2)) + c2*s2*(2*c44 - 2*c66 + (c13 - c33 + 2*c55)*pow(s2,2))) + pow(c1,2)*s1*(-2*pow(c2,3)*c56 + 2*pow(c2,2)*(-c16 + c36 + c45)*s2 + s2*(c26 - 2*c45*pow(s2,2)) + c2*(c24 + 2*(c14 - c34 + c56)*pow(s2,2)));
    
    c26r = c1*(pow(c1,3)*(c2*c26 - c24*s2) - pow(c1,2)*s1*(c12*pow(c2,2) - c22 + 2*pow(c2,2)*c66 - 2*c2*c25*s2 - 4*c2*c46*s2 + c23*pow(s2,2) + 2*c44*pow(s2,2)) + c1*pow(s1,2)*(3*c16*pow(c2,3) + 2*c24*s2 - 3*pow(c2,2)*(c14 + 2*c56)*s2 - 3*c34*pow(s2,3) + c2*(-2*c26 + 3*(c36 + 2*c45)*pow(s2,2))) + pow(s1,3)*(c12*pow(c2,2) - c11*pow(c2,4) + s2*(4*c15*pow(c2,3) - 2*pow(c2,2)*(c13 + 2*c55)*s2 + s2*(c23 - c33*pow(s2,2)) - 2*c2*(c25 - 2*c35*pow(s2,2)))));
    
    c33r = pow(c2,4)*c33 + 4*pow(c2,3)*c35*s2 + 2*pow(c2,2)*(c13 + 2*c55)*pow(s2,2) + 4*c15*c2*pow(s2,3) + c11*pow(s2,4);
    
    c34r = c1*(pow(c2,3)*c34 + pow(c2,2)*(c36 + 2*c45)*s2 + c2*(c14 + 2*c56)*pow(s2,2) + c16*pow(s2,3)) + s1*(-(pow(c2,4)*c35) + pow(c2,3)*(-c13 + c33 - 2*c55)*s2 + 3*pow(c2,2)*(-c15 + c35)*pow(s2,2) + c2*(-c11 + c13 + 2*c55)*pow(s2,3) + c15*pow(s2,4));
    
    c35r = s1*(pow(c2,3)*c34 + pow(c2,2)*(c36 + 2*c45)*s2 + c2*(c14 + 2*c56)*pow(s2,2) + c16*pow(s2,3)) + c1*(pow(c2,4)*c35 + pow(c2,3)*(c13 - c33 + 2*c55)*s2 + 3*pow(c2,2)*(c15 - c35)*pow(s2,2) + c2*(c11 - c13 - 2*c55)*pow(s2,3) - c15*pow(s2,4));
    
    c36r = c1*(s1*(pow(c2,2)*c23 + 2*c2*c25*s2 + c12*pow(s2,2)) + c1*(pow(c2,3)*c36 - pow(c2,2)*(c34 - 2*c56)*s2 + c2*(c16 - 2*c45)*pow(s2,2) - c14*pow(s2,3)) - s1*(c13*(pow(c2,4) + pow(s2,4)) + c2*s2*(-2*pow(c2,2)*c35 + c2*(c11 + c33 - 4*c55)*s2 + 2*c35*pow(s2,2) + 2*c15*(pow(c2,2) - pow(s2,2)))));
    
    c44r = pow(c1,2)*(pow(c2,2)*c44 + 2*c2*c46*s2 + c66*pow(s2,2)) - 2*c1*s1*(pow(c2,3)*c45 + pow(c2,2)*(c14 - c34 + c56)*s2 + c2*(c16 - c36 - c45)*pow(s2,2) - c56*pow(s2,3)) + pow(s1,2)*(pow(c2,4)*c55 + 2*pow(c2,3)*(c15 - c35)*s2 + pow(c2,2)*(c11 - 2*c13 + c33 - 2*c55)*pow(s2,2) + 2*c2*(-c15 + c35)*pow(s2,3) + c55*pow(s2,4));
    
    c45r = pow(c1,2)*(pow(c2,3)*c45 + pow(c2,2)*(c14 - c34 + c56)*s2 + c2*(c16 - c36 - c45)*pow(s2,2) - c56*pow(s2,3)) + pow(s1,2)*(-(pow(c2,3)*c45) - pow(c2,2)*(c14 - c34 + c56)*s2 + c2*(-c16 + c36 + c45)*pow(s2,2) + c56*pow(s2,3)) + c1*s1*(-(pow(c2,4)*c55) + 2*pow(c2,3)*(-c15 + c35)*s2 + 2*c2*s2*(c46 + (c15 - c35)*pow(s2,2)) + pow(c2,2)*(c44 - (c11 - 2*c13 + c33 - 2*c55)*pow(s2,2)) + pow(s2,2)*(c66 - c55*pow(s2,2)));
    
    c46r = c1*(pow(c1,2)*(pow(c2,2)*c46 + c2*(-c44 + c66)*s2 - c46*pow(s2,2)) + pow(s1,2)*(pow(c2,3)*(c11 - c13 - 2*c55)*s2 + c15*(pow(c2,4) - 3*pow(c2,2)*pow(s2,2)) - pow(c2,2)*(c25 - 3*c35*pow(s2,2)) + pow(s2,2)*(c25 - c35*pow(s2,2)) + c2*s2*(-c12 + c23 + (c13 - c33 + 2*c55)*pow(s2,2))) + c1*s1*(-(pow(c2,3)*c56) + pow(c2,2)*(-2*c16 + c36 + 3*c45)*s2 + c14*(-pow(c2,3) + c2*pow(s2,2)) + s2*(c26 - (c36 + c45)*pow(s2,2)) + c2*(c24 + (-2*c34 + 3*c56)*pow(s2,2))));
    
    c55r = pow(s1,2)*(pow(c2,2)*c44 + 2*c2*c46*s2 + c66*pow(s2,2)) + 2*c1*s1*(pow(c2,3)*c45 + pow(c2,2)*(c14 - c34 + c56)*s2 + c2*(c16 - c36 - c45)*pow(s2,2) - c56*pow(s2,3)) + pow(c1,2)*(pow(c2,4)*c55 + 2*pow(c2,3)*(c15 - c35)*s2 + pow(c2,2)*(c11 - 2*c13 + c33 - 2*c55)*pow(s2,2) + 2*c2*(-c15 + c35)*pow(s2,3) + c55*pow(s2,4));
    
    c56r = c1*(pow(c1,2)*(pow(c2,3)*c56 + pow(c2,2)*(c16 - c36 - c45)*s2 - c2*(c14 - c34 + c56)*pow(s2,2) + c45*pow(s2,3)) + c1*s1*(pow(c2,3)*(-c11 + c13 + 2*c55)*s2 - c15*(pow(c2,4) - 3*pow(c2,2)*pow(s2,2)) + pow(c2,2)*(c25 + c46 - 3*c35*pow(s2,2)) + pow(s2,2)*(-c25 - c46 + c35*pow(s2,2)) + c2*s2*(c12 - c23 - c44 + c66 - c13*pow(s2,2) + c33*pow(s2,2) - 2*c55*pow(s2,2))) - pow(s1,2)*(c14*pow(c2,3) - c26*s2 + pow(c2,2)*(c16 - 2*c45)*s2 + c36*pow(s2,3) - c2*(c24 - (c34 - 2*c56)*pow(s2,2))));
    
    c66r = pow(c1,2)*(pow(c1,2)*(pow(c2,2)*c66 - 2*c2*c46*s2 + c44*pow(s2,2)) + pow(s1,2)*(-2*c12*pow(c2,2) + c11*pow(c2,4) + c22 - 4*c15*pow(c2,3)*s2 + 4*c2*c25*s2 + 2*c13*pow(c2,2)*pow(s2,2) - 2*c23*pow(s2,2) + 4*pow(c2,2)*c55*pow(s2,2) - 4*c2*c35*pow(s2,3) + c33*pow(s2,4)) - 2*c1*s1*(c16*pow(c2,3) - pow(c2,2)*(c14 + 2*c56)*s2 + s2*(c24 - c34*pow(s2,2)) + c2*(-c26 + (c36 + 2*c45)*pow(s2,2))));
    
    *c11o=c11r; *c12o=c12r; *c13o=c13r;
    *c14o=c14r; *c15o=c15r; *c16o=c16r;
    *c22o=c22r; *c23o=c23r; *c24o=c24r;
    *c25o=c25r; *c26o=c26r; *c33o=c33r;
    *c34o=c34r; *c35o=c35r; *c36o=c36r;
    *c44o=c44r; *c45o=c45r; *c46o=c46r;
    *c55o=c55r; *c56o=c56r; *c66o=c66r;
    
   return;
}
