/*Computation for orthorhombic group velocity for qP in layered media*/
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


float xPvelo(float c11, float c22,float c33,float c44,float c55,float c66,float c12,float c13,float c23, float p1, float p2) {
/*<Compute qP-wave group velocity (x) given cij and px (p1) and py (p2) for layered orthorhombic media>*/

    float complex p3, dx, dt;


/*p3 ###################################################################################################################*/

    p3 = (c33*c44 + c33*c55 + c44*c55 + cpow(c13,2)*c44*cpow(p1,2) - c11*c33*c44*cpow(p1,2) + 2*c13*c44*c55*cpow(p1,2) - c33*c55*c66*cpow(p1,2) + 
	  cpow(c23,2)*c55*cpow(p2,2) - c22*c33*c55*cpow(p2,2) + 2*c23*c44*c55*cpow(p2,2) - c33*c44*c66*cpow(p2,2))/(3.*c33*c44*c55) + 
	(cpow(2,0.3333333333333333)*(-cpow(c33*c44 + c33*c55 + c44*c55 + cpow(c13,2)*c44*cpow(p1,2) - c11*c33*c44*cpow(p1,2) + 2*c13*c44*c55*cpow(p1,2) - 
					   c33*c55*c66*cpow(p1,2) + cpow(c23,2)*c55*cpow(p2,2) - c22*c33*c55*cpow(p2,2) + 2*c23*c44*c55*cpow(p2,2) - c33*c44*c66*cpow(p2,2),2) - 
				     3*c33*c44*c55*(-c33 - c44 - c55 - cpow(c13,2)*cpow(p1,2) + c11*c33*cpow(p1,2) + c11*c44*cpow(p1,2) - 2*c13*c55*cpow(p1,2) + c44*c55*cpow(p1,2) + 
						    c33*c66*cpow(p1,2) + c55*c66*cpow(p1,2) - c11*c44*c55*cpow(p1,4) + cpow(c13,2)*c66*cpow(p1,4) - c11*c33*c66*cpow(p1,4) + 2*c13*c55*c66*cpow(p1,4) - 
						    cpow(c23,2)*cpow(p2,2) + c22*c33*cpow(p2,2) - 2*c23*c44*cpow(p2,2) + c22*c55*cpow(p2,2) + c44*c55*cpow(p2,2) + c33*c66*cpow(p2,2) + 
						    c44*c66*cpow(p2,2) + cpow(c13,2)*c22*cpow(p1,2)*cpow(p2,2) - 2*c12*c13*c23*cpow(p1,2)*cpow(p2,2) + c11*cpow(c23,2)*cpow(p1,2)*cpow(p2,2) + 
						    cpow(c12,2)*c33*cpow(p1,2)*cpow(p2,2) - c11*c22*c33*cpow(p1,2)*cpow(p2,2) - 2*c12*c13*c44*cpow(p1,2)*cpow(p2,2) + 
						    2*c11*c23*c44*cpow(p1,2)*cpow(p2,2) + 2*c13*c22*c55*cpow(p1,2)*cpow(p2,2) - 2*c12*c23*c55*cpow(p1,2)*cpow(p2,2) - 
						    2*c12*c44*c55*cpow(p1,2)*cpow(p2,2) - 2*c13*c23*c66*cpow(p1,2)*cpow(p2,2) + 2*c12*c33*c66*cpow(p1,2)*cpow(p2,2) - 
						    2*c13*c44*c66*cpow(p1,2)*cpow(p2,2) - 2*c23*c55*c66*cpow(p1,2)*cpow(p2,2) - 4*c44*c55*c66*cpow(p1,2)*cpow(p2,2) - c22*c44*c55*cpow(p2,4) + 
						    cpow(c23,2)*c66*cpow(p2,4) - c22*c33*c66*cpow(p2,4) + 2*c23*c44*c66*cpow(p2,4))))/
	(3.*c33*c44*c55*cpow(-2*cpow(c33,3)*cpow(c44,3) + 3*cpow(c33,3)*cpow(c44,2)*c55 + 3*cpow(c33,2)*cpow(c44,3)*c55 + 3*cpow(c33,3)*c44*cpow(c55,2) - 
			     12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2) + 3*c33*cpow(c44,3)*cpow(c55,2) - 2*cpow(c33,3)*cpow(c55,3) + 3*cpow(c33,2)*c44*cpow(c55,3) + 
			     3*c33*cpow(c44,2)*cpow(c55,3) - 2*cpow(c44,3)*cpow(c55,3) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,2) + 
			     6*c11*cpow(c33,3)*cpow(c44,3)*cpow(p1,2) + 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2) - 6*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2) - 
			     3*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,2) - 6*c11*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) - 12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) + 
			     3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2) - 3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2) + 
			     6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 
			     12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) - 6*cpow(c13,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 
			     3*c11*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 6*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 9*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) + 
			     6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2) + 12*c13*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) + 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) - 
			     12*c13*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 3*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2) - 
			     6*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2) + 6*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2) - 
			     6*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2) - 3*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2) - 6*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,4) + 
			     12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,4) - 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*cpow(p1,4) + 
			     3*cpow(c13,4)*c33*cpow(c44,2)*c55*cpow(p1,4) - 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4) + 
			     3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4) - 6*cpow(c13,4)*cpow(c44,3)*c55*cpow(p1,4) + 3*c11*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,4) - 
			     24*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4) + 3*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 
			     24*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 
			     12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 24*cpow(c13,3)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
			     6*c11*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) - 33*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
			     18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 
			     18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 24*cpow(c13,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) + 
			     9*c11*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 
			     6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4) + 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4) - 
			     6*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4) + 6*c11*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4) - 
			     6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 
			     12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c13*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4) - 
			     12*c13*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) - 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) + 
			     3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4) - 6*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) + 
			     3*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) - 2*cpow(c13,6)*cpow(c44,3)*cpow(p1,6) + 6*c11*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,6) - 
			     6*cpow(c11,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,6) + 2*cpow(c11,3)*cpow(c33,3)*cpow(c44,3)*cpow(p1,6) - 
			     12*cpow(c13,5)*cpow(c44,3)*c55*cpow(p1,6) + 24*c11*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,6) - 
			     12*cpow(c11,2)*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,6) - 24*cpow(c13,4)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) + 
			     33*c11*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 9*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 
			     16*cpow(c13,3)*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) + 18*c11*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) - 
			     3*cpow(c13,4)*c33*cpow(c44,2)*c55*c66*cpow(p1,6) + 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,6) - 
			     3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,6) - 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) + 
			     12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) - 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 
			     18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) - 
			     3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) + 6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,6) + 
			     2*cpow(c33,3)*cpow(c55,3)*cpow(c66,3)*cpow(p1,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p2,2) - 
			     3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p2,2) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p2,2) + 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,2) - 
			     6*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,2) + 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 
			     12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 
			     12*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) + 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,2) + 
			     6*c22*cpow(c33,3)*cpow(c55,3)*cpow(p2,2) - 3*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,2) - 6*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 
			     12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 3*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 
			     6*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 9*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 12*c23*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) - 
			     9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) + 6*cpow(c33,3)*cpow(c44,3)*c66*cpow(p2,2) - 6*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,2) - 
			     6*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,2) - 3*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,2) - 
			     3*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,2) - 3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
			     6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
			     6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 
			     6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) + 
			     18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 
			     3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
			     18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     12*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
			     18*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     6*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
			     18*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     12*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     12*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     24*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     9*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c11*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
			     18*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
			     12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
			     24*c13*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 12*c13*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
			     18*c12*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
			     9*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
			     18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 48*c13*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
			     18*c12*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
			     18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 12*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 
			     12*c11*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
			     18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) - 
			     18*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 3*cpow(c13,2)*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
			     6*c11*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 42*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
			     18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 
			     18*c12*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
			     18*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 36*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
			     6*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 6*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
			     24*c13*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
			     12*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 12*c22*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
			     3*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 6*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
			     42*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 24*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
			     18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 36*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 
			     3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) + 
			     6*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,4)*cpow(c23,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
			     3*cpow(c13,4)*c22*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 18*c12*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
			     3*c11*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
			     9*cpow(c12,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
			     6*c11*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
			     3*cpow(c11,2)*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 9*c11*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
			     3*cpow(c11,2)*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,4)*c23*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 
			     18*c12*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 6*c11*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 
			     18*c11*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 6*cpow(c11,2)*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 
			     24*cpow(c13,3)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,3)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
			     54*c12*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 6*c11*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
			     18*cpow(c12,2)*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
			     12*c11*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
			     48*cpow(c13,3)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 54*c12*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
			     12*c11*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
			     24*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
			     36*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 9*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
			     27*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 18*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
			     48*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 36*c12*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
			     18*c11*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 6*cpow(c13,4)*c33*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) - 
			     12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 
			     18*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 
			     18*c11*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 18*c11*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
			     42*cpow(c13,3)*c33*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) - 42*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
			     3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
			     18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
			     9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
			     60*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 54*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
			     6*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 96*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
			     18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
			     12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
			     48*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 72*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
			     72*c13*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 3*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
			     3*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
			     18*c12*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 24*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
			     6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
			     30*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 36*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
			     6*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,3)*cpow(p1,4)*cpow(p2,2) + 3*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p2,4) - 
			     6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,4) + 3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,4) + 
			     12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) - 12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) + 
			     12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 
			     6*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,4) + 12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,4) - 6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*cpow(p2,4) - 
			     6*cpow(c23,4)*c44*cpow(c55,3)*cpow(p2,4) + 3*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,4) - 24*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,4) + 
			     3*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) + 24*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) - 
			     24*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 6*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 
			     33*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 
			     24*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) + 9*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 
			     6*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p2,4) + 6*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,4) - 
			     12*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,4) - 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,4) + 
			     6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,4) - 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
			     12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
			     12*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 
			     6*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p2,4) + 3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,4) + 
			     3*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c23,4)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
			     3*cpow(c13,2)*c22*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 18*c12*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
			     3*c11*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 3*cpow(c13,2)*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
			     18*c12*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 9*cpow(c12,2)*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
			     6*c11*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 9*cpow(c12,2)*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
			     3*c11*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 24*cpow(c13,2)*cpow(c23,3)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
			     6*cpow(c13,2)*c22*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 54*c12*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
			     12*c11*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
			     18*cpow(c12,2)*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
			     12*c11*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 24*cpow(c13,2)*cpow(c23,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
			     9*cpow(c13,2)*c22*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 36*c12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
			     12*c11*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 27*cpow(c12,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
			     18*c11*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 12*c13*cpow(c23,4)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
			     6*c13*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 18*c12*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
			     6*c13*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 18*c12*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
			     48*c13*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 12*c13*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
			     54*c12*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 18*c12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
			     48*c13*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 18*c13*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
			     36*c12*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
			     6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
			     6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 
			     6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) - 
			     18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
			     18*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 18*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
			     18*c12*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 18*c12*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 
			     60*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 6*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
			     54*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 48*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
			     72*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c23,4)*c33*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 
			     12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
			     42*cpow(c23,3)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 42*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
			     96*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
			     72*c23*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
			     6*c11*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
			     18*c12*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 30*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
			     3*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 3*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
			     24*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
			     36*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 6*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,3)*cpow(p1,2)*cpow(p2,4) - 
			     2*cpow(c23,6)*cpow(c55,3)*cpow(p2,6) + 6*c22*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,6) - 6*cpow(c22,2)*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,6) + 
			     2*cpow(c22,3)*cpow(c33,3)*cpow(c55,3)*cpow(p2,6) - 12*cpow(c23,5)*c44*cpow(c55,3)*cpow(p2,6) + 24*c22*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,6) - 
			     12*cpow(c22,2)*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,6) - 24*cpow(c23,4)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) + 
			     33*c22*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 9*cpow(c22,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 
			     16*cpow(c23,3)*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) + 18*c22*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) - 
			     3*cpow(c23,4)*c33*c44*cpow(c55,2)*c66*cpow(p2,6) + 6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,6) - 
			     3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) + 
			     12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 
			     18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) - 
			     3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,6) + 
			     2*cpow(c33,3)*cpow(c44,3)*cpow(c66,3)*cpow(p2,6) + csqrt(cpow(-2*cpow(c33,3)*cpow(c44,3) + 3*cpow(c33,3)*cpow(c44,2)*c55 + 
											   3*cpow(c33,2)*cpow(c44,3)*c55 + 3*cpow(c33,3)*c44*cpow(c55,2) - 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2) + 3*c33*cpow(c44,3)*cpow(c55,2) - 
											   2*cpow(c33,3)*cpow(c55,3) + 3*cpow(c33,2)*c44*cpow(c55,3) + 3*c33*cpow(c44,2)*cpow(c55,3) - 2*cpow(c44,3)*cpow(c55,3) - 
											   6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,2) + 6*c11*cpow(c33,3)*cpow(c44,3)*cpow(p1,2) + 
											   6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2) - 6*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2) - 
											   3*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,2) - 6*c11*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) - 12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) + 
											   3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2) - 3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2) + 
											   6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 
											   12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) - 6*cpow(c13,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 
											   3*c11*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 6*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 9*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) + 
											   6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2) + 12*c13*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) + 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) - 
											   12*c13*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 3*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2) - 
											   6*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2) + 6*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2) - 
											   6*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2) - 3*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2) - 6*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,4) + 
											   12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,4) - 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*cpow(p1,4) + 
											   3*cpow(c13,4)*c33*cpow(c44,2)*c55*cpow(p1,4) - 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4) + 
											   3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4) - 6*cpow(c13,4)*cpow(c44,3)*c55*cpow(p1,4) + 
											   3*c11*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,4) - 24*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4) + 
											   3*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 24*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 
											   12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 
											   24*cpow(c13,3)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 6*c11*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) - 
											   33*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
											   12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 
											   24*cpow(c13,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) + 9*c11*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 
											   6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4) + 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4) - 
											   6*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4) + 6*c11*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4) - 
											   6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 
											   12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c13*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4) - 
											   12*c13*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) - 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) + 
											   3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4) - 6*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) + 
											   3*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) - 2*cpow(c13,6)*cpow(c44,3)*cpow(p1,6) + 6*c11*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,6) - 
											   6*cpow(c11,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,6) + 2*cpow(c11,3)*cpow(c33,3)*cpow(c44,3)*cpow(p1,6) - 
											   12*cpow(c13,5)*cpow(c44,3)*c55*cpow(p1,6) + 24*c11*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,6) - 
											   12*cpow(c11,2)*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,6) - 24*cpow(c13,4)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) + 
											   33*c11*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 9*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 
											   16*cpow(c13,3)*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) + 18*c11*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) - 
											   3*cpow(c13,4)*c33*cpow(c44,2)*c55*c66*cpow(p1,6) + 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,6) - 
											   3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,6) - 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) + 
											   12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) - 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 
											   18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) - 
											   3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) + 6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,6) + 
											   2*cpow(c33,3)*cpow(c55,3)*cpow(c66,3)*cpow(p1,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p2,2) - 
											   3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p2,2) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p2,2) + 
											   6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,2) - 6*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,2) + 
											   6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 
											   12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 12*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) + 
											   18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(p2,2) - 
											   3*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,2) - 6*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 
											   6*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 3*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 6*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 
											   9*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 12*c23*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) + 
											   6*cpow(c33,3)*cpow(c44,3)*c66*cpow(p2,2) - 6*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,2) - 6*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,2) - 
											   3*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,2) - 
											   3*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,2) - 3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
											   6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
											   6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 
											   6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) + 
											   18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 
											   3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
											   18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   12*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
											   18*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   6*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
											   18*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   12*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   12*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   24*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   9*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c11*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
											   18*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
											   12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
											   24*c13*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 12*c13*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
											   18*c12*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
											   9*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
											   18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 48*c13*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
											   18*c12*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
											   18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 12*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 
											   12*c11*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
											   18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) - 
											   18*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 3*cpow(c13,2)*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
											   6*c11*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 42*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
											   18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 
											   18*c12*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
											   18*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 36*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
											   6*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 6*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
											   24*c13*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
											   12*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 12*c22*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
											   3*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 6*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
											   42*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 24*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
											   18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 36*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 
											   3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) + 
											   6*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,4)*cpow(c23,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
											   3*cpow(c13,4)*c22*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 18*c12*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
											   3*c11*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
											   9*cpow(c12,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
											   6*c11*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
											   3*cpow(c11,2)*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
											   9*c11*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 3*cpow(c11,2)*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
											   12*cpow(c13,4)*c23*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 18*c12*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 
											   6*c11*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 
											   6*cpow(c11,2)*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 24*cpow(c13,3)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
											   12*cpow(c13,3)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 54*c12*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
											   6*c11*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
											   18*cpow(c12,2)*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
											   12*c11*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
											   18*c11*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 48*cpow(c13,3)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
											   54*c12*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 12*c11*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
											   18*c11*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 24*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
											   12*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 36*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
											   9*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 27*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
											   18*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 48*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
											   36*c12*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 18*c11*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
											   6*cpow(c13,4)*c33*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) - 12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 
											   6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 18*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 
											   18*c12*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 18*c11*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
											   18*c11*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 42*cpow(c13,3)*c33*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) - 
											   42*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) + 3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
											   6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
											   6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
											   6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 60*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
											   54*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
											   6*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 96*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
											   18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
											   12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
											   48*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 72*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
											   72*c13*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 3*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
											   3*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
											   18*c12*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
											   24*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
											   6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
											   30*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 36*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
											   6*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,3)*cpow(p1,4)*cpow(p2,2) + 3*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p2,4) - 
											   6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,4) + 3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,4) + 
											   12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) - 12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) + 
											   12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 
											   6*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,4) + 12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,4) - 
											   6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*cpow(p2,4) - 6*cpow(c23,4)*c44*cpow(c55,3)*cpow(p2,4) + 3*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,4) - 
											   24*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,4) + 3*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) + 
											   24*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) - 24*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 
											   6*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 33*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 
											   18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 24*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) + 
											   9*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 
											   6*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p2,4) + 6*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,4) - 
											   12*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,4) - 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,4) + 
											   6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,4) - 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
											   12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
											   12*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 
											   6*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p2,4) + 3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,4) + 
											   3*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c23,4)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
											   3*cpow(c13,2)*c22*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 18*c12*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   3*c11*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 3*cpow(c13,2)*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   18*c12*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   9*cpow(c12,2)*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
											   6*c11*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 9*cpow(c12,2)*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   3*c11*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 24*cpow(c13,2)*cpow(c23,3)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
											   6*cpow(c13,2)*c22*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
											   54*c12*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 12*c11*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   18*c12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   18*cpow(c12,2)*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
											   12*c11*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   24*cpow(c13,2)*cpow(c23,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 9*cpow(c13,2)*c22*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
											   36*c12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 12*c11*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   27*cpow(c12,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 18*c11*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   12*c13*cpow(c23,4)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 6*c13*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
											   18*c12*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 6*c13*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
											   18*c12*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 48*c13*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
											   12*c13*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 54*c12*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
											   18*c12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 48*c13*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
											   18*c13*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 36*c12*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
											   3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 
											   18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
											   9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
											   6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
											   12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 18*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
											   18*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 18*c12*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 
											   18*c12*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 60*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
											   6*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
											   54*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 48*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
											   72*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c23,4)*c33*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 
											   12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
											   42*cpow(c23,3)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 42*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
											   96*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
											   72*c23*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
											   6*c11*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
											   18*c12*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 30*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
											   3*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
											   3*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
											   24*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
											   36*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 6*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,3)*cpow(p1,2)*cpow(p2,4) - 
											   2*cpow(c23,6)*cpow(c55,3)*cpow(p2,6) + 6*c22*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,6) - 
											   6*cpow(c22,2)*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,6) + 2*cpow(c22,3)*cpow(c33,3)*cpow(c55,3)*cpow(p2,6) - 
											   12*cpow(c23,5)*c44*cpow(c55,3)*cpow(p2,6) + 24*c22*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,6) - 
											   12*cpow(c22,2)*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,6) - 24*cpow(c23,4)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) + 
											   33*c22*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 9*cpow(c22,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 
											   16*cpow(c23,3)*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) + 18*c22*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) - 
											   3*cpow(c23,4)*c33*c44*cpow(c55,2)*c66*cpow(p2,6) + 6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,6) - 
											   3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) + 
											   12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 
											   18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) - 
											   3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,6) + 
											   2*cpow(c33,3)*cpow(c44,3)*cpow(c66,3)*cpow(p2,6),2) + 
										      4*cpow(-cpow(c33*c44 + c33*c55 + c44*c55 + cpow(c13,2)*c44*cpow(p1,2) - c11*c33*c44*cpow(p1,2) + 2*c13*c44*c55*cpow(p1,2) - c33*c55*c66*cpow(p1,2) + 
												   cpow(c23,2)*c55*cpow(p2,2) - c22*c33*c55*cpow(p2,2) + 2*c23*c44*c55*cpow(p2,2) - c33*c44*c66*cpow(p2,2),2) - 
											     3*c33*c44*c55*(-c33 - c44 - c55 - cpow(c13,2)*cpow(p1,2) + c11*c33*cpow(p1,2) + c11*c44*cpow(p1,2) - 2*c13*c55*cpow(p1,2) + c44*c55*cpow(p1,2) + 
													    c33*c66*cpow(p1,2) + c55*c66*cpow(p1,2) - c11*c44*c55*cpow(p1,4) + cpow(c13,2)*c66*cpow(p1,4) - c11*c33*c66*cpow(p1,4) + 
													    2*c13*c55*c66*cpow(p1,4) - cpow(c23,2)*cpow(p2,2) + c22*c33*cpow(p2,2) - 2*c23*c44*cpow(p2,2) + c22*c55*cpow(p2,2) + c44*c55*cpow(p2,2) + 
													    c33*c66*cpow(p2,2) + c44*c66*cpow(p2,2) + cpow(c13,2)*c22*cpow(p1,2)*cpow(p2,2) - 2*c12*c13*c23*cpow(p1,2)*cpow(p2,2) + 
													    c11*cpow(c23,2)*cpow(p1,2)*cpow(p2,2) + cpow(c12,2)*c33*cpow(p1,2)*cpow(p2,2) - c11*c22*c33*cpow(p1,2)*cpow(p2,2) - 
													    2*c12*c13*c44*cpow(p1,2)*cpow(p2,2) + 2*c11*c23*c44*cpow(p1,2)*cpow(p2,2) + 2*c13*c22*c55*cpow(p1,2)*cpow(p2,2) - 
													    2*c12*c23*c55*cpow(p1,2)*cpow(p2,2) - 2*c12*c44*c55*cpow(p1,2)*cpow(p2,2) - 2*c13*c23*c66*cpow(p1,2)*cpow(p2,2) + 
													    2*c12*c33*c66*cpow(p1,2)*cpow(p2,2) - 2*c13*c44*c66*cpow(p1,2)*cpow(p2,2) - 2*c23*c55*c66*cpow(p1,2)*cpow(p2,2) - 
													    4*c44*c55*c66*cpow(p1,2)*cpow(p2,2) - c22*c44*c55*cpow(p2,4) + cpow(c23,2)*c66*cpow(p2,4) - c22*c33*c66*cpow(p2,4) + 2*c23*c44*c66*cpow(p2,4)),3)),
			     0.3333333333333333)) - cpow(-2*cpow(c33,3)*cpow(c44,3) + 3*cpow(c33,3)*cpow(c44,2)*c55 + 3*cpow(c33,2)*cpow(c44,3)*c55 + 3*cpow(c33,3)*c44*cpow(c55,2) - 
							 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2) + 3*c33*cpow(c44,3)*cpow(c55,2) - 2*cpow(c33,3)*cpow(c55,3) + 3*cpow(c33,2)*c44*cpow(c55,3) + 
							 3*c33*cpow(c44,2)*cpow(c55,3) - 2*cpow(c44,3)*cpow(c55,3) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,2) + 
							 6*c11*cpow(c33,3)*cpow(c44,3)*cpow(p1,2) + 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2) - 6*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2) - 
							 3*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,2) - 6*c11*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) - 12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) + 
							 3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2) - 3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2) + 
							 6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 
							 12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) - 6*cpow(c13,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 
							 3*c11*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 6*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 9*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) + 
							 6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2) + 12*c13*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) + 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) - 
							 12*c13*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 3*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2) - 
							 6*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2) + 6*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2) - 
							 6*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2) - 3*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2) - 6*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,4) + 
							 12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,4) - 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*cpow(p1,4) + 
							 3*cpow(c13,4)*c33*cpow(c44,2)*c55*cpow(p1,4) - 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4) + 
							 3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4) - 6*cpow(c13,4)*cpow(c44,3)*c55*cpow(p1,4) + 3*c11*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,4) - 
							 24*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4) + 3*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 
							 24*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 
							 12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 24*cpow(c13,3)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
							 6*c11*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) - 33*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
							 18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 
							 18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 24*cpow(c13,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) + 
							 9*c11*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4) + 
							 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4) - 6*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4) + 
							 6*c11*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4) - 6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 
							 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 
							 12*c13*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4) - 12*c13*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) - 
							 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) + 3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4) - 
							 6*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) + 3*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) - 2*cpow(c13,6)*cpow(c44,3)*cpow(p1,6) + 
							 6*c11*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,6) - 6*cpow(c11,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,6) + 
							 2*cpow(c11,3)*cpow(c33,3)*cpow(c44,3)*cpow(p1,6) - 12*cpow(c13,5)*cpow(c44,3)*c55*cpow(p1,6) + 24*c11*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,6) - 
							 12*cpow(c11,2)*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,6) - 24*cpow(c13,4)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) + 
							 33*c11*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 9*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 
							 16*cpow(c13,3)*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) + 18*c11*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) - 
							 3*cpow(c13,4)*c33*cpow(c44,2)*c55*c66*cpow(p1,6) + 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,6) - 
							 3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,6) - 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) + 
							 12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) - 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 
							 18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) - 
							 3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) + 6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,6) + 
							 2*cpow(c33,3)*cpow(c55,3)*cpow(c66,3)*cpow(p1,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p2,2) - 3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p2,2) + 
							 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p2,2) + 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,2) - 6*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,2) + 
							 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 
							 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 12*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) + 
							 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(p2,2) - 
							 3*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,2) - 6*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 
							 6*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 3*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 6*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 
							 9*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 12*c23*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) + 
							 6*cpow(c33,3)*cpow(c44,3)*c66*cpow(p2,2) - 6*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,2) - 6*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,2) - 
							 3*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,2) - 3*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,2) - 
							 3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 
							 18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
							 9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 
							 12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
							 6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
							 18*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
							 18*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 24*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
							 18*c12*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 9*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 12*c11*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
							 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 24*c13*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
							 12*c13*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
							 12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 9*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
							 18*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
							 48*c13*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
							 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
							 12*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 12*c11*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
							 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) - 18*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
							 3*cpow(c13,2)*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 6*c11*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
							 42*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 18*c12*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
							 6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 18*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
							 36*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 6*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
							 6*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 24*c13*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
							 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 12*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 
							 12*c22*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 3*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
							 6*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 42*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
							 24*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
							 36*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 
							 3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) + 6*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c13,4)*cpow(c23,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 3*cpow(c13,4)*c22*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
							 18*c12*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 3*c11*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
							 9*cpow(c12,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 6*c11*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
							 18*c11*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 3*cpow(c11,2)*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
							 9*c11*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 3*cpow(c11,2)*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
							 12*cpow(c13,4)*c23*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 18*c12*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 
							 6*c11*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 
							 6*cpow(c11,2)*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 24*cpow(c13,3)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
							 12*cpow(c13,3)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 54*c12*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
							 6*c11*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*cpow(c12,2)*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
							 12*c11*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
							 48*cpow(c13,3)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 54*c12*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
							 12*c11*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
							 24*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
							 36*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 9*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
							 27*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 18*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
							 48*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 36*c12*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
							 18*c11*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 6*cpow(c13,4)*c33*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) - 
							 12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 
							 18*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 
							 18*c11*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 18*c11*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
							 42*cpow(c13,3)*c33*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) - 42*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
							 3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
							 18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
							 9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
							 60*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 54*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
							 6*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 96*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
							 18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
							 12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
							 48*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 72*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
							 72*c13*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 3*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
							 3*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
							 18*c12*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 24*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
							 6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
							 30*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 36*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
							 6*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,3)*cpow(p1,4)*cpow(p2,2) + 3*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p2,4) - 
							 6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,4) + 3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,4) + 
							 12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) - 12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) + 
							 12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 
							 6*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,4) + 12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,4) - 6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*cpow(p2,4) - 
							 6*cpow(c23,4)*c44*cpow(c55,3)*cpow(p2,4) + 3*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,4) - 24*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,4) + 
							 3*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) + 24*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) - 
							 24*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 6*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 
							 33*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 
							 24*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) + 9*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 
							 6*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p2,4) + 6*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,4) - 
							 12*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,4) - 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,4) + 
							 6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,4) - 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
							 12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
							 12*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 
							 6*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p2,4) + 3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,4) + 
							 3*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c23,4)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
							 3*cpow(c13,2)*c22*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 18*c12*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
							 3*c11*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 3*cpow(c13,2)*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
							 18*c12*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 9*cpow(c12,2)*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
							 6*c11*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 9*cpow(c12,2)*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
							 3*c11*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 24*cpow(c13,2)*cpow(c23,3)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
							 6*cpow(c13,2)*c22*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 54*c12*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
							 12*c11*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
							 18*cpow(c12,2)*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
							 12*c11*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 24*cpow(c13,2)*cpow(c23,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
							 9*cpow(c13,2)*c22*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 36*c12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
							 12*c11*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 27*cpow(c12,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
							 18*c11*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 12*c13*cpow(c23,4)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
							 6*c13*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 18*c12*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
							 6*c13*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 18*c12*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
							 48*c13*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 12*c13*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
							 54*c12*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 18*c12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
							 48*c13*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 18*c13*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
							 36*c12*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
							 6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
							 6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 
							 6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) - 
							 18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
							 18*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 18*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
							 18*c12*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 18*c12*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 
							 60*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 6*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
							 54*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 48*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
							 72*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c23,4)*c33*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 
							 12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
							 42*cpow(c23,3)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 42*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
							 96*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
							 72*c23*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
							 6*c11*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
							 18*c12*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 30*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
							 3*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 3*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
							 24*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
							 36*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 6*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,3)*cpow(p1,2)*cpow(p2,4) - 
							 2*cpow(c23,6)*cpow(c55,3)*cpow(p2,6) + 6*c22*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,6) - 6*cpow(c22,2)*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,6) + 
							 2*cpow(c22,3)*cpow(c33,3)*cpow(c55,3)*cpow(p2,6) - 12*cpow(c23,5)*c44*cpow(c55,3)*cpow(p2,6) + 24*c22*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,6) - 
							 12*cpow(c22,2)*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,6) - 24*cpow(c23,4)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) + 
							 33*c22*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 9*cpow(c22,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 
							 16*cpow(c23,3)*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) + 18*c22*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) - 
							 3*cpow(c23,4)*c33*c44*cpow(c55,2)*c66*cpow(p2,6) + 6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,6) - 
							 3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) + 
							 12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 
							 18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) - 
							 3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,6) + 
							 2*cpow(c33,3)*cpow(c44,3)*cpow(c66,3)*cpow(p2,6) + csqrt(cpow(-2*cpow(c33,3)*cpow(c44,3) + 3*cpow(c33,3)*cpow(c44,2)*c55 + 
														       3*cpow(c33,2)*cpow(c44,3)*c55 + 3*cpow(c33,3)*c44*cpow(c55,2) - 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2) + 3*c33*cpow(c44,3)*cpow(c55,2) - 
														       2*cpow(c33,3)*cpow(c55,3) + 3*cpow(c33,2)*c44*cpow(c55,3) + 3*c33*cpow(c44,2)*cpow(c55,3) - 2*cpow(c44,3)*cpow(c55,3) - 
														       6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,2) + 6*c11*cpow(c33,3)*cpow(c44,3)*cpow(p1,2) + 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2) - 
														       6*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2) - 3*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,2) - 6*c11*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) - 
														       12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) + 3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2) - 3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2) + 
														       6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 
														       12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) - 6*cpow(c13,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 
														       3*c11*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 6*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 9*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) + 
														       6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2) + 12*c13*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) + 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) - 
														       12*c13*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 3*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2) - 
														       6*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2) + 6*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2) - 
														       6*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2) - 3*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2) - 6*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,4) + 
														       12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,4) - 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*cpow(p1,4) + 
														       3*cpow(c13,4)*c33*cpow(c44,2)*c55*cpow(p1,4) - 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4) + 
														       3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4) - 6*cpow(c13,4)*cpow(c44,3)*c55*cpow(p1,4) + 3*c11*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,4) - 
														       24*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4) + 3*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 
														       24*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 
														       12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 24*cpow(c13,3)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
														       6*c11*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) - 33*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
														       18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 
														       18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 24*cpow(c13,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) + 
														       9*c11*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 
														       6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4) + 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4) - 
														       6*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4) + 6*c11*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4) - 
														       6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 
														       12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c13*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4) - 
														       12*c13*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) - 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) + 
														       3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4) - 6*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) + 
														       3*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) - 2*cpow(c13,6)*cpow(c44,3)*cpow(p1,6) + 6*c11*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,6) - 
														       6*cpow(c11,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,6) + 2*cpow(c11,3)*cpow(c33,3)*cpow(c44,3)*cpow(p1,6) - 
														       12*cpow(c13,5)*cpow(c44,3)*c55*cpow(p1,6) + 24*c11*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,6) - 
														       12*cpow(c11,2)*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,6) - 24*cpow(c13,4)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) + 
														       33*c11*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 9*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 
														       16*cpow(c13,3)*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) + 18*c11*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) - 
														       3*cpow(c13,4)*c33*cpow(c44,2)*c55*c66*cpow(p1,6) + 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,6) - 
														       3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,6) - 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) + 
														       12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) - 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 
														       18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) - 
														       3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) + 6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,6) + 
														       2*cpow(c33,3)*cpow(c55,3)*cpow(c66,3)*cpow(p1,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p2,2) - 
														       3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p2,2) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p2,2) + 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,2) - 
														       6*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,2) + 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 
														       12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 
														       12*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) + 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) - 
														       6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(p2,2) - 3*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,2) - 
														       6*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 
														       3*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 6*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 9*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 
														       12*c23*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) + 6*cpow(c33,3)*cpow(c44,3)*c66*cpow(p2,2) - 
														       6*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,2) - 6*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,2) - 3*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,2) + 
														       12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,2) - 3*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,2) - 
														       3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 
														       18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
														       9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 
														       12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
														       6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
														       18*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
														       18*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 24*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
														       18*c12*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 9*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       12*c11*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
														       18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 24*c13*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
														       12*c13*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
														       12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 9*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
														       18*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
														       48*c13*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
														       18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
														       12*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 12*c11*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
														       6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) - 18*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
														       3*cpow(c13,2)*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 6*c11*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
														       42*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 18*c12*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
														       6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 18*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
														       36*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 6*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
														       6*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 24*c13*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
														       18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 12*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 
														       12*c22*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 3*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
														       6*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 42*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
														       24*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
														       36*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 
														       3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) + 6*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c13,4)*cpow(c23,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 3*cpow(c13,4)*c22*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
														       18*c12*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 3*c11*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
														       9*cpow(c12,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
														       6*c11*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
														       3*cpow(c11,2)*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 9*c11*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
														       3*cpow(c11,2)*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,4)*c23*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 
														       18*c12*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 6*c11*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 
														       18*c11*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 6*cpow(c11,2)*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 
														       24*cpow(c13,3)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,3)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
														       54*c12*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 6*c11*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
														       18*cpow(c12,2)*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
														       12*c11*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
														       48*cpow(c13,3)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 54*c12*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
														       12*c11*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
														       24*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
														       36*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 9*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
														       27*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 18*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
														       48*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 36*c12*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
														       18*c11*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 6*cpow(c13,4)*c33*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) - 
														       12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 
														       18*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 
														       18*c11*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 18*c11*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
														       42*cpow(c13,3)*c33*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) - 42*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
														       3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
														       18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
														       9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
														       60*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 54*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
														       6*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 96*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
														       18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
														       12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
														       48*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 72*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
														       72*c13*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 3*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
														       3*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
														       18*c12*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
														       24*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
														       6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
														       30*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 36*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
														       6*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,3)*cpow(p1,4)*cpow(p2,2) + 3*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p2,4) - 
														       6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,4) + 3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,4) + 
														       12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) - 12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) + 
														       12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 
														       6*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,4) + 12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,4) - 
														       6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*cpow(p2,4) - 6*cpow(c23,4)*c44*cpow(c55,3)*cpow(p2,4) + 3*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,4) - 
														       24*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,4) + 3*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) + 
														       24*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) - 24*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 
														       6*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 33*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 
														       18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 24*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) + 
														       9*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 
														       6*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p2,4) + 6*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,4) - 
														       12*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,4) - 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,4) + 
														       6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,4) - 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
														       12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
														       12*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 
														       6*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p2,4) + 3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,4) + 
														       3*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c23,4)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
														       3*cpow(c13,2)*c22*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 18*c12*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       3*c11*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 3*cpow(c13,2)*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       18*c12*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 9*cpow(c12,2)*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
														       6*c11*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 9*cpow(c12,2)*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       3*c11*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 24*cpow(c13,2)*cpow(c23,3)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
														       6*cpow(c13,2)*c22*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 54*c12*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       12*c11*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       18*cpow(c12,2)*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
														       12*c11*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       24*cpow(c13,2)*cpow(c23,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 9*cpow(c13,2)*c22*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
														       36*c12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 12*c11*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       27*cpow(c12,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 18*c11*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       12*c13*cpow(c23,4)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 6*c13*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
														       18*c12*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 6*c13*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
														       18*c12*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 48*c13*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
														       12*c13*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 54*c12*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
														       18*c12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 48*c13*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
														       18*c13*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 36*c12*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
														       3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 
														       18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
														       9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
														       6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
														       12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 18*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
														       18*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 18*c12*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 
														       18*c12*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 60*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
														       6*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 54*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 
														       48*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 72*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 
														       6*cpow(c23,4)*c33*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
														       6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 42*cpow(c23,3)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 
														       42*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 96*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 
														       18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 72*c23*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 
														       6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 6*c11*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
														       18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 18*c12*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
														       30*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 3*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
														       3*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 24*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
														       36*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 6*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,3)*cpow(p1,2)*cpow(p2,4) - 
														       2*cpow(c23,6)*cpow(c55,3)*cpow(p2,6) + 6*c22*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,6) - 6*cpow(c22,2)*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,6) + 
														       2*cpow(c22,3)*cpow(c33,3)*cpow(c55,3)*cpow(p2,6) - 12*cpow(c23,5)*c44*cpow(c55,3)*cpow(p2,6) + 24*c22*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,6) - 
														       12*cpow(c22,2)*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,6) - 24*cpow(c23,4)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) + 
														       33*c22*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 9*cpow(c22,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 
														       16*cpow(c23,3)*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) + 18*c22*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) - 
														       3*cpow(c23,4)*c33*c44*cpow(c55,2)*c66*cpow(p2,6) + 6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,6) - 
														       3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) + 
														       12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 
														       18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) - 
														       3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,6) + 
														       2*cpow(c33,3)*cpow(c44,3)*cpow(c66,3)*cpow(p2,6),2) + 4*
														  cpow(-cpow(c33*c44 + c33*c55 + c44*c55 + cpow(c13,2)*c44*cpow(p1,2) - c11*c33*c44*cpow(p1,2) + 2*c13*c44*c55*cpow(p1,2) - c33*c55*c66*cpow(p1,2) + 
															     cpow(c23,2)*c55*cpow(p2,2) - c22*c33*c55*cpow(p2,2) + 2*c23*c44*c55*cpow(p2,2) - c33*c44*c66*cpow(p2,2),2) - 
														       3*c33*c44*c55*(-c33 - c44 - c55 - cpow(c13,2)*cpow(p1,2) + c11*c33*cpow(p1,2) + c11*c44*cpow(p1,2) - 2*c13*c55*cpow(p1,2) + c44*c55*cpow(p1,2) + 
																      c33*c66*cpow(p1,2) + c55*c66*cpow(p1,2) - c11*c44*c55*cpow(p1,4) + cpow(c13,2)*c66*cpow(p1,4) - c11*c33*c66*cpow(p1,4) + 2*c13*c55*c66*cpow(p1,4) - 
																      cpow(c23,2)*cpow(p2,2) + c22*c33*cpow(p2,2) - 2*c23*c44*cpow(p2,2) + c22*c55*cpow(p2,2) + c44*c55*cpow(p2,2) + c33*c66*cpow(p2,2) + 
																      c44*c66*cpow(p2,2) + cpow(c13,2)*c22*cpow(p1,2)*cpow(p2,2) - 2*c12*c13*c23*cpow(p1,2)*cpow(p2,2) + c11*cpow(c23,2)*cpow(p1,2)*cpow(p2,2) + 
																      cpow(c12,2)*c33*cpow(p1,2)*cpow(p2,2) - c11*c22*c33*cpow(p1,2)*cpow(p2,2) - 2*c12*c13*c44*cpow(p1,2)*cpow(p2,2) + 
																      2*c11*c23*c44*cpow(p1,2)*cpow(p2,2) + 2*c13*c22*c55*cpow(p1,2)*cpow(p2,2) - 2*c12*c23*c55*cpow(p1,2)*cpow(p2,2) - 
																      2*c12*c44*c55*cpow(p1,2)*cpow(p2,2) - 2*c13*c23*c66*cpow(p1,2)*cpow(p2,2) + 2*c12*c33*c66*cpow(p1,2)*cpow(p2,2) - 
																      2*c13*c44*c66*cpow(p1,2)*cpow(p2,2) - 2*c23*c55*c66*cpow(p1,2)*cpow(p2,2) - 4*c44*c55*c66*cpow(p1,2)*cpow(p2,2) - c22*c44*c55*cpow(p2,4) + 
																      cpow(c23,2)*c66*cpow(p2,4) - c22*c33*c66*cpow(p2,4) + 2*c23*c44*c66*cpow(p2,4)),3)),0.3333333333333333)/(3.*cpow(2,0.3333333333333333)*c33*c44*c55);
              

/*dx ###################################################################################################################*/

    dx = p1*(2*(c23 + c44)*(-(c11*(c23 + c44)) + c12*(c13 + c55) + (c13 + c55)*c66)*cpow(p2,2)*cpow(p3,2) + (c13 + c55)*cpow(p3,2)*((c23 + c44)*(c12 + c66)*cpow(p2,2) - (c13 + c55)*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))) + (c13 + c55)*cpow(p3,2)*(-2*(c13 + c55)*c66*cpow(p1,2) + (c23 + c44)*(c12 + c66)*cpow(p2,2) - (c13 + c55)*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))) + 2*(-1 + c55*cpow(p1,2) + c44*cpow(p2,2) + c33*cpow(p3,2))*(-(cpow(c12,2)*cpow(p2,2)) + c11*(-1 + 2*c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2)) + c66*(-1 - 2*c12*cpow(p2,2) + c55*cpow(p3,2))) + 2*c55*(-(cpow(c12 + c66,2)*cpow(p1,2)*cpow(p2,2)) + (-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))));

/*dt ###################################################################################################################*/

    dt = cpow(p1,2)*(2*(c23 + c44)*(-(c11*(c23 + c44)) + c12*(c13 + c55) + (c13 + c55)*c66)*cpow(p2,2)*cpow(p3,2) + 
		     (c13 + c55)*cpow(p3,2)*((c23 + c44)*(c12 + c66)*cpow(p2,2) - (c13 + c55)*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))) + 
		     (c13 + c55)*cpow(p3,2)*(-2*(c13 + c55)*c66*cpow(p1,2) + (c23 + c44)*(c12 + c66)*cpow(p2,2) - 
					     (c13 + c55)*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))) + 
		     2*(-1 + c55*cpow(p1,2) + c44*cpow(p2,2) + c33*cpow(p3,2))*(-(cpow(c12,2)*cpow(p2,2)) + c11*(-1 + 2*c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2)) + 
										c66*(-1 - 2*c12*cpow(p2,2) + c55*cpow(p3,2))) + 2*c55*(-(cpow(c12 + c66,2)*cpow(p1,2)*cpow(p2,2)) + 
																       (-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2)))) + 
	cpow(p2,2)*(-2*(c13 + c55)*(c13*c22 - c12*(c23 + c44) + c22*c55 - c23*c66 - c44*c66)*cpow(p1,2)*cpow(p3,2) - 
		    (c23 + c44)*cpow(p3,2)*(-((c13 + c55)*(c12 + c66)*cpow(p1,2)) + (c23 + c44)*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) - 
		    (c23 + c44)*cpow(p3,2)*(-((c13 + c55)*(c12 + c66)*cpow(p1,2)) + 2*(c23 + c44)*c66*cpow(p2,2) + 
					    (c23 + c44)*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) + 
		    2*c44*(-(cpow(c12 + c66,2)*cpow(p1,2)*cpow(p2,2)) + (-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))*
			   (-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) + 
		    2*(-1 + c55*cpow(p1,2) + c44*cpow(p2,2) + c33*cpow(p3,2))*(-(cpow(c12,2)*cpow(p1,2)) + c66*(-1 - 2*c12*cpow(p1,2) + c44*cpow(p3,2)) + 
									       c22*(-1 + c11*cpow(p1,2) + 2*c66*cpow(p2,2) + c55*cpow(p3,2)))) + 
	cpow(p3,2)*((c13 + c55)*cpow(p1,2)*((c23 + c44)*(c12 + c66)*cpow(p2,2) - (c13 + c55)*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))) + 
		    (c13 + c55)*cpow(p1,2)*((c23 + c44)*(c12 + c66)*cpow(p2,2) - 2*c44*(c13 + c55)*cpow(p3,2) - 
					    (c13 + c55)*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))) - 
		    (c23 + c44)*cpow(p2,2)*(-((c13 + c55)*(c12 + c66)*cpow(p1,2)) + (c23 + c44)*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) - 
		    (c23 + c44)*cpow(p2,2)*(-((c13 + c55)*(c12 + c66)*cpow(p1,2)) + 2*(c23 + c44)*c55*cpow(p3,2) + 
					    (c23 + c44)*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) + 
		    2*c33*(-(cpow(c12 + c66,2)*cpow(p1,2)*cpow(p2,2)) + (-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))*
			   (-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) + 
		    2*(-1 + c55*cpow(p1,2) + c44*cpow(p2,2) + c33*cpow(p3,2))*(c55*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2)) + 
									       c44*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + 2*c55*cpow(p3,2))));


/*sf_warning("p1 %f p2 %f p3 %f dx %f dt %f",p1,p2,creal(p3),creal(dx),creal(dt));*/

    return creal(dx)/creal(dt);

}


float yPvelo(float c11, float c22,float c33,float c44,float c55,float c66,float c12,float c13,float c23, float p1, float p2) {
/*<Compute qP-wave group velocity (y) given cij and px (p1) and py (p2) for layered orthorhombic media>*/

    float complex p3, dy, dt;


/*p3 ###################################################################################################################*/


    p3 = (c33*c44 + c33*c55 + c44*c55 + cpow(c13,2)*c44*cpow(p1,2) - c11*c33*c44*cpow(p1,2) + 2*c13*c44*c55*cpow(p1,2) - c33*c55*c66*cpow(p1,2) + 
	  cpow(c23,2)*c55*cpow(p2,2) - c22*c33*c55*cpow(p2,2) + 2*c23*c44*c55*cpow(p2,2) - c33*c44*c66*cpow(p2,2))/(3.*c33*c44*c55) + 
	(cpow(2,0.3333333333333333)*(-cpow(c33*c44 + c33*c55 + c44*c55 + cpow(c13,2)*c44*cpow(p1,2) - c11*c33*c44*cpow(p1,2) + 2*c13*c44*c55*cpow(p1,2) - 
					   c33*c55*c66*cpow(p1,2) + cpow(c23,2)*c55*cpow(p2,2) - c22*c33*c55*cpow(p2,2) + 2*c23*c44*c55*cpow(p2,2) - c33*c44*c66*cpow(p2,2),2) - 
				     3*c33*c44*c55*(-c33 - c44 - c55 - cpow(c13,2)*cpow(p1,2) + c11*c33*cpow(p1,2) + c11*c44*cpow(p1,2) - 2*c13*c55*cpow(p1,2) + c44*c55*cpow(p1,2) + 
						    c33*c66*cpow(p1,2) + c55*c66*cpow(p1,2) - c11*c44*c55*cpow(p1,4) + cpow(c13,2)*c66*cpow(p1,4) - c11*c33*c66*cpow(p1,4) + 2*c13*c55*c66*cpow(p1,4) - 
						    cpow(c23,2)*cpow(p2,2) + c22*c33*cpow(p2,2) - 2*c23*c44*cpow(p2,2) + c22*c55*cpow(p2,2) + c44*c55*cpow(p2,2) + c33*c66*cpow(p2,2) + 
						    c44*c66*cpow(p2,2) + cpow(c13,2)*c22*cpow(p1,2)*cpow(p2,2) - 2*c12*c13*c23*cpow(p1,2)*cpow(p2,2) + c11*cpow(c23,2)*cpow(p1,2)*cpow(p2,2) + 
						    cpow(c12,2)*c33*cpow(p1,2)*cpow(p2,2) - c11*c22*c33*cpow(p1,2)*cpow(p2,2) - 2*c12*c13*c44*cpow(p1,2)*cpow(p2,2) + 
						    2*c11*c23*c44*cpow(p1,2)*cpow(p2,2) + 2*c13*c22*c55*cpow(p1,2)*cpow(p2,2) - 2*c12*c23*c55*cpow(p1,2)*cpow(p2,2) - 
						    2*c12*c44*c55*cpow(p1,2)*cpow(p2,2) - 2*c13*c23*c66*cpow(p1,2)*cpow(p2,2) + 2*c12*c33*c66*cpow(p1,2)*cpow(p2,2) - 
						    2*c13*c44*c66*cpow(p1,2)*cpow(p2,2) - 2*c23*c55*c66*cpow(p1,2)*cpow(p2,2) - 4*c44*c55*c66*cpow(p1,2)*cpow(p2,2) - c22*c44*c55*cpow(p2,4) + 
						    cpow(c23,2)*c66*cpow(p2,4) - c22*c33*c66*cpow(p2,4) + 2*c23*c44*c66*cpow(p2,4))))/
	(3.*c33*c44*c55*cpow(-2*cpow(c33,3)*cpow(c44,3) + 3*cpow(c33,3)*cpow(c44,2)*c55 + 3*cpow(c33,2)*cpow(c44,3)*c55 + 3*cpow(c33,3)*c44*cpow(c55,2) - 
			     12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2) + 3*c33*cpow(c44,3)*cpow(c55,2) - 2*cpow(c33,3)*cpow(c55,3) + 3*cpow(c33,2)*c44*cpow(c55,3) + 
			     3*c33*cpow(c44,2)*cpow(c55,3) - 2*cpow(c44,3)*cpow(c55,3) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,2) + 
			     6*c11*cpow(c33,3)*cpow(c44,3)*cpow(p1,2) + 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2) - 6*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2) - 
			     3*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,2) - 6*c11*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) - 12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) + 
			     3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2) - 3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2) + 
			     6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 
			     12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) - 6*cpow(c13,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 
			     3*c11*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 6*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 9*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) + 
			     6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2) + 12*c13*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) + 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) - 
			     12*c13*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 3*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2) - 
			     6*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2) + 6*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2) - 
			     6*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2) - 3*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2) - 6*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,4) + 
			     12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,4) - 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*cpow(p1,4) + 
			     3*cpow(c13,4)*c33*cpow(c44,2)*c55*cpow(p1,4) - 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4) + 
			     3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4) - 6*cpow(c13,4)*cpow(c44,3)*c55*cpow(p1,4) + 3*c11*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,4) - 
			     24*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4) + 3*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 
			     24*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 
			     12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 24*cpow(c13,3)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
			     6*c11*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) - 33*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
			     18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 
			     18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 24*cpow(c13,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) + 
			     9*c11*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 
			     6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4) + 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4) - 
			     6*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4) + 6*c11*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4) - 
			     6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 
			     12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c13*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4) - 
			     12*c13*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) - 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) + 
			     3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4) - 6*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) + 
			     3*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) - 2*cpow(c13,6)*cpow(c44,3)*cpow(p1,6) + 6*c11*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,6) - 
			     6*cpow(c11,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,6) + 2*cpow(c11,3)*cpow(c33,3)*cpow(c44,3)*cpow(p1,6) - 
			     12*cpow(c13,5)*cpow(c44,3)*c55*cpow(p1,6) + 24*c11*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,6) - 
			     12*cpow(c11,2)*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,6) - 24*cpow(c13,4)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) + 
			     33*c11*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 9*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 
			     16*cpow(c13,3)*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) + 18*c11*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) - 
			     3*cpow(c13,4)*c33*cpow(c44,2)*c55*c66*cpow(p1,6) + 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,6) - 
			     3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,6) - 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) + 
			     12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) - 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 
			     18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) - 
			     3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) + 6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,6) + 
			     2*cpow(c33,3)*cpow(c55,3)*cpow(c66,3)*cpow(p1,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p2,2) - 
			     3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p2,2) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p2,2) + 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,2) - 
			     6*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,2) + 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 
			     12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 
			     12*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) + 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,2) + 
			     6*c22*cpow(c33,3)*cpow(c55,3)*cpow(p2,2) - 3*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,2) - 6*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 
			     12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 3*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 
			     6*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 9*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 12*c23*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) - 
			     9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) + 6*cpow(c33,3)*cpow(c44,3)*c66*cpow(p2,2) - 6*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,2) - 
			     6*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,2) - 3*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,2) - 
			     3*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,2) - 3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
			     6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
			     6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 
			     6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) + 
			     18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 
			     3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
			     18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     12*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
			     18*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     6*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
			     18*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     12*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     12*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     24*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     9*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c11*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
			     18*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
			     12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
			     24*c13*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 12*c13*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
			     18*c12*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
			     9*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
			     18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 48*c13*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
			     18*c12*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
			     18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 12*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 
			     12*c11*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
			     18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) - 
			     18*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 3*cpow(c13,2)*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
			     6*c11*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 42*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
			     18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 
			     18*c12*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
			     18*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 36*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
			     6*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 6*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
			     24*c13*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
			     12*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 12*c22*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
			     3*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 6*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
			     42*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 24*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
			     18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 36*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 
			     3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) + 
			     6*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,4)*cpow(c23,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
			     3*cpow(c13,4)*c22*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 18*c12*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
			     3*c11*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
			     9*cpow(c12,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
			     6*c11*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
			     3*cpow(c11,2)*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 9*c11*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
			     3*cpow(c11,2)*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,4)*c23*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 
			     18*c12*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 6*c11*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 
			     18*c11*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 6*cpow(c11,2)*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 
			     24*cpow(c13,3)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,3)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
			     54*c12*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 6*c11*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
			     18*cpow(c12,2)*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
			     12*c11*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
			     48*cpow(c13,3)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 54*c12*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
			     12*c11*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
			     24*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
			     36*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 9*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
			     27*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 18*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
			     48*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 36*c12*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
			     18*c11*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 6*cpow(c13,4)*c33*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) - 
			     12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 
			     18*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 
			     18*c11*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 18*c11*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
			     42*cpow(c13,3)*c33*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) - 42*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
			     3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
			     18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
			     9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
			     60*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 54*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
			     6*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 96*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
			     18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
			     12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
			     48*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 72*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
			     72*c13*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 3*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
			     3*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
			     18*c12*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 24*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
			     6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
			     30*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 36*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
			     6*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,3)*cpow(p1,4)*cpow(p2,2) + 3*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p2,4) - 
			     6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,4) + 3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,4) + 
			     12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) - 12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) + 
			     12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 
			     6*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,4) + 12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,4) - 6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*cpow(p2,4) - 
			     6*cpow(c23,4)*c44*cpow(c55,3)*cpow(p2,4) + 3*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,4) - 24*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,4) + 
			     3*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) + 24*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) - 
			     24*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 6*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 
			     33*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 
			     24*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) + 9*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 
			     6*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p2,4) + 6*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,4) - 
			     12*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,4) - 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,4) + 
			     6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,4) - 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
			     12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
			     12*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 
			     6*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p2,4) + 3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,4) + 
			     3*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c23,4)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
			     3*cpow(c13,2)*c22*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 18*c12*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
			     3*c11*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 3*cpow(c13,2)*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
			     18*c12*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 9*cpow(c12,2)*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
			     6*c11*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 9*cpow(c12,2)*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
			     3*c11*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 24*cpow(c13,2)*cpow(c23,3)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
			     6*cpow(c13,2)*c22*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 54*c12*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
			     12*c11*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
			     18*cpow(c12,2)*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
			     12*c11*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 24*cpow(c13,2)*cpow(c23,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
			     9*cpow(c13,2)*c22*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 36*c12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
			     12*c11*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 27*cpow(c12,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
			     18*c11*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 12*c13*cpow(c23,4)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
			     6*c13*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 18*c12*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
			     6*c13*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 18*c12*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
			     48*c13*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 12*c13*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
			     54*c12*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 18*c12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
			     48*c13*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 18*c13*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
			     36*c12*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
			     6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
			     6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 
			     6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) - 
			     18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
			     18*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 18*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
			     18*c12*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 18*c12*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 
			     60*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 6*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
			     54*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 48*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
			     72*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c23,4)*c33*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 
			     12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
			     42*cpow(c23,3)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 42*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
			     96*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
			     72*c23*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
			     6*c11*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
			     18*c12*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 30*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
			     3*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 3*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
			     24*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
			     36*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 6*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,3)*cpow(p1,2)*cpow(p2,4) - 
			     2*cpow(c23,6)*cpow(c55,3)*cpow(p2,6) + 6*c22*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,6) - 6*cpow(c22,2)*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,6) + 
			     2*cpow(c22,3)*cpow(c33,3)*cpow(c55,3)*cpow(p2,6) - 12*cpow(c23,5)*c44*cpow(c55,3)*cpow(p2,6) + 24*c22*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,6) - 
			     12*cpow(c22,2)*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,6) - 24*cpow(c23,4)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) + 
			     33*c22*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 9*cpow(c22,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 
			     16*cpow(c23,3)*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) + 18*c22*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) - 
			     3*cpow(c23,4)*c33*c44*cpow(c55,2)*c66*cpow(p2,6) + 6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,6) - 
			     3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) + 
			     12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 
			     18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) - 
			     3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,6) + 
			     2*cpow(c33,3)*cpow(c44,3)*cpow(c66,3)*cpow(p2,6) + csqrt(cpow(-2*cpow(c33,3)*cpow(c44,3) + 3*cpow(c33,3)*cpow(c44,2)*c55 + 
											   3*cpow(c33,2)*cpow(c44,3)*c55 + 3*cpow(c33,3)*c44*cpow(c55,2) - 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2) + 3*c33*cpow(c44,3)*cpow(c55,2) - 
											   2*cpow(c33,3)*cpow(c55,3) + 3*cpow(c33,2)*c44*cpow(c55,3) + 3*c33*cpow(c44,2)*cpow(c55,3) - 2*cpow(c44,3)*cpow(c55,3) - 
											   6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,2) + 6*c11*cpow(c33,3)*cpow(c44,3)*cpow(p1,2) + 
											   6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2) - 6*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2) - 
											   3*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,2) - 6*c11*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) - 12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) + 
											   3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2) - 3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2) + 
											   6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 
											   12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) - 6*cpow(c13,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 
											   3*c11*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 6*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 9*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) + 
											   6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2) + 12*c13*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) + 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) - 
											   12*c13*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 3*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2) - 
											   6*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2) + 6*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2) - 
											   6*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2) - 3*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2) - 6*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,4) + 
											   12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,4) - 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*cpow(p1,4) + 
											   3*cpow(c13,4)*c33*cpow(c44,2)*c55*cpow(p1,4) - 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4) + 
											   3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4) - 6*cpow(c13,4)*cpow(c44,3)*c55*cpow(p1,4) + 
											   3*c11*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,4) - 24*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4) + 
											   3*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 24*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 
											   12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 
											   24*cpow(c13,3)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 6*c11*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) - 
											   33*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
											   12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 
											   24*cpow(c13,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) + 9*c11*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 
											   6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4) + 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4) - 
											   6*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4) + 6*c11*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4) - 
											   6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 
											   12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c13*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4) - 
											   12*c13*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) - 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) + 
											   3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4) - 6*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) + 
											   3*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) - 2*cpow(c13,6)*cpow(c44,3)*cpow(p1,6) + 6*c11*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,6) - 
											   6*cpow(c11,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,6) + 2*cpow(c11,3)*cpow(c33,3)*cpow(c44,3)*cpow(p1,6) - 
											   12*cpow(c13,5)*cpow(c44,3)*c55*cpow(p1,6) + 24*c11*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,6) - 
											   12*cpow(c11,2)*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,6) - 24*cpow(c13,4)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) + 
											   33*c11*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 9*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 
											   16*cpow(c13,3)*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) + 18*c11*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) - 
											   3*cpow(c13,4)*c33*cpow(c44,2)*c55*c66*cpow(p1,6) + 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,6) - 
											   3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,6) - 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) + 
											   12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) - 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 
											   18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) - 
											   3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) + 6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,6) + 
											   2*cpow(c33,3)*cpow(c55,3)*cpow(c66,3)*cpow(p1,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p2,2) - 
											   3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p2,2) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p2,2) + 
											   6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,2) - 6*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,2) + 
											   6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 
											   12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 12*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) + 
											   18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(p2,2) - 
											   3*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,2) - 6*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 
											   6*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 3*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 6*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 
											   9*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 12*c23*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) + 
											   6*cpow(c33,3)*cpow(c44,3)*c66*cpow(p2,2) - 6*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,2) - 6*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,2) - 
											   3*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,2) - 
											   3*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,2) - 3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
											   6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
											   6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 
											   6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) + 
											   18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 
											   3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
											   18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   12*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
											   18*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   6*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
											   18*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   12*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   12*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   24*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   9*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c11*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
											   18*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
											   12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
											   24*c13*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 12*c13*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
											   18*c12*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
											   9*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
											   18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 48*c13*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
											   18*c12*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
											   18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 12*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 
											   12*c11*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
											   18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) - 
											   18*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 3*cpow(c13,2)*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
											   6*c11*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 42*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
											   18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 
											   18*c12*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
											   18*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 36*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
											   6*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 6*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
											   24*c13*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
											   12*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 12*c22*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
											   3*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 6*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
											   42*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 24*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
											   18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 36*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 
											   3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) + 
											   6*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,4)*cpow(c23,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
											   3*cpow(c13,4)*c22*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 18*c12*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
											   3*c11*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
											   9*cpow(c12,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
											   6*c11*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
											   3*cpow(c11,2)*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
											   9*c11*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 3*cpow(c11,2)*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
											   12*cpow(c13,4)*c23*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 18*c12*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 
											   6*c11*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 
											   6*cpow(c11,2)*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 24*cpow(c13,3)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
											   12*cpow(c13,3)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 54*c12*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
											   6*c11*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
											   18*cpow(c12,2)*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
											   12*c11*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
											   18*c11*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 48*cpow(c13,3)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
											   54*c12*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 12*c11*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
											   18*c11*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 24*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
											   12*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 36*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
											   9*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 27*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
											   18*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 48*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
											   36*c12*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 18*c11*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
											   6*cpow(c13,4)*c33*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) - 12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 
											   6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 18*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 
											   18*c12*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 18*c11*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
											   18*c11*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 42*cpow(c13,3)*c33*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) - 
											   42*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) + 3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
											   6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
											   6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
											   6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 60*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
											   54*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
											   6*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 96*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
											   18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
											   12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
											   48*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 72*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
											   72*c13*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 3*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
											   3*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
											   18*c12*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
											   24*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
											   6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
											   30*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 36*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
											   6*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,3)*cpow(p1,4)*cpow(p2,2) + 3*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p2,4) - 
											   6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,4) + 3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,4) + 
											   12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) - 12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) + 
											   12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 
											   6*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,4) + 12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,4) - 
											   6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*cpow(p2,4) - 6*cpow(c23,4)*c44*cpow(c55,3)*cpow(p2,4) + 3*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,4) - 
											   24*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,4) + 3*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) + 
											   24*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) - 24*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 
											   6*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 33*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 
											   18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 24*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) + 
											   9*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 
											   6*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p2,4) + 6*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,4) - 
											   12*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,4) - 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,4) + 
											   6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,4) - 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
											   12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
											   12*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 
											   6*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p2,4) + 3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,4) + 
											   3*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c23,4)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
											   3*cpow(c13,2)*c22*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 18*c12*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   3*c11*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 3*cpow(c13,2)*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   18*c12*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   9*cpow(c12,2)*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
											   6*c11*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 9*cpow(c12,2)*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   3*c11*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 24*cpow(c13,2)*cpow(c23,3)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
											   6*cpow(c13,2)*c22*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
											   54*c12*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 12*c11*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   18*c12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   18*cpow(c12,2)*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
											   12*c11*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   24*cpow(c13,2)*cpow(c23,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 9*cpow(c13,2)*c22*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
											   36*c12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 12*c11*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   27*cpow(c12,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 18*c11*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   12*c13*cpow(c23,4)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 6*c13*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
											   18*c12*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 6*c13*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
											   18*c12*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 48*c13*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
											   12*c13*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 54*c12*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
											   18*c12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 48*c13*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
											   18*c13*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 36*c12*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
											   3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 
											   18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
											   9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
											   6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
											   12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 18*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
											   18*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 18*c12*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 
											   18*c12*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 60*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
											   6*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
											   54*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 48*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
											   72*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c23,4)*c33*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 
											   12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
											   42*cpow(c23,3)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 42*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
											   96*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
											   72*c23*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
											   6*c11*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
											   18*c12*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 30*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
											   3*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
											   3*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
											   24*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
											   36*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 6*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,3)*cpow(p1,2)*cpow(p2,4) - 
											   2*cpow(c23,6)*cpow(c55,3)*cpow(p2,6) + 6*c22*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,6) - 
											   6*cpow(c22,2)*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,6) + 2*cpow(c22,3)*cpow(c33,3)*cpow(c55,3)*cpow(p2,6) - 
											   12*cpow(c23,5)*c44*cpow(c55,3)*cpow(p2,6) + 24*c22*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,6) - 
											   12*cpow(c22,2)*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,6) - 24*cpow(c23,4)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) + 
											   33*c22*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 9*cpow(c22,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 
											   16*cpow(c23,3)*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) + 18*c22*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) - 
											   3*cpow(c23,4)*c33*c44*cpow(c55,2)*c66*cpow(p2,6) + 6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,6) - 
											   3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) + 
											   12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 
											   18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) - 
											   3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,6) + 
											   2*cpow(c33,3)*cpow(c44,3)*cpow(c66,3)*cpow(p2,6),2) + 
										      4*cpow(-cpow(c33*c44 + c33*c55 + c44*c55 + cpow(c13,2)*c44*cpow(p1,2) - c11*c33*c44*cpow(p1,2) + 2*c13*c44*c55*cpow(p1,2) - c33*c55*c66*cpow(p1,2) + 
												   cpow(c23,2)*c55*cpow(p2,2) - c22*c33*c55*cpow(p2,2) + 2*c23*c44*c55*cpow(p2,2) - c33*c44*c66*cpow(p2,2),2) - 
											     3*c33*c44*c55*(-c33 - c44 - c55 - cpow(c13,2)*cpow(p1,2) + c11*c33*cpow(p1,2) + c11*c44*cpow(p1,2) - 2*c13*c55*cpow(p1,2) + c44*c55*cpow(p1,2) + 
													    c33*c66*cpow(p1,2) + c55*c66*cpow(p1,2) - c11*c44*c55*cpow(p1,4) + cpow(c13,2)*c66*cpow(p1,4) - c11*c33*c66*cpow(p1,4) + 
													    2*c13*c55*c66*cpow(p1,4) - cpow(c23,2)*cpow(p2,2) + c22*c33*cpow(p2,2) - 2*c23*c44*cpow(p2,2) + c22*c55*cpow(p2,2) + c44*c55*cpow(p2,2) + 
													    c33*c66*cpow(p2,2) + c44*c66*cpow(p2,2) + cpow(c13,2)*c22*cpow(p1,2)*cpow(p2,2) - 2*c12*c13*c23*cpow(p1,2)*cpow(p2,2) + 
													    c11*cpow(c23,2)*cpow(p1,2)*cpow(p2,2) + cpow(c12,2)*c33*cpow(p1,2)*cpow(p2,2) - c11*c22*c33*cpow(p1,2)*cpow(p2,2) - 
													    2*c12*c13*c44*cpow(p1,2)*cpow(p2,2) + 2*c11*c23*c44*cpow(p1,2)*cpow(p2,2) + 2*c13*c22*c55*cpow(p1,2)*cpow(p2,2) - 
													    2*c12*c23*c55*cpow(p1,2)*cpow(p2,2) - 2*c12*c44*c55*cpow(p1,2)*cpow(p2,2) - 2*c13*c23*c66*cpow(p1,2)*cpow(p2,2) + 
													    2*c12*c33*c66*cpow(p1,2)*cpow(p2,2) - 2*c13*c44*c66*cpow(p1,2)*cpow(p2,2) - 2*c23*c55*c66*cpow(p1,2)*cpow(p2,2) - 
													    4*c44*c55*c66*cpow(p1,2)*cpow(p2,2) - c22*c44*c55*cpow(p2,4) + cpow(c23,2)*c66*cpow(p2,4) - c22*c33*c66*cpow(p2,4) + 2*c23*c44*c66*cpow(p2,4)),3)),
			     0.3333333333333333)) - cpow(-2*cpow(c33,3)*cpow(c44,3) + 3*cpow(c33,3)*cpow(c44,2)*c55 + 3*cpow(c33,2)*cpow(c44,3)*c55 + 3*cpow(c33,3)*c44*cpow(c55,2) - 
							 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2) + 3*c33*cpow(c44,3)*cpow(c55,2) - 2*cpow(c33,3)*cpow(c55,3) + 3*cpow(c33,2)*c44*cpow(c55,3) + 
							 3*c33*cpow(c44,2)*cpow(c55,3) - 2*cpow(c44,3)*cpow(c55,3) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,2) + 
							 6*c11*cpow(c33,3)*cpow(c44,3)*cpow(p1,2) + 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2) - 6*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2) - 
							 3*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,2) - 6*c11*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) - 12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) + 
							 3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2) - 3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2) + 
							 6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 
							 12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) - 6*cpow(c13,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 
							 3*c11*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 6*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 9*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) + 
							 6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2) + 12*c13*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) + 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) - 
							 12*c13*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 3*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2) - 
							 6*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2) + 6*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2) - 
							 6*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2) - 3*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2) - 6*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,4) + 
							 12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,4) - 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*cpow(p1,4) + 
							 3*cpow(c13,4)*c33*cpow(c44,2)*c55*cpow(p1,4) - 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4) + 
							 3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4) - 6*cpow(c13,4)*cpow(c44,3)*c55*cpow(p1,4) + 3*c11*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,4) - 
							 24*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4) + 3*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 
							 24*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 
							 12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 24*cpow(c13,3)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
							 6*c11*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) - 33*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
							 18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 
							 18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 24*cpow(c13,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) + 
							 9*c11*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4) + 
							 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4) - 6*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4) + 
							 6*c11*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4) - 6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 
							 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 
							 12*c13*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4) - 12*c13*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) - 
							 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) + 3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4) - 
							 6*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) + 3*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) - 2*cpow(c13,6)*cpow(c44,3)*cpow(p1,6) + 
							 6*c11*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,6) - 6*cpow(c11,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,6) + 
							 2*cpow(c11,3)*cpow(c33,3)*cpow(c44,3)*cpow(p1,6) - 12*cpow(c13,5)*cpow(c44,3)*c55*cpow(p1,6) + 24*c11*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,6) - 
							 12*cpow(c11,2)*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,6) - 24*cpow(c13,4)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) + 
							 33*c11*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 9*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 
							 16*cpow(c13,3)*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) + 18*c11*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) - 
							 3*cpow(c13,4)*c33*cpow(c44,2)*c55*c66*cpow(p1,6) + 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,6) - 
							 3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,6) - 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) + 
							 12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) - 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 
							 18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) - 
							 3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) + 6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,6) + 
							 2*cpow(c33,3)*cpow(c55,3)*cpow(c66,3)*cpow(p1,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p2,2) - 3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p2,2) + 
							 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p2,2) + 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,2) - 6*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,2) + 
							 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 
							 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 12*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) + 
							 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(p2,2) - 
							 3*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,2) - 6*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 
							 6*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 3*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 6*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 
							 9*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 12*c23*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) + 
							 6*cpow(c33,3)*cpow(c44,3)*c66*cpow(p2,2) - 6*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,2) - 6*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,2) - 
							 3*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,2) - 3*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,2) - 
							 3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 
							 18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
							 9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 
							 12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
							 6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
							 18*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
							 18*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 24*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
							 18*c12*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 9*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 12*c11*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
							 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 24*c13*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
							 12*c13*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
							 12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 9*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
							 18*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
							 48*c13*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
							 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
							 12*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 12*c11*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
							 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) - 18*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
							 3*cpow(c13,2)*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 6*c11*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
							 42*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 18*c12*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
							 6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 18*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
							 36*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 6*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
							 6*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 24*c13*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
							 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 12*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 
							 12*c22*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 3*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
							 6*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 42*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
							 24*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
							 36*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 
							 3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) + 6*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c13,4)*cpow(c23,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 3*cpow(c13,4)*c22*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
							 18*c12*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 3*c11*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
							 9*cpow(c12,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 6*c11*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
							 18*c11*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 3*cpow(c11,2)*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
							 9*c11*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 3*cpow(c11,2)*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
							 12*cpow(c13,4)*c23*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 18*c12*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 
							 6*c11*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 
							 6*cpow(c11,2)*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 24*cpow(c13,3)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
							 12*cpow(c13,3)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 54*c12*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
							 6*c11*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*cpow(c12,2)*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
							 12*c11*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
							 48*cpow(c13,3)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 54*c12*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
							 12*c11*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
							 24*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
							 36*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 9*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
							 27*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 18*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
							 48*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 36*c12*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
							 18*c11*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 6*cpow(c13,4)*c33*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) - 
							 12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 
							 18*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 
							 18*c11*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 18*c11*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
							 42*cpow(c13,3)*c33*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) - 42*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
							 3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
							 18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
							 9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
							 60*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 54*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
							 6*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 96*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
							 18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
							 12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
							 48*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 72*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
							 72*c13*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 3*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
							 3*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
							 18*c12*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 24*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
							 6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
							 30*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 36*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
							 6*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,3)*cpow(p1,4)*cpow(p2,2) + 3*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p2,4) - 
							 6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,4) + 3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,4) + 
							 12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) - 12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) + 
							 12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 
							 6*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,4) + 12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,4) - 6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*cpow(p2,4) - 
							 6*cpow(c23,4)*c44*cpow(c55,3)*cpow(p2,4) + 3*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,4) - 24*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,4) + 
							 3*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) + 24*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) - 
							 24*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 6*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 
							 33*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 
							 24*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) + 9*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 
							 6*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p2,4) + 6*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,4) - 
							 12*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,4) - 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,4) + 
							 6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,4) - 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
							 12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
							 12*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 
							 6*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p2,4) + 3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,4) + 
							 3*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c23,4)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
							 3*cpow(c13,2)*c22*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 18*c12*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
							 3*c11*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 3*cpow(c13,2)*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
							 18*c12*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 9*cpow(c12,2)*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
							 6*c11*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 9*cpow(c12,2)*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
							 3*c11*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 24*cpow(c13,2)*cpow(c23,3)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
							 6*cpow(c13,2)*c22*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 54*c12*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
							 12*c11*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
							 18*cpow(c12,2)*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
							 12*c11*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 24*cpow(c13,2)*cpow(c23,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
							 9*cpow(c13,2)*c22*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 36*c12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
							 12*c11*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 27*cpow(c12,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
							 18*c11*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 12*c13*cpow(c23,4)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
							 6*c13*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 18*c12*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
							 6*c13*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 18*c12*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
							 48*c13*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 12*c13*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
							 54*c12*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 18*c12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
							 48*c13*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 18*c13*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
							 36*c12*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
							 6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
							 6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 
							 6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) - 
							 18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
							 18*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 18*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
							 18*c12*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 18*c12*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 
							 60*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 6*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
							 54*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 48*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
							 72*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c23,4)*c33*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 
							 12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
							 42*cpow(c23,3)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 42*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
							 96*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
							 72*c23*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
							 6*c11*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
							 18*c12*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 30*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
							 3*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 3*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
							 24*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
							 36*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 6*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,3)*cpow(p1,2)*cpow(p2,4) - 
							 2*cpow(c23,6)*cpow(c55,3)*cpow(p2,6) + 6*c22*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,6) - 6*cpow(c22,2)*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,6) + 
							 2*cpow(c22,3)*cpow(c33,3)*cpow(c55,3)*cpow(p2,6) - 12*cpow(c23,5)*c44*cpow(c55,3)*cpow(p2,6) + 24*c22*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,6) - 
							 12*cpow(c22,2)*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,6) - 24*cpow(c23,4)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) + 
							 33*c22*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 9*cpow(c22,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 
							 16*cpow(c23,3)*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) + 18*c22*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) - 
							 3*cpow(c23,4)*c33*c44*cpow(c55,2)*c66*cpow(p2,6) + 6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,6) - 
							 3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) + 
							 12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 
							 18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) - 
							 3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,6) + 
							 2*cpow(c33,3)*cpow(c44,3)*cpow(c66,3)*cpow(p2,6) + csqrt(cpow(-2*cpow(c33,3)*cpow(c44,3) + 3*cpow(c33,3)*cpow(c44,2)*c55 + 
														       3*cpow(c33,2)*cpow(c44,3)*c55 + 3*cpow(c33,3)*c44*cpow(c55,2) - 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2) + 3*c33*cpow(c44,3)*cpow(c55,2) - 
														       2*cpow(c33,3)*cpow(c55,3) + 3*cpow(c33,2)*c44*cpow(c55,3) + 3*c33*cpow(c44,2)*cpow(c55,3) - 2*cpow(c44,3)*cpow(c55,3) - 
														       6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,2) + 6*c11*cpow(c33,3)*cpow(c44,3)*cpow(p1,2) + 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2) - 
														       6*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2) - 3*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,2) - 6*c11*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) - 
														       12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) + 3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2) - 3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2) + 
														       6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 
														       12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) - 6*cpow(c13,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 
														       3*c11*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 6*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 9*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) + 
														       6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2) + 12*c13*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) + 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) - 
														       12*c13*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 3*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2) - 
														       6*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2) + 6*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2) - 
														       6*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2) - 3*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2) - 6*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,4) + 
														       12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,4) - 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*cpow(p1,4) + 
														       3*cpow(c13,4)*c33*cpow(c44,2)*c55*cpow(p1,4) - 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4) + 
														       3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4) - 6*cpow(c13,4)*cpow(c44,3)*c55*cpow(p1,4) + 3*c11*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,4) - 
														       24*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4) + 3*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 
														       24*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 
														       12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 24*cpow(c13,3)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
														       6*c11*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) - 33*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
														       18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 
														       18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 24*cpow(c13,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) + 
														       9*c11*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 
														       6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4) + 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4) - 
														       6*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4) + 6*c11*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4) - 
														       6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 
														       12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c13*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4) - 
														       12*c13*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) - 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) + 
														       3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4) - 6*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) + 
														       3*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) - 2*cpow(c13,6)*cpow(c44,3)*cpow(p1,6) + 6*c11*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,6) - 
														       6*cpow(c11,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,6) + 2*cpow(c11,3)*cpow(c33,3)*cpow(c44,3)*cpow(p1,6) - 
														       12*cpow(c13,5)*cpow(c44,3)*c55*cpow(p1,6) + 24*c11*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,6) - 
														       12*cpow(c11,2)*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,6) - 24*cpow(c13,4)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) + 
														       33*c11*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 9*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 
														       16*cpow(c13,3)*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) + 18*c11*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) - 
														       3*cpow(c13,4)*c33*cpow(c44,2)*c55*c66*cpow(p1,6) + 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,6) - 
														       3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,6) - 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) + 
														       12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) - 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 
														       18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) - 
														       3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) + 6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,6) + 
														       2*cpow(c33,3)*cpow(c55,3)*cpow(c66,3)*cpow(p1,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p2,2) - 
														       3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p2,2) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p2,2) + 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,2) - 
														       6*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,2) + 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 
														       12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 
														       12*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) + 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) - 
														       6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(p2,2) - 3*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,2) - 
														       6*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 
														       3*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 6*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 9*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 
														       12*c23*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) + 6*cpow(c33,3)*cpow(c44,3)*c66*cpow(p2,2) - 
														       6*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,2) - 6*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,2) - 3*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,2) + 
														       12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,2) - 3*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,2) - 
														       3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 
														       18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
														       9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 
														       12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
														       6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
														       18*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
														       18*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 24*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
														       18*c12*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 9*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       12*c11*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
														       18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 24*c13*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
														       12*c13*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
														       12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 9*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
														       18*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
														       48*c13*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
														       18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
														       12*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 12*c11*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
														       6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) - 18*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
														       3*cpow(c13,2)*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 6*c11*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
														       42*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 18*c12*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
														       6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 18*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
														       36*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 6*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
														       6*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 24*c13*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
														       18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 12*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 
														       12*c22*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 3*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
														       6*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 42*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
														       24*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
														       36*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 
														       3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) + 6*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c13,4)*cpow(c23,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 3*cpow(c13,4)*c22*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
														       18*c12*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 3*c11*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
														       9*cpow(c12,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
														       6*c11*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
														       3*cpow(c11,2)*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 9*c11*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
														       3*cpow(c11,2)*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,4)*c23*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 
														       18*c12*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 6*c11*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 
														       18*c11*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 6*cpow(c11,2)*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 
														       24*cpow(c13,3)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,3)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
														       54*c12*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 6*c11*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
														       18*cpow(c12,2)*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
														       12*c11*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
														       48*cpow(c13,3)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 54*c12*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
														       12*c11*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
														       24*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
														       36*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 9*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
														       27*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 18*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
														       48*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 36*c12*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
														       18*c11*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 6*cpow(c13,4)*c33*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) - 
														       12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 
														       18*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 
														       18*c11*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 18*c11*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
														       42*cpow(c13,3)*c33*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) - 42*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
														       3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
														       18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
														       9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
														       60*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 54*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
														       6*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 96*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
														       18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
														       12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
														       48*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 72*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
														       72*c13*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 3*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
														       3*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
														       18*c12*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
														       24*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
														       6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
														       30*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 36*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
														       6*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,3)*cpow(p1,4)*cpow(p2,2) + 3*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p2,4) - 
														       6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,4) + 3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,4) + 
														       12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) - 12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) + 
														       12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 
														       6*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,4) + 12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,4) - 
														       6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*cpow(p2,4) - 6*cpow(c23,4)*c44*cpow(c55,3)*cpow(p2,4) + 3*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,4) - 
														       24*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,4) + 3*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) + 
														       24*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) - 24*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 
														       6*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 33*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 
														       18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 24*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) + 
														       9*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 
														       6*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p2,4) + 6*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,4) - 
														       12*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,4) - 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,4) + 
														       6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,4) - 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
														       12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
														       12*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 
														       6*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p2,4) + 3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,4) + 
														       3*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c23,4)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
														       3*cpow(c13,2)*c22*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 18*c12*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       3*c11*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 3*cpow(c13,2)*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       18*c12*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 9*cpow(c12,2)*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
														       6*c11*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 9*cpow(c12,2)*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       3*c11*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 24*cpow(c13,2)*cpow(c23,3)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
														       6*cpow(c13,2)*c22*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 54*c12*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       12*c11*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       18*cpow(c12,2)*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
														       12*c11*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       24*cpow(c13,2)*cpow(c23,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 9*cpow(c13,2)*c22*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
														       36*c12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 12*c11*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       27*cpow(c12,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 18*c11*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       12*c13*cpow(c23,4)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 6*c13*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
														       18*c12*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 6*c13*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
														       18*c12*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 48*c13*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
														       12*c13*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 54*c12*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
														       18*c12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 48*c13*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
														       18*c13*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 36*c12*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
														       3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 
														       18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
														       9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
														       6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
														       12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 18*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
														       18*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 18*c12*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 
														       18*c12*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 60*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
														       6*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 54*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 
														       48*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 72*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 
														       6*cpow(c23,4)*c33*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
														       6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 42*cpow(c23,3)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 
														       42*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 96*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 
														       18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 72*c23*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 
														       6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 6*c11*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
														       18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 18*c12*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
														       30*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 3*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
														       3*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 24*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
														       36*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 6*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,3)*cpow(p1,2)*cpow(p2,4) - 
														       2*cpow(c23,6)*cpow(c55,3)*cpow(p2,6) + 6*c22*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,6) - 6*cpow(c22,2)*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,6) + 
														       2*cpow(c22,3)*cpow(c33,3)*cpow(c55,3)*cpow(p2,6) - 12*cpow(c23,5)*c44*cpow(c55,3)*cpow(p2,6) + 24*c22*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,6) - 
														       12*cpow(c22,2)*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,6) - 24*cpow(c23,4)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) + 
														       33*c22*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 9*cpow(c22,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 
														       16*cpow(c23,3)*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) + 18*c22*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) - 
														       3*cpow(c23,4)*c33*c44*cpow(c55,2)*c66*cpow(p2,6) + 6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,6) - 
														       3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) + 
														       12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 
														       18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) - 
														       3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,6) + 
														       2*cpow(c33,3)*cpow(c44,3)*cpow(c66,3)*cpow(p2,6),2) + 4*
														  cpow(-cpow(c33*c44 + c33*c55 + c44*c55 + cpow(c13,2)*c44*cpow(p1,2) - c11*c33*c44*cpow(p1,2) + 2*c13*c44*c55*cpow(p1,2) - c33*c55*c66*cpow(p1,2) + 
															     cpow(c23,2)*c55*cpow(p2,2) - c22*c33*c55*cpow(p2,2) + 2*c23*c44*c55*cpow(p2,2) - c33*c44*c66*cpow(p2,2),2) - 
														       3*c33*c44*c55*(-c33 - c44 - c55 - cpow(c13,2)*cpow(p1,2) + c11*c33*cpow(p1,2) + c11*c44*cpow(p1,2) - 2*c13*c55*cpow(p1,2) + c44*c55*cpow(p1,2) + 
																      c33*c66*cpow(p1,2) + c55*c66*cpow(p1,2) - c11*c44*c55*cpow(p1,4) + cpow(c13,2)*c66*cpow(p1,4) - c11*c33*c66*cpow(p1,4) + 2*c13*c55*c66*cpow(p1,4) - 
																      cpow(c23,2)*cpow(p2,2) + c22*c33*cpow(p2,2) - 2*c23*c44*cpow(p2,2) + c22*c55*cpow(p2,2) + c44*c55*cpow(p2,2) + c33*c66*cpow(p2,2) + 
																      c44*c66*cpow(p2,2) + cpow(c13,2)*c22*cpow(p1,2)*cpow(p2,2) - 2*c12*c13*c23*cpow(p1,2)*cpow(p2,2) + c11*cpow(c23,2)*cpow(p1,2)*cpow(p2,2) + 
																      cpow(c12,2)*c33*cpow(p1,2)*cpow(p2,2) - c11*c22*c33*cpow(p1,2)*cpow(p2,2) - 2*c12*c13*c44*cpow(p1,2)*cpow(p2,2) + 
																      2*c11*c23*c44*cpow(p1,2)*cpow(p2,2) + 2*c13*c22*c55*cpow(p1,2)*cpow(p2,2) - 2*c12*c23*c55*cpow(p1,2)*cpow(p2,2) - 
																      2*c12*c44*c55*cpow(p1,2)*cpow(p2,2) - 2*c13*c23*c66*cpow(p1,2)*cpow(p2,2) + 2*c12*c33*c66*cpow(p1,2)*cpow(p2,2) - 
																      2*c13*c44*c66*cpow(p1,2)*cpow(p2,2) - 2*c23*c55*c66*cpow(p1,2)*cpow(p2,2) - 4*c44*c55*c66*cpow(p1,2)*cpow(p2,2) - c22*c44*c55*cpow(p2,4) + 
																      cpow(c23,2)*c66*cpow(p2,4) - c22*c33*c66*cpow(p2,4) + 2*c23*c44*c66*cpow(p2,4)),3)),0.3333333333333333)/(3.*cpow(2,0.3333333333333333)*c33*c44*c55);

/*dy ###################################################################################################################*/

    dy = p2*(-2*(c13 + c55)*(c13*c22 - c12*(c23 + c44) + c22*c55 - c23*c66 - c44*c66)*cpow(p1,2)*cpow(p3,2) - 
	     (c23 + c44)*cpow(p3,2)*(-((c13 + c55)*(c12 + c66)*cpow(p1,2)) + (c23 + c44)*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) - 
	     (c23 + c44)*cpow(p3,2)*(-((c13 + c55)*(c12 + c66)*cpow(p1,2)) + 2*(c23 + c44)*c66*cpow(p2,2) + 
				     (c23 + c44)*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) + 
	     2*c44*(-(cpow(c12 + c66,2)*cpow(p1,2)*cpow(p2,2)) + (-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))*
		    (-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) + 
	     2*(-1 + c55*cpow(p1,2) + c44*cpow(p2,2) + c33*cpow(p3,2))*(-(cpow(c12,2)*cpow(p1,2)) + c66*(-1 - 2*c12*cpow(p1,2) + c44*cpow(p3,2)) + 
									c22*(-1 + c11*cpow(p1,2) + 2*c66*cpow(p2,2) + c55*cpow(p3,2))));

/*dt ###################################################################################################################*/

    dt = cpow(p1,2)*(2*(c23 + c44)*(-(c11*(c23 + c44)) + c12*(c13 + c55) + (c13 + c55)*c66)*cpow(p2,2)*cpow(p3,2) + 
		     (c13 + c55)*cpow(p3,2)*((c23 + c44)*(c12 + c66)*cpow(p2,2) - (c13 + c55)*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))) + 
		     (c13 + c55)*cpow(p3,2)*(-2*(c13 + c55)*c66*cpow(p1,2) + (c23 + c44)*(c12 + c66)*cpow(p2,2) - 
					     (c13 + c55)*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))) + 
		     2*(-1 + c55*cpow(p1,2) + c44*cpow(p2,2) + c33*cpow(p3,2))*(-(cpow(c12,2)*cpow(p2,2)) + c11*(-1 + 2*c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2)) + 
										c66*(-1 - 2*c12*cpow(p2,2) + c55*cpow(p3,2))) + 2*c55*(-(cpow(c12 + c66,2)*cpow(p1,2)*cpow(p2,2)) + 
																       (-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2)))) + 
	cpow(p2,2)*(-2*(c13 + c55)*(c13*c22 - c12*(c23 + c44) + c22*c55 - c23*c66 - c44*c66)*cpow(p1,2)*cpow(p3,2) - 
		    (c23 + c44)*cpow(p3,2)*(-((c13 + c55)*(c12 + c66)*cpow(p1,2)) + (c23 + c44)*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) - 
		    (c23 + c44)*cpow(p3,2)*(-((c13 + c55)*(c12 + c66)*cpow(p1,2)) + 2*(c23 + c44)*c66*cpow(p2,2) + 
					    (c23 + c44)*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) + 
		    2*c44*(-(cpow(c12 + c66,2)*cpow(p1,2)*cpow(p2,2)) + (-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))*
			   (-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) + 
		    2*(-1 + c55*cpow(p1,2) + c44*cpow(p2,2) + c33*cpow(p3,2))*(-(cpow(c12,2)*cpow(p1,2)) + c66*(-1 - 2*c12*cpow(p1,2) + c44*cpow(p3,2)) + 
									       c22*(-1 + c11*cpow(p1,2) + 2*c66*cpow(p2,2) + c55*cpow(p3,2)))) + 
	cpow(p3,2)*((c13 + c55)*cpow(p1,2)*((c23 + c44)*(c12 + c66)*cpow(p2,2) - (c13 + c55)*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))) + 
		    (c13 + c55)*cpow(p1,2)*((c23 + c44)*(c12 + c66)*cpow(p2,2) - 2*c44*(c13 + c55)*cpow(p3,2) - 
					    (c13 + c55)*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))) - 
		    (c23 + c44)*cpow(p2,2)*(-((c13 + c55)*(c12 + c66)*cpow(p1,2)) + (c23 + c44)*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) - 
		    (c23 + c44)*cpow(p2,2)*(-((c13 + c55)*(c12 + c66)*cpow(p1,2)) + 2*(c23 + c44)*c55*cpow(p3,2) + 
					    (c23 + c44)*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) + 
		    2*c33*(-(cpow(c12 + c66,2)*cpow(p1,2)*cpow(p2,2)) + (-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))*
			   (-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) + 
		    2*(-1 + c55*cpow(p1,2) + c44*cpow(p2,2) + c33*cpow(p3,2))*(c55*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2)) + 
									       c44*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + 2*c55*cpow(p3,2))));
         
         
    return creal(dy)/creal(dt);

}


float zPvelo(float c11, float c22,float c33,float c44,float c55,float c66,float c12,float c13,float c23, float p1, float p2) {
/*<Compute qP-wave group velocity (z) given cij and px (p1) and py (p2) for layered orthorhombic media>*/

    float complex p3, dz, dt;


/*p3 ###################################################################################################################*/

    p3 = (c33*c44 + c33*c55 + c44*c55 + cpow(c13,2)*c44*cpow(p1,2) - c11*c33*c44*cpow(p1,2) + 2*c13*c44*c55*cpow(p1,2) - c33*c55*c66*cpow(p1,2) + 
	  cpow(c23,2)*c55*cpow(p2,2) - c22*c33*c55*cpow(p2,2) + 2*c23*c44*c55*cpow(p2,2) - c33*c44*c66*cpow(p2,2))/(3.*c33*c44*c55) + 
	(cpow(2,0.3333333333333333)*(-cpow(c33*c44 + c33*c55 + c44*c55 + cpow(c13,2)*c44*cpow(p1,2) - c11*c33*c44*cpow(p1,2) + 2*c13*c44*c55*cpow(p1,2) - 
					   c33*c55*c66*cpow(p1,2) + cpow(c23,2)*c55*cpow(p2,2) - c22*c33*c55*cpow(p2,2) + 2*c23*c44*c55*cpow(p2,2) - c33*c44*c66*cpow(p2,2),2) - 
				     3*c33*c44*c55*(-c33 - c44 - c55 - cpow(c13,2)*cpow(p1,2) + c11*c33*cpow(p1,2) + c11*c44*cpow(p1,2) - 2*c13*c55*cpow(p1,2) + c44*c55*cpow(p1,2) + 
						    c33*c66*cpow(p1,2) + c55*c66*cpow(p1,2) - c11*c44*c55*cpow(p1,4) + cpow(c13,2)*c66*cpow(p1,4) - c11*c33*c66*cpow(p1,4) + 2*c13*c55*c66*cpow(p1,4) - 
						    cpow(c23,2)*cpow(p2,2) + c22*c33*cpow(p2,2) - 2*c23*c44*cpow(p2,2) + c22*c55*cpow(p2,2) + c44*c55*cpow(p2,2) + c33*c66*cpow(p2,2) + 
						    c44*c66*cpow(p2,2) + cpow(c13,2)*c22*cpow(p1,2)*cpow(p2,2) - 2*c12*c13*c23*cpow(p1,2)*cpow(p2,2) + c11*cpow(c23,2)*cpow(p1,2)*cpow(p2,2) + 
						    cpow(c12,2)*c33*cpow(p1,2)*cpow(p2,2) - c11*c22*c33*cpow(p1,2)*cpow(p2,2) - 2*c12*c13*c44*cpow(p1,2)*cpow(p2,2) + 
						    2*c11*c23*c44*cpow(p1,2)*cpow(p2,2) + 2*c13*c22*c55*cpow(p1,2)*cpow(p2,2) - 2*c12*c23*c55*cpow(p1,2)*cpow(p2,2) - 
						    2*c12*c44*c55*cpow(p1,2)*cpow(p2,2) - 2*c13*c23*c66*cpow(p1,2)*cpow(p2,2) + 2*c12*c33*c66*cpow(p1,2)*cpow(p2,2) - 
						    2*c13*c44*c66*cpow(p1,2)*cpow(p2,2) - 2*c23*c55*c66*cpow(p1,2)*cpow(p2,2) - 4*c44*c55*c66*cpow(p1,2)*cpow(p2,2) - c22*c44*c55*cpow(p2,4) + 
						    cpow(c23,2)*c66*cpow(p2,4) - c22*c33*c66*cpow(p2,4) + 2*c23*c44*c66*cpow(p2,4))))/
	(3.*c33*c44*c55*cpow(-2*cpow(c33,3)*cpow(c44,3) + 3*cpow(c33,3)*cpow(c44,2)*c55 + 3*cpow(c33,2)*cpow(c44,3)*c55 + 3*cpow(c33,3)*c44*cpow(c55,2) - 
			     12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2) + 3*c33*cpow(c44,3)*cpow(c55,2) - 2*cpow(c33,3)*cpow(c55,3) + 3*cpow(c33,2)*c44*cpow(c55,3) + 
			     3*c33*cpow(c44,2)*cpow(c55,3) - 2*cpow(c44,3)*cpow(c55,3) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,2) + 
			     6*c11*cpow(c33,3)*cpow(c44,3)*cpow(p1,2) + 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2) - 6*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2) - 
			     3*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,2) - 6*c11*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) - 12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) + 
			     3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2) - 3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2) + 
			     6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 
			     12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) - 6*cpow(c13,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 
			     3*c11*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 6*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 9*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) + 
			     6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2) + 12*c13*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) + 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) - 
			     12*c13*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 3*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2) - 
			     6*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2) + 6*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2) - 
			     6*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2) - 3*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2) - 6*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,4) + 
			     12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,4) - 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*cpow(p1,4) + 
			     3*cpow(c13,4)*c33*cpow(c44,2)*c55*cpow(p1,4) - 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4) + 
			     3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4) - 6*cpow(c13,4)*cpow(c44,3)*c55*cpow(p1,4) + 3*c11*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,4) - 
			     24*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4) + 3*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 
			     24*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 
			     12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 24*cpow(c13,3)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
			     6*c11*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) - 33*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
			     18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 
			     18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 24*cpow(c13,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) + 
			     9*c11*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 
			     6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4) + 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4) - 
			     6*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4) + 6*c11*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4) - 
			     6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 
			     12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c13*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4) - 
			     12*c13*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) - 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) + 
			     3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4) - 6*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) + 
			     3*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) - 2*cpow(c13,6)*cpow(c44,3)*cpow(p1,6) + 6*c11*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,6) - 
			     6*cpow(c11,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,6) + 2*cpow(c11,3)*cpow(c33,3)*cpow(c44,3)*cpow(p1,6) - 
			     12*cpow(c13,5)*cpow(c44,3)*c55*cpow(p1,6) + 24*c11*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,6) - 
			     12*cpow(c11,2)*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,6) - 24*cpow(c13,4)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) + 
			     33*c11*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 9*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 
			     16*cpow(c13,3)*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) + 18*c11*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) - 
			     3*cpow(c13,4)*c33*cpow(c44,2)*c55*c66*cpow(p1,6) + 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,6) - 
			     3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,6) - 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) + 
			     12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) - 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 
			     18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) - 
			     3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) + 6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,6) + 
			     2*cpow(c33,3)*cpow(c55,3)*cpow(c66,3)*cpow(p1,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p2,2) - 
			     3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p2,2) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p2,2) + 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,2) - 
			     6*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,2) + 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 
			     12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 
			     12*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) + 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,2) + 
			     6*c22*cpow(c33,3)*cpow(c55,3)*cpow(p2,2) - 3*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,2) - 6*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 
			     12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 3*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 
			     6*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 9*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 12*c23*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) - 
			     9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) + 6*cpow(c33,3)*cpow(c44,3)*c66*cpow(p2,2) - 6*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,2) - 
			     6*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,2) - 3*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,2) - 
			     3*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,2) - 3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
			     6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
			     6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 
			     6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) + 
			     18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 
			     3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
			     18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     12*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
			     18*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     6*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
			     18*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     12*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     12*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     24*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     9*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c11*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
			     12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
			     18*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
			     12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
			     24*c13*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 12*c13*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
			     18*c12*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
			     9*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
			     18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 48*c13*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
			     18*c12*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
			     18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 12*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 
			     12*c11*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
			     18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) - 
			     18*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 3*cpow(c13,2)*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
			     6*c11*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 42*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
			     18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 
			     18*c12*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
			     18*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 36*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
			     6*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 6*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
			     24*c13*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
			     12*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 12*c22*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
			     3*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 6*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
			     42*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 24*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
			     18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 36*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 
			     3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) + 
			     6*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,4)*cpow(c23,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
			     3*cpow(c13,4)*c22*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 18*c12*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
			     3*c11*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
			     9*cpow(c12,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
			     6*c11*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
			     3*cpow(c11,2)*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 9*c11*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
			     3*cpow(c11,2)*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,4)*c23*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 
			     18*c12*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 6*c11*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 
			     18*c11*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 6*cpow(c11,2)*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 
			     24*cpow(c13,3)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,3)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
			     54*c12*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 6*c11*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
			     18*cpow(c12,2)*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
			     12*c11*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
			     48*cpow(c13,3)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 54*c12*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
			     12*c11*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
			     24*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
			     36*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 9*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
			     27*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 18*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
			     48*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 36*c12*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
			     18*c11*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 6*cpow(c13,4)*c33*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) - 
			     12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 
			     18*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 
			     18*c11*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 18*c11*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
			     42*cpow(c13,3)*c33*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) - 42*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
			     3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
			     18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
			     9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
			     60*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 54*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
			     6*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 96*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
			     18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
			     12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
			     48*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 72*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
			     72*c13*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 3*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
			     3*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
			     18*c12*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 24*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
			     6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
			     30*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 36*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
			     6*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,3)*cpow(p1,4)*cpow(p2,2) + 3*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p2,4) - 
			     6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,4) + 3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,4) + 
			     12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) - 12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) + 
			     12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 
			     6*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,4) + 12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,4) - 6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*cpow(p2,4) - 
			     6*cpow(c23,4)*c44*cpow(c55,3)*cpow(p2,4) + 3*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,4) - 24*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,4) + 
			     3*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) + 24*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) - 
			     24*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 6*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 
			     33*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 
			     24*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) + 9*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 
			     6*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p2,4) + 6*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,4) - 
			     12*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,4) - 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,4) + 
			     6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,4) - 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
			     12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
			     12*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 
			     6*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p2,4) + 3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,4) + 
			     3*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c23,4)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
			     3*cpow(c13,2)*c22*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 18*c12*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
			     3*c11*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 3*cpow(c13,2)*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
			     18*c12*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 9*cpow(c12,2)*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
			     6*c11*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 9*cpow(c12,2)*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
			     3*c11*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 24*cpow(c13,2)*cpow(c23,3)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
			     6*cpow(c13,2)*c22*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 54*c12*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
			     12*c11*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
			     18*cpow(c12,2)*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
			     12*c11*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 24*cpow(c13,2)*cpow(c23,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
			     9*cpow(c13,2)*c22*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 36*c12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
			     12*c11*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 27*cpow(c12,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
			     18*c11*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 12*c13*cpow(c23,4)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
			     6*c13*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 18*c12*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
			     6*c13*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 18*c12*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
			     48*c13*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 12*c13*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
			     54*c12*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 18*c12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
			     48*c13*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 18*c13*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
			     36*c12*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
			     6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
			     6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 
			     6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) - 
			     18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
			     18*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 18*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
			     18*c12*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 18*c12*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 
			     60*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 6*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
			     54*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 48*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
			     72*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c23,4)*c33*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 
			     12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
			     42*cpow(c23,3)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 42*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
			     96*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
			     72*c23*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
			     6*c11*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
			     18*c12*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 30*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
			     3*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 3*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
			     24*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
			     36*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 6*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,3)*cpow(p1,2)*cpow(p2,4) - 
			     2*cpow(c23,6)*cpow(c55,3)*cpow(p2,6) + 6*c22*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,6) - 6*cpow(c22,2)*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,6) + 
			     2*cpow(c22,3)*cpow(c33,3)*cpow(c55,3)*cpow(p2,6) - 12*cpow(c23,5)*c44*cpow(c55,3)*cpow(p2,6) + 24*c22*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,6) - 
			     12*cpow(c22,2)*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,6) - 24*cpow(c23,4)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) + 
			     33*c22*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 9*cpow(c22,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 
			     16*cpow(c23,3)*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) + 18*c22*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) - 
			     3*cpow(c23,4)*c33*c44*cpow(c55,2)*c66*cpow(p2,6) + 6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,6) - 
			     3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) + 
			     12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 
			     18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) - 
			     3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,6) + 
			     2*cpow(c33,3)*cpow(c44,3)*cpow(c66,3)*cpow(p2,6) + csqrt(cpow(-2*cpow(c33,3)*cpow(c44,3) + 3*cpow(c33,3)*cpow(c44,2)*c55 + 
											   3*cpow(c33,2)*cpow(c44,3)*c55 + 3*cpow(c33,3)*c44*cpow(c55,2) - 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2) + 3*c33*cpow(c44,3)*cpow(c55,2) - 
											   2*cpow(c33,3)*cpow(c55,3) + 3*cpow(c33,2)*c44*cpow(c55,3) + 3*c33*cpow(c44,2)*cpow(c55,3) - 2*cpow(c44,3)*cpow(c55,3) - 
											   6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,2) + 6*c11*cpow(c33,3)*cpow(c44,3)*cpow(p1,2) + 
											   6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2) - 6*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2) - 
											   3*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,2) - 6*c11*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) - 12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) + 
											   3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2) - 3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2) + 
											   6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 
											   12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) - 6*cpow(c13,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 
											   3*c11*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 6*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 9*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) + 
											   6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2) + 12*c13*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) + 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) - 
											   12*c13*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 3*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2) - 
											   6*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2) + 6*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2) - 
											   6*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2) - 3*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2) - 6*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,4) + 
											   12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,4) - 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*cpow(p1,4) + 
											   3*cpow(c13,4)*c33*cpow(c44,2)*c55*cpow(p1,4) - 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4) + 
											   3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4) - 6*cpow(c13,4)*cpow(c44,3)*c55*cpow(p1,4) + 
											   3*c11*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,4) - 24*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4) + 
											   3*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 24*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 
											   12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 
											   24*cpow(c13,3)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 6*c11*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) - 
											   33*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
											   12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 
											   24*cpow(c13,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) + 9*c11*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 
											   6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4) + 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4) - 
											   6*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4) + 6*c11*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4) - 
											   6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 
											   12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c13*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4) - 
											   12*c13*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) - 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) + 
											   3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4) - 6*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) + 
											   3*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) - 2*cpow(c13,6)*cpow(c44,3)*cpow(p1,6) + 6*c11*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,6) - 
											   6*cpow(c11,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,6) + 2*cpow(c11,3)*cpow(c33,3)*cpow(c44,3)*cpow(p1,6) - 
											   12*cpow(c13,5)*cpow(c44,3)*c55*cpow(p1,6) + 24*c11*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,6) - 
											   12*cpow(c11,2)*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,6) - 24*cpow(c13,4)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) + 
											   33*c11*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 9*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 
											   16*cpow(c13,3)*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) + 18*c11*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) - 
											   3*cpow(c13,4)*c33*cpow(c44,2)*c55*c66*cpow(p1,6) + 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,6) - 
											   3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,6) - 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) + 
											   12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) - 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 
											   18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) - 
											   3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) + 6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,6) + 
											   2*cpow(c33,3)*cpow(c55,3)*cpow(c66,3)*cpow(p1,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p2,2) - 
											   3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p2,2) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p2,2) + 
											   6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,2) - 6*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,2) + 
											   6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 
											   12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 12*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) + 
											   18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(p2,2) - 
											   3*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,2) - 6*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 
											   6*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 3*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 6*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 
											   9*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 12*c23*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) + 
											   6*cpow(c33,3)*cpow(c44,3)*c66*cpow(p2,2) - 6*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,2) - 6*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,2) - 
											   3*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,2) - 
											   3*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,2) - 3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
											   6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
											   6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 
											   6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) + 
											   18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 
											   3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
											   18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   12*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
											   18*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   6*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
											   18*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   12*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   12*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   24*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   9*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c11*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
											   12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
											   18*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
											   12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
											   24*c13*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 12*c13*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
											   18*c12*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
											   9*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
											   18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 48*c13*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
											   18*c12*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
											   18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 12*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 
											   12*c11*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
											   18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) - 
											   18*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 3*cpow(c13,2)*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
											   6*c11*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 42*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
											   18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 
											   18*c12*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
											   18*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 36*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
											   6*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 6*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
											   24*c13*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
											   12*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 12*c22*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
											   3*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 6*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
											   42*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 24*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
											   18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 36*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 
											   3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) + 
											   6*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,4)*cpow(c23,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
											   3*cpow(c13,4)*c22*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 18*c12*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
											   3*c11*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
											   9*cpow(c12,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
											   6*c11*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
											   3*cpow(c11,2)*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
											   9*c11*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 3*cpow(c11,2)*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
											   12*cpow(c13,4)*c23*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 18*c12*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 
											   6*c11*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 
											   6*cpow(c11,2)*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 24*cpow(c13,3)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
											   12*cpow(c13,3)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 54*c12*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
											   6*c11*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
											   18*cpow(c12,2)*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
											   12*c11*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
											   18*c11*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 48*cpow(c13,3)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
											   54*c12*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 12*c11*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
											   18*c11*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 24*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
											   12*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 36*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
											   9*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 27*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
											   18*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 48*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
											   36*c12*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 18*c11*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
											   6*cpow(c13,4)*c33*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) - 12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 
											   6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 18*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 
											   18*c12*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 18*c11*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
											   18*c11*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 42*cpow(c13,3)*c33*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) - 
											   42*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) + 3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
											   6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
											   6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
											   6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 60*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
											   54*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
											   6*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 96*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
											   18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
											   12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
											   48*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 72*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
											   72*c13*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 3*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
											   3*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
											   18*c12*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
											   24*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
											   6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
											   30*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 36*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
											   6*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,3)*cpow(p1,4)*cpow(p2,2) + 3*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p2,4) - 
											   6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,4) + 3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,4) + 
											   12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) - 12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) + 
											   12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 
											   6*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,4) + 12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,4) - 
											   6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*cpow(p2,4) - 6*cpow(c23,4)*c44*cpow(c55,3)*cpow(p2,4) + 3*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,4) - 
											   24*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,4) + 3*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) + 
											   24*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) - 24*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 
											   6*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 33*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 
											   18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 24*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) + 
											   9*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 
											   6*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p2,4) + 6*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,4) - 
											   12*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,4) - 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,4) + 
											   6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,4) - 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
											   12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
											   12*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 
											   6*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p2,4) + 3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,4) + 
											   3*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c23,4)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
											   3*cpow(c13,2)*c22*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 18*c12*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   3*c11*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 3*cpow(c13,2)*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   18*c12*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   9*cpow(c12,2)*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
											   6*c11*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 9*cpow(c12,2)*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   3*c11*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 24*cpow(c13,2)*cpow(c23,3)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
											   6*cpow(c13,2)*c22*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
											   54*c12*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 12*c11*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   18*c12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   18*cpow(c12,2)*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
											   12*c11*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   24*cpow(c13,2)*cpow(c23,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 9*cpow(c13,2)*c22*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
											   36*c12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 12*c11*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   27*cpow(c12,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 18*c11*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
											   12*c13*cpow(c23,4)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 6*c13*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
											   18*c12*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 6*c13*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
											   18*c12*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 48*c13*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
											   12*c13*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 54*c12*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
											   18*c12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 48*c13*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
											   18*c13*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 36*c12*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
											   3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 
											   18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
											   9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
											   6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
											   12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 18*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
											   18*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 18*c12*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 
											   18*c12*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 60*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
											   6*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
											   54*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 48*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
											   72*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c23,4)*c33*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 
											   12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
											   42*cpow(c23,3)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 42*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
											   96*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
											   72*c23*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
											   6*c11*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
											   18*c12*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 30*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
											   3*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
											   3*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
											   24*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
											   36*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 6*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,3)*cpow(p1,2)*cpow(p2,4) - 
											   2*cpow(c23,6)*cpow(c55,3)*cpow(p2,6) + 6*c22*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,6) - 
											   6*cpow(c22,2)*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,6) + 2*cpow(c22,3)*cpow(c33,3)*cpow(c55,3)*cpow(p2,6) - 
											   12*cpow(c23,5)*c44*cpow(c55,3)*cpow(p2,6) + 24*c22*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,6) - 
											   12*cpow(c22,2)*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,6) - 24*cpow(c23,4)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) + 
											   33*c22*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 9*cpow(c22,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 
											   16*cpow(c23,3)*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) + 18*c22*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) - 
											   3*cpow(c23,4)*c33*c44*cpow(c55,2)*c66*cpow(p2,6) + 6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,6) - 
											   3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) + 
											   12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 
											   18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) - 
											   3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,6) + 
											   2*cpow(c33,3)*cpow(c44,3)*cpow(c66,3)*cpow(p2,6),2) + 
										      4*cpow(-cpow(c33*c44 + c33*c55 + c44*c55 + cpow(c13,2)*c44*cpow(p1,2) - c11*c33*c44*cpow(p1,2) + 2*c13*c44*c55*cpow(p1,2) - c33*c55*c66*cpow(p1,2) + 
												   cpow(c23,2)*c55*cpow(p2,2) - c22*c33*c55*cpow(p2,2) + 2*c23*c44*c55*cpow(p2,2) - c33*c44*c66*cpow(p2,2),2) - 
											     3*c33*c44*c55*(-c33 - c44 - c55 - cpow(c13,2)*cpow(p1,2) + c11*c33*cpow(p1,2) + c11*c44*cpow(p1,2) - 2*c13*c55*cpow(p1,2) + c44*c55*cpow(p1,2) + 
													    c33*c66*cpow(p1,2) + c55*c66*cpow(p1,2) - c11*c44*c55*cpow(p1,4) + cpow(c13,2)*c66*cpow(p1,4) - c11*c33*c66*cpow(p1,4) + 
													    2*c13*c55*c66*cpow(p1,4) - cpow(c23,2)*cpow(p2,2) + c22*c33*cpow(p2,2) - 2*c23*c44*cpow(p2,2) + c22*c55*cpow(p2,2) + c44*c55*cpow(p2,2) + 
													    c33*c66*cpow(p2,2) + c44*c66*cpow(p2,2) + cpow(c13,2)*c22*cpow(p1,2)*cpow(p2,2) - 2*c12*c13*c23*cpow(p1,2)*cpow(p2,2) + 
													    c11*cpow(c23,2)*cpow(p1,2)*cpow(p2,2) + cpow(c12,2)*c33*cpow(p1,2)*cpow(p2,2) - c11*c22*c33*cpow(p1,2)*cpow(p2,2) - 
													    2*c12*c13*c44*cpow(p1,2)*cpow(p2,2) + 2*c11*c23*c44*cpow(p1,2)*cpow(p2,2) + 2*c13*c22*c55*cpow(p1,2)*cpow(p2,2) - 
													    2*c12*c23*c55*cpow(p1,2)*cpow(p2,2) - 2*c12*c44*c55*cpow(p1,2)*cpow(p2,2) - 2*c13*c23*c66*cpow(p1,2)*cpow(p2,2) + 
													    2*c12*c33*c66*cpow(p1,2)*cpow(p2,2) - 2*c13*c44*c66*cpow(p1,2)*cpow(p2,2) - 2*c23*c55*c66*cpow(p1,2)*cpow(p2,2) - 
													    4*c44*c55*c66*cpow(p1,2)*cpow(p2,2) - c22*c44*c55*cpow(p2,4) + cpow(c23,2)*c66*cpow(p2,4) - c22*c33*c66*cpow(p2,4) + 2*c23*c44*c66*cpow(p2,4)),3)),
			     0.3333333333333333)) - cpow(-2*cpow(c33,3)*cpow(c44,3) + 3*cpow(c33,3)*cpow(c44,2)*c55 + 3*cpow(c33,2)*cpow(c44,3)*c55 + 3*cpow(c33,3)*c44*cpow(c55,2) - 
							 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2) + 3*c33*cpow(c44,3)*cpow(c55,2) - 2*cpow(c33,3)*cpow(c55,3) + 3*cpow(c33,2)*c44*cpow(c55,3) + 
							 3*c33*cpow(c44,2)*cpow(c55,3) - 2*cpow(c44,3)*cpow(c55,3) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,2) + 
							 6*c11*cpow(c33,3)*cpow(c44,3)*cpow(p1,2) + 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2) - 6*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2) - 
							 3*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,2) - 6*c11*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) - 12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) + 
							 3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2) - 3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2) + 
							 6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 
							 12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) - 6*cpow(c13,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 
							 3*c11*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 6*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 9*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) + 
							 6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2) + 12*c13*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) + 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) - 
							 12*c13*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 3*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2) - 
							 6*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2) + 6*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2) - 
							 6*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2) - 3*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2) - 6*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,4) + 
							 12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,4) - 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*cpow(p1,4) + 
							 3*cpow(c13,4)*c33*cpow(c44,2)*c55*cpow(p1,4) - 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4) + 
							 3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4) - 6*cpow(c13,4)*cpow(c44,3)*c55*cpow(p1,4) + 3*c11*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,4) - 
							 24*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4) + 3*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 
							 24*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 
							 12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 24*cpow(c13,3)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
							 6*c11*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) - 33*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
							 18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 
							 18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 24*cpow(c13,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) + 
							 9*c11*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4) + 
							 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4) - 6*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4) + 
							 6*c11*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4) - 6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 
							 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 
							 12*c13*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4) - 12*c13*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) - 
							 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) + 3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4) - 
							 6*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) + 3*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) - 2*cpow(c13,6)*cpow(c44,3)*cpow(p1,6) + 
							 6*c11*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,6) - 6*cpow(c11,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,6) + 
							 2*cpow(c11,3)*cpow(c33,3)*cpow(c44,3)*cpow(p1,6) - 12*cpow(c13,5)*cpow(c44,3)*c55*cpow(p1,6) + 24*c11*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,6) - 
							 12*cpow(c11,2)*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,6) - 24*cpow(c13,4)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) + 
							 33*c11*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 9*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 
							 16*cpow(c13,3)*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) + 18*c11*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) - 
							 3*cpow(c13,4)*c33*cpow(c44,2)*c55*c66*cpow(p1,6) + 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,6) - 
							 3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,6) - 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) + 
							 12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) - 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 
							 18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) - 
							 3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) + 6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,6) + 
							 2*cpow(c33,3)*cpow(c55,3)*cpow(c66,3)*cpow(p1,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p2,2) - 3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p2,2) + 
							 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p2,2) + 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,2) - 6*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,2) + 
							 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 
							 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 12*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) + 
							 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(p2,2) - 
							 3*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,2) - 6*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 
							 6*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 3*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 6*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 
							 9*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 12*c23*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) + 
							 6*cpow(c33,3)*cpow(c44,3)*c66*cpow(p2,2) - 6*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,2) - 6*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,2) - 
							 3*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,2) - 3*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,2) - 
							 3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 
							 18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
							 9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 
							 12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
							 6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
							 18*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
							 18*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 24*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
							 18*c12*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 9*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 12*c11*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
							 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 24*c13*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
							 12*c13*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
							 12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 9*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
							 18*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
							 48*c13*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
							 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
							 12*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 12*c11*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
							 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) - 18*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
							 3*cpow(c13,2)*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 6*c11*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
							 42*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 18*c12*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
							 6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 18*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
							 36*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 6*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
							 6*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 24*c13*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
							 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 12*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 
							 12*c22*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 3*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
							 6*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 42*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
							 24*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
							 36*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 
							 3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) + 6*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 
							 6*cpow(c13,4)*cpow(c23,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 3*cpow(c13,4)*c22*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
							 18*c12*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 3*c11*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
							 9*cpow(c12,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 6*c11*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
							 18*c11*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 3*cpow(c11,2)*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
							 9*c11*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 3*cpow(c11,2)*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
							 12*cpow(c13,4)*c23*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 18*c12*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 
							 6*c11*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 
							 6*cpow(c11,2)*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 24*cpow(c13,3)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
							 12*cpow(c13,3)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 54*c12*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
							 6*c11*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*cpow(c12,2)*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
							 12*c11*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
							 48*cpow(c13,3)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 54*c12*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
							 12*c11*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
							 24*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
							 36*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 9*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
							 27*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 18*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
							 48*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 36*c12*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
							 18*c11*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 6*cpow(c13,4)*c33*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) - 
							 12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 
							 18*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 
							 18*c11*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 18*c11*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
							 42*cpow(c13,3)*c33*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) - 42*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
							 3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
							 18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
							 9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
							 60*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 54*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
							 6*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 96*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
							 18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
							 12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
							 48*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 72*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
							 72*c13*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 3*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
							 3*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
							 18*c12*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 24*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
							 6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
							 30*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 36*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
							 6*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,3)*cpow(p1,4)*cpow(p2,2) + 3*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p2,4) - 
							 6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,4) + 3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,4) + 
							 12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) - 12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) + 
							 12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 
							 6*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,4) + 12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,4) - 6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*cpow(p2,4) - 
							 6*cpow(c23,4)*c44*cpow(c55,3)*cpow(p2,4) + 3*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,4) - 24*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,4) + 
							 3*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) + 24*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) - 
							 24*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 6*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 
							 33*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 
							 24*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) + 9*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 
							 6*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p2,4) + 6*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,4) - 
							 12*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,4) - 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,4) + 
							 6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,4) - 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
							 12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
							 12*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 
							 6*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p2,4) + 3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,4) + 
							 3*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c23,4)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
							 3*cpow(c13,2)*c22*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 18*c12*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
							 3*c11*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 3*cpow(c13,2)*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
							 18*c12*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 9*cpow(c12,2)*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
							 6*c11*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 9*cpow(c12,2)*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
							 3*c11*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 24*cpow(c13,2)*cpow(c23,3)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
							 6*cpow(c13,2)*c22*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 54*c12*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
							 12*c11*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
							 18*cpow(c12,2)*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
							 12*c11*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 24*cpow(c13,2)*cpow(c23,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
							 9*cpow(c13,2)*c22*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 36*c12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
							 12*c11*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 27*cpow(c12,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
							 18*c11*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 12*c13*cpow(c23,4)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
							 6*c13*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 18*c12*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
							 6*c13*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 18*c12*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
							 48*c13*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 12*c13*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
							 54*c12*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 18*c12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
							 48*c13*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 18*c13*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
							 36*c12*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
							 6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
							 6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 
							 6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) - 
							 18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
							 18*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 18*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
							 18*c12*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 18*c12*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 
							 60*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 6*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
							 54*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 48*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
							 72*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c23,4)*c33*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 
							 12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
							 42*cpow(c23,3)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 42*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
							 96*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
							 72*c23*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
							 6*c11*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
							 18*c12*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 30*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
							 3*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 3*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
							 24*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
							 36*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 6*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,3)*cpow(p1,2)*cpow(p2,4) - 
							 2*cpow(c23,6)*cpow(c55,3)*cpow(p2,6) + 6*c22*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,6) - 6*cpow(c22,2)*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,6) + 
							 2*cpow(c22,3)*cpow(c33,3)*cpow(c55,3)*cpow(p2,6) - 12*cpow(c23,5)*c44*cpow(c55,3)*cpow(p2,6) + 24*c22*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,6) - 
							 12*cpow(c22,2)*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,6) - 24*cpow(c23,4)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) + 
							 33*c22*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 9*cpow(c22,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 
							 16*cpow(c23,3)*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) + 18*c22*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) - 
							 3*cpow(c23,4)*c33*c44*cpow(c55,2)*c66*cpow(p2,6) + 6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,6) - 
							 3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) + 
							 12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 
							 18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) - 
							 3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,6) + 
							 2*cpow(c33,3)*cpow(c44,3)*cpow(c66,3)*cpow(p2,6) + csqrt(cpow(-2*cpow(c33,3)*cpow(c44,3) + 3*cpow(c33,3)*cpow(c44,2)*c55 + 
														       3*cpow(c33,2)*cpow(c44,3)*c55 + 3*cpow(c33,3)*c44*cpow(c55,2) - 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2) + 3*c33*cpow(c44,3)*cpow(c55,2) - 
														       2*cpow(c33,3)*cpow(c55,3) + 3*cpow(c33,2)*c44*cpow(c55,3) + 3*c33*cpow(c44,2)*cpow(c55,3) - 2*cpow(c44,3)*cpow(c55,3) - 
														       6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,2) + 6*c11*cpow(c33,3)*cpow(c44,3)*cpow(p1,2) + 6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2) - 
														       6*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2) - 3*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,2) - 6*c11*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) - 
														       12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2) + 3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2) - 3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2) + 
														       6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) + 
														       12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2) - 6*cpow(c13,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 
														       3*c11*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 6*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) - 9*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2) + 
														       6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2) + 12*c13*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) + 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2) - 
														       12*c13*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2) - 3*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2) - 
														       6*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2) + 12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2) + 6*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2) - 
														       6*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2) - 3*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2) - 6*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,4) + 
														       12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,4) - 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*cpow(p1,4) + 
														       3*cpow(c13,4)*c33*cpow(c44,2)*c55*cpow(p1,4) - 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4) + 
														       3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4) - 6*cpow(c13,4)*cpow(c44,3)*c55*cpow(p1,4) + 3*c11*cpow(c13,2)*c33*cpow(c44,3)*c55*cpow(p1,4) - 
														       24*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4) + 3*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 
														       24*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4) + 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 
														       12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4) - 24*cpow(c13,3)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
														       6*c11*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) - 33*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 
														       18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4) + 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 
														       18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4) - 24*cpow(c13,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) + 
														       9*c11*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4) - 
														       6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4) + 6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4) - 
														       6*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4) + 6*c11*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4) - 
														       6*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 
														       12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4) - 12*c13*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4) - 
														       12*c13*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) - 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4) + 
														       3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4) - 6*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) + 
														       3*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4) - 2*cpow(c13,6)*cpow(c44,3)*cpow(p1,6) + 6*c11*cpow(c13,4)*c33*cpow(c44,3)*cpow(p1,6) - 
														       6*cpow(c11,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(p1,6) + 2*cpow(c11,3)*cpow(c33,3)*cpow(c44,3)*cpow(p1,6) - 
														       12*cpow(c13,5)*cpow(c44,3)*c55*cpow(p1,6) + 24*c11*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,6) - 
														       12*cpow(c11,2)*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,6) - 24*cpow(c13,4)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) + 
														       33*c11*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 9*cpow(c11,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,6) - 
														       16*cpow(c13,3)*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) + 18*c11*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,6) - 
														       3*cpow(c13,4)*c33*cpow(c44,2)*c55*c66*cpow(p1,6) + 6*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,6) - 
														       3*cpow(c11,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,6) - 12*cpow(c13,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) + 
														       12*c11*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,6) - 12*cpow(c13,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 
														       18*c11*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,6) + 3*cpow(c13,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) - 
														       3*c11*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,6) + 6*c13*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,6) + 
														       2*cpow(c33,3)*cpow(c55,3)*cpow(c66,3)*cpow(p1,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p2,2) - 
														       3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p2,2) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p2,2) + 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,2) - 
														       6*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,2) + 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 
														       12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,2) + 
														       12*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) + 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,2) - 
														       6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(p2,2) - 3*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,2) - 
														       6*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,2) - 6*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 
														       3*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 6*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 9*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,2) - 
														       12*c23*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) - 9*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,2) + 6*cpow(c33,3)*cpow(c44,3)*c66*cpow(p2,2) - 
														       6*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,2) - 6*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,2) - 3*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,2) + 
														       12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,2) - 3*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,2) - 
														       3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 
														       18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
														       9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) + 6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 
														       12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,2)*cpow(p2,2) - 3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
														       6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 6*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
														       18*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
														       18*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 24*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 
														       18*c12*c13*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 9*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       12*c11*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) + 18*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
														       18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 24*c13*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
														       12*c13*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
														       12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 9*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
														       18*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
														       48*c13*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 18*c12*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 
														       18*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) - 18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,2) + 
														       12*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 12*c11*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
														       6*c11*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) - 18*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
														       3*cpow(c13,2)*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 6*c11*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 
														       42*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,2) + 18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) - 18*c12*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
														       6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 18*c13*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
														       36*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 6*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
														       6*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 24*c13*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 
														       18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,2) + 12*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 
														       12*c22*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 3*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
														       6*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 42*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
														       24*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 18*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) + 
														       36*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,2) - 3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 
														       3*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) + 6*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,2) - 
														       6*cpow(c13,4)*cpow(c23,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 3*cpow(c13,4)*c22*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
														       18*c12*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 3*c11*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
														       9*cpow(c12,2)*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
														       6*c11*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 
														       3*cpow(c11,2)*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) + 9*c11*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 
														       3*cpow(c11,2)*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,4)*c23*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 
														       18*c12*cpow(c13,3)*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 6*c11*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 
														       18*c11*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) + 6*cpow(c11,2)*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(p1,4)*cpow(p2,2) - 
														       24*cpow(c13,3)*cpow(c23,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,3)*c22*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
														       54*c12*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 6*c11*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
														       18*cpow(c12,2)*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
														       12*c11*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
														       48*cpow(c13,3)*c23*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 54*c12*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) + 
														       12*c11*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 18*c11*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,4)*cpow(p2,2) - 
														       24*cpow(c13,2)*cpow(c23,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 12*cpow(c13,2)*c22*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
														       36*c12*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 9*c11*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
														       27*cpow(c12,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 18*c11*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) - 
														       48*cpow(c13,2)*c23*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 36*c12*c13*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 
														       18*c11*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,4)*cpow(p2,2) + 6*cpow(c13,4)*c33*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) - 
														       12*c11*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 6*cpow(c11,2)*cpow(c33,3)*cpow(c44,3)*c66*cpow(p1,4)*cpow(p2,2) + 
														       18*cpow(c13,3)*c23*c33*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) - 
														       18*c11*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 18*c11*c12*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
														       42*cpow(c13,3)*c33*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) - 42*c11*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,4)*cpow(p2,2) + 
														       3*cpow(c13,2)*cpow(c23,2)*c33*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*cpow(c13,2)*c22*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
														       18*c12*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c11*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
														       9*cpow(c12,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 6*c11*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 
														       60*cpow(c13,2)*c23*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 54*c12*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
														       6*c11*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 96*cpow(c13,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) - 
														       18*c11*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,4)*cpow(p2,2) + 6*c13*cpow(c23,2)*c33*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
														       12*c13*c22*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 18*c12*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
														       48*c13*c23*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 72*c12*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) + 
														       72*c13*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,4)*cpow(p2,2) - 3*cpow(c13,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
														       3*c11*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 18*c13*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
														       18*c12*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
														       24*c13*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
														       6*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 6*c22*cpow(c33,3)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 
														       30*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) - 36*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(c66,2)*cpow(p1,4)*cpow(p2,2) + 
														       6*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,3)*cpow(p1,4)*cpow(p2,2) + 3*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p2,4) - 
														       6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p2,4) + 3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p2,4) + 
														       12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) - 12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p2,4) + 
														       12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p2,4) - 
														       6*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,4) + 12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,4) - 
														       6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*cpow(p2,4) - 6*cpow(c23,4)*c44*cpow(c55,3)*cpow(p2,4) + 3*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p2,4) - 
														       24*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,4) + 3*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) + 
														       24*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,4) - 24*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 
														       6*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 33*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) + 
														       18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,4) - 24*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) + 
														       9*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 18*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,4) - 
														       6*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p2,4) + 6*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p2,4) - 
														       12*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p2,4) - 6*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,4) + 
														       6*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,4) - 6*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
														       12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,4) - 
														       12*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 18*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,4) - 
														       6*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p2,4) + 3*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,4) + 
														       3*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,4) - 6*cpow(c13,2)*cpow(c23,4)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
														       3*cpow(c13,2)*c22*cpow(c23,2)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 18*c12*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       3*c11*cpow(c23,4)*c33*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 3*cpow(c13,2)*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       18*c12*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 9*cpow(c12,2)*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
														       6*c11*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 9*cpow(c12,2)*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       3*c11*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 24*cpow(c13,2)*cpow(c23,3)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
														       6*cpow(c13,2)*c22*c23*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 54*c12*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       12*c11*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       18*cpow(c12,2)*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
														       12*c11*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       24*cpow(c13,2)*cpow(c23,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 9*cpow(c13,2)*c22*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 
														       36*c12*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 12*c11*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       27*cpow(c12,2)*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) + 18*c11*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(p1,2)*cpow(p2,4) - 
														       12*c13*cpow(c23,4)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 6*c13*c22*cpow(c23,2)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
														       18*c12*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 6*c13*cpow(c22,2)*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
														       18*c12*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 48*c13*cpow(c23,3)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
														       12*c13*c22*c23*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 54*c12*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 
														       18*c12*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) - 48*c13*cpow(c23,2)*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
														       18*c13*c22*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 36*c12*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p1,2)*cpow(p2,4) + 
														       3*cpow(c13,2)*cpow(c23,2)*c33*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 6*cpow(c13,2)*c22*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 
														       18*c12*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 6*c11*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
														       9*cpow(c12,2)*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) - 6*c11*c22*cpow(c33,3)*cpow(c44,2)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
														       6*cpow(c13,2)*c23*c33*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) - 18*c12*c13*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 
														       12*c11*c23*cpow(c33,2)*cpow(c44,3)*c55*c66*cpow(p1,2)*cpow(p2,4) + 18*c13*cpow(c23,3)*c33*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
														       18*c13*c22*c23*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 18*c12*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 
														       18*c12*c22*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 60*c13*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 
														       6*c13*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 54*c12*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 
														       48*c13*c23*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) - 72*c12*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p1,2)*cpow(p2,4) + 
														       6*cpow(c23,4)*c33*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 12*c22*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 
														       6*cpow(c22,2)*cpow(c33,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 42*cpow(c23,3)*c33*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 
														       42*c22*c23*cpow(c33,2)*c44*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 96*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 
														       18*c22*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) + 72*c23*c33*cpow(c44,3)*cpow(c55,3)*c66*cpow(p1,2)*cpow(p2,4) - 
														       6*cpow(c13,2)*cpow(c33,2)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 6*c11*cpow(c33,3)*cpow(c44,3)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
														       18*c13*c23*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 18*c12*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
														       30*c13*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 3*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 
														       3*c22*cpow(c33,3)*c44*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 24*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) - 
														       36*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*cpow(c66,2)*cpow(p1,2)*cpow(p2,4) + 6*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,3)*cpow(p1,2)*cpow(p2,4) - 
														       2*cpow(c23,6)*cpow(c55,3)*cpow(p2,6) + 6*c22*cpow(c23,4)*c33*cpow(c55,3)*cpow(p2,6) - 6*cpow(c22,2)*cpow(c23,2)*cpow(c33,2)*cpow(c55,3)*cpow(p2,6) + 
														       2*cpow(c22,3)*cpow(c33,3)*cpow(c55,3)*cpow(p2,6) - 12*cpow(c23,5)*c44*cpow(c55,3)*cpow(p2,6) + 24*c22*cpow(c23,3)*c33*c44*cpow(c55,3)*cpow(p2,6) - 
														       12*cpow(c22,2)*c23*cpow(c33,2)*c44*cpow(c55,3)*cpow(p2,6) - 24*cpow(c23,4)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) + 
														       33*c22*cpow(c23,2)*c33*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 9*cpow(c22,2)*cpow(c33,2)*cpow(c44,2)*cpow(c55,3)*cpow(p2,6) - 
														       16*cpow(c23,3)*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) + 18*c22*c23*c33*cpow(c44,3)*cpow(c55,3)*cpow(p2,6) - 
														       3*cpow(c23,4)*c33*c44*cpow(c55,2)*c66*cpow(p2,6) + 6*c22*cpow(c23,2)*cpow(c33,2)*c44*cpow(c55,2)*c66*cpow(p2,6) - 
														       3*cpow(c22,2)*cpow(c33,3)*c44*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,3)*c33*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) + 
														       12*c22*c23*cpow(c33,2)*cpow(c44,2)*cpow(c55,2)*c66*cpow(p2,6) - 12*cpow(c23,2)*c33*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 
														       18*c22*cpow(c33,2)*cpow(c44,3)*cpow(c55,2)*c66*cpow(p2,6) + 3*cpow(c23,2)*cpow(c33,2)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) - 
														       3*c22*cpow(c33,3)*cpow(c44,2)*c55*cpow(c66,2)*cpow(p2,6) + 6*c23*cpow(c33,2)*cpow(c44,3)*c55*cpow(c66,2)*cpow(p2,6) + 
														       2*cpow(c33,3)*cpow(c44,3)*cpow(c66,3)*cpow(p2,6),2) + 4*
														  cpow(-cpow(c33*c44 + c33*c55 + c44*c55 + cpow(c13,2)*c44*cpow(p1,2) - c11*c33*c44*cpow(p1,2) + 2*c13*c44*c55*cpow(p1,2) - c33*c55*c66*cpow(p1,2) + 
															     cpow(c23,2)*c55*cpow(p2,2) - c22*c33*c55*cpow(p2,2) + 2*c23*c44*c55*cpow(p2,2) - c33*c44*c66*cpow(p2,2),2) - 
														       3*c33*c44*c55*(-c33 - c44 - c55 - cpow(c13,2)*cpow(p1,2) + c11*c33*cpow(p1,2) + c11*c44*cpow(p1,2) - 2*c13*c55*cpow(p1,2) + c44*c55*cpow(p1,2) + 
																      c33*c66*cpow(p1,2) + c55*c66*cpow(p1,2) - c11*c44*c55*cpow(p1,4) + cpow(c13,2)*c66*cpow(p1,4) - c11*c33*c66*cpow(p1,4) + 2*c13*c55*c66*cpow(p1,4) - 
																      cpow(c23,2)*cpow(p2,2) + c22*c33*cpow(p2,2) - 2*c23*c44*cpow(p2,2) + c22*c55*cpow(p2,2) + c44*c55*cpow(p2,2) + c33*c66*cpow(p2,2) + 
																      c44*c66*cpow(p2,2) + cpow(c13,2)*c22*cpow(p1,2)*cpow(p2,2) - 2*c12*c13*c23*cpow(p1,2)*cpow(p2,2) + c11*cpow(c23,2)*cpow(p1,2)*cpow(p2,2) + 
																      cpow(c12,2)*c33*cpow(p1,2)*cpow(p2,2) - c11*c22*c33*cpow(p1,2)*cpow(p2,2) - 2*c12*c13*c44*cpow(p1,2)*cpow(p2,2) + 
																      2*c11*c23*c44*cpow(p1,2)*cpow(p2,2) + 2*c13*c22*c55*cpow(p1,2)*cpow(p2,2) - 2*c12*c23*c55*cpow(p1,2)*cpow(p2,2) - 
																      2*c12*c44*c55*cpow(p1,2)*cpow(p2,2) - 2*c13*c23*c66*cpow(p1,2)*cpow(p2,2) + 2*c12*c33*c66*cpow(p1,2)*cpow(p2,2) - 
																      2*c13*c44*c66*cpow(p1,2)*cpow(p2,2) - 2*c23*c55*c66*cpow(p1,2)*cpow(p2,2) - 4*c44*c55*c66*cpow(p1,2)*cpow(p2,2) - c22*c44*c55*cpow(p2,4) + 
																      cpow(c23,2)*c66*cpow(p2,4) - c22*c33*c66*cpow(p2,4) + 2*c23*c44*c66*cpow(p2,4)),3)),0.3333333333333333)/(3.*cpow(2,0.3333333333333333)*c33*c44*c55);

/*dy ###################################################################################################################*/

    dz = p3*((c13 + c55)*cpow(p1,2)*((c23 + c44)*(c12 + c66)*cpow(p2,2) - (c13 + c55)*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))) + 
	     (c13 + c55)*cpow(p1,2)*((c23 + c44)*(c12 + c66)*cpow(p2,2) - 2*c44*(c13 + c55)*cpow(p3,2) - 
				     (c13 + c55)*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))) - 
	     (c23 + c44)*cpow(p2,2)*(-((c13 + c55)*(c12 + c66)*cpow(p1,2)) + (c23 + c44)*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) - 
	     (c23 + c44)*cpow(p2,2)*(-((c13 + c55)*(c12 + c66)*cpow(p1,2)) + 2*(c23 + c44)*c55*cpow(p3,2) + 
				     (c23 + c44)*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) + 
	     2*c33*(-(cpow(c12 + c66,2)*cpow(p1,2)*cpow(p2,2)) + (-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))*
		    (-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) + 
	     2*(-1 + c55*cpow(p1,2) + c44*cpow(p2,2) + c33*cpow(p3,2))*(c55*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2)) + 
									c44*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + 2*c55*cpow(p3,2))));

/*dt ###################################################################################################################*/

    dt = cpow(p1,2)*(2*(c23 + c44)*(-(c11*(c23 + c44)) + c12*(c13 + c55) + (c13 + c55)*c66)*cpow(p2,2)*cpow(p3,2) + 
		     (c13 + c55)*cpow(p3,2)*((c23 + c44)*(c12 + c66)*cpow(p2,2) - (c13 + c55)*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))) + 
		     (c13 + c55)*cpow(p3,2)*(-2*(c13 + c55)*c66*cpow(p1,2) + (c23 + c44)*(c12 + c66)*cpow(p2,2) - 
					     (c13 + c55)*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))) + 
		     2*(-1 + c55*cpow(p1,2) + c44*cpow(p2,2) + c33*cpow(p3,2))*(-(cpow(c12,2)*cpow(p2,2)) + c11*(-1 + 2*c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2)) + 
										c66*(-1 - 2*c12*cpow(p2,2) + c55*cpow(p3,2))) + 2*c55*(-(cpow(c12 + c66,2)*cpow(p1,2)*cpow(p2,2)) + 
																       (-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2)))) + 
	cpow(p2,2)*(-2*(c13 + c55)*(c13*c22 - c12*(c23 + c44) + c22*c55 - c23*c66 - c44*c66)*cpow(p1,2)*cpow(p3,2) - 
		    (c23 + c44)*cpow(p3,2)*(-((c13 + c55)*(c12 + c66)*cpow(p1,2)) + (c23 + c44)*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) - 
		    (c23 + c44)*cpow(p3,2)*(-((c13 + c55)*(c12 + c66)*cpow(p1,2)) + 2*(c23 + c44)*c66*cpow(p2,2) + 
					    (c23 + c44)*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) + 
		    2*c44*(-(cpow(c12 + c66,2)*cpow(p1,2)*cpow(p2,2)) + (-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))*
			   (-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) + 
		    2*(-1 + c55*cpow(p1,2) + c44*cpow(p2,2) + c33*cpow(p3,2))*(-(cpow(c12,2)*cpow(p1,2)) + c66*(-1 - 2*c12*cpow(p1,2) + c44*cpow(p3,2)) + 
									       c22*(-1 + c11*cpow(p1,2) + 2*c66*cpow(p2,2) + c55*cpow(p3,2)))) + 
	cpow(p3,2)*((c13 + c55)*cpow(p1,2)*((c23 + c44)*(c12 + c66)*cpow(p2,2) - (c13 + c55)*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))) + 
		    (c13 + c55)*cpow(p1,2)*((c23 + c44)*(c12 + c66)*cpow(p2,2) - 2*c44*(c13 + c55)*cpow(p3,2) - 
					    (c13 + c55)*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))) - 
		    (c23 + c44)*cpow(p2,2)*(-((c13 + c55)*(c12 + c66)*cpow(p1,2)) + (c23 + c44)*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) - 
		    (c23 + c44)*cpow(p2,2)*(-((c13 + c55)*(c12 + c66)*cpow(p1,2)) + 2*(c23 + c44)*c55*cpow(p3,2) + 
					    (c23 + c44)*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) + 
		    2*c33*(-(cpow(c12 + c66,2)*cpow(p1,2)*cpow(p2,2)) + (-1 + c66*cpow(p1,2) + c22*cpow(p2,2) + c44*cpow(p3,2))*
			   (-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + c55*cpow(p3,2))) + 
		    2*(-1 + c55*cpow(p1,2) + c44*cpow(p2,2) + c33*cpow(p3,2))*(c55*(-1 + c66*cpow(p1,2) + c22*cpow(p2,2)) + 
									       c44*(-1 + c11*cpow(p1,2) + c66*cpow(p2,2) + 2*c55*cpow(p3,2))));
         
         
    return creal(dz)/creal(dt);

}


