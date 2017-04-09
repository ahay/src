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

#ifndef _ml_traveltime_h

typedef float (*func1)(int, float);
/* function pointer for int and float -> float */
/*^*/

#endif


float y_1k, y_k1,v_k1,v_k;
int j;

float traveltime_k(int j,float v_k1,float y_k,float y_k1,func1 z, func1 zd1, func1 zd2)
/*<Traveltime from layer k_th to k+1_th>*/
{
	float t_k;
	
	t_k = hypotf(y_k-y_k1,z(j,y_k)-z(j+1,y_k1))/v_k1;
	
	return t_k;

}


float traveltime_k_k(int j, float v_k1,float y_1k,float y_k,float y_k1,func1 z, func1 zd1, func1 zd2)
/*<The derivative with respect to y_k of the traveltime from layer k_th to k+1_th >*/
{

	float t_k_k;
	
	t_k_k = (z(j,y_k)*zd1(j,y_k)+y_k-z(j+1,y_k1)*zd1(j,y_k)-y_k1)/(v_k1*hypotf(z(j+1,y_k1)-z(j,y_k), y_k1-y_k));
	
	return t_k_k;
}

float traveltime_1k_k(int j, float v_k,float y_1k,float y_k,float y_k1,func1 z, func1 zd1, func1 zd2)
/*<The derivative with respect to y_k of the traveltime from layer k-1_th to k_th >*/
{
	
	float t_1k_k;
	
	t_1k_k = (z(j,y_k)*zd1(j,y_k)+y_k-z(j-1,y_1k)*zd1(j,y_k)-y_1k)/(v_k*hypotf(z(j,y_k)-z(j-1,y_1k), y_k-y_1k));
	
	return t_1k_k;
}

float traveltime_1k_k_k(int j,float v_k,float y_1k,float y_k,float y_k1,func1 z, func1 zd1, func1 zd2)
/*<The second derivative with respect to y_k of the traveltime from layer k-1_th to k_th >*/
{
	
	float t_1k_k_k;
	
	t_1k_k_k = (4*pow((y_k-y_1k)*zd1(j,y_k)+z(j-1,y_1k)-z(j,y_k),2)-4*(pow(hypotf(z(j-1,y_1k)-z(j,y_k), y_1k-y_k),2))*(z(j-1,y_1k)-z(j,y_k))*zd2(j,y_k))/(4*v_k*pow(hypotf(z(j,y_k)-z(j-1,y_1k), y_k-y_1k),3));
	
	return t_1k_k_k;
}

float traveltime_k_k_k(int j,float v_k1,float y_1k,float y_k,float y_k1,func1 z, func1 zd1, func1 zd2)
/*<The second derivative with respect to y_k of the traveltime from layer k_th to k+1_th >*/
{
	
	float t_k_k_k;
	
	t_k_k_k = ((4*(pow(hypotf(z(j+1,y_k1)-z(j,y_k), y_k1-y_k),2))*((z(j,y_k)-z(j+1,y_k1))*zd2(j,y_k)+pow(zd1(j,y_k), 2)+1))-(4*pow(z(j+1,y_k1)*zd1(j,y_k)+y_k1-z(j,y_k)*zd1(j,y_k)-y_k, 2)))/(4*v_k1*pow(hypotf(z(j+1,y_k1)-z(j,y_k), y_k1-y_k), 3));
	
	return t_k_k_k;
}

float traveltime_1k_1k_k(int j,float v_k,float y_1k,float y_k,float y_k1,func1 z, func1 zd1, func1 zd2)
/*<The second derivative with respect to y_k-1 and y_k of the traveltime from layer k-1_th to k_th >*/
{
	
	float t_1k_1k_k;
	
	t_1k_1k_k = (-1)*(((y_1k-y_k)*zd1(j,y_k)+z(j,y_k)-z(j-1,y_1k))*((y_1k-y_k)*zd1(j-1,y_1k)+z(j,y_k)-z(j-1,y_1k)))/(v_k*pow(hypotf(z(j,y_k)-z(j-1,y_1k), y_k-y_1k),3));
	
	return t_1k_1k_k;
}

float traveltime_k_k_k1(int j,float v_k1,float y_1k,float y_k,float y_k1,func1 z, func1 zd1, func1 zd2)
/*<The second derivative with respect to y_k and y_k+1 of the traveltime from layer k_th to k+1_th >*/
{
	
	float t_k_k_k1;
	
	t_k_k_k1 = (-1)*(((y_k-y_k1)*zd1(j+1,y_k1)+z(j+1,y_k1)-z(j,y_k))*((y_k-y_k1)*zd1(j,y_k)+z(j+1,y_k1)-z(j,y_k)))/(v_k1*pow(hypotf(z(j+1,y_k1)-z(j,y_k), y_k1-y_k),3));
	
	return t_k_k_k1;
}
