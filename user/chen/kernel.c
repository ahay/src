/* nonlinear kernel functions */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <rsf.h>


float kernel_gaussian(float t)
/*< normalized Gaussian function >*/
{
	return expf(-t*t/2.0);
}

float kernel_sigmoid(float t)
/*< sigmoid/logistic function >*/
{
	return 1.0 / ( 1 + expf(-t) );
}

float kernel_func4(float t)
/*< 1/(1+t^4) function >*/
{
	return 1.0 / ( 1 + powf(t, 4.0) );
}

float kernel_dsigmoid(float t)
/*< sigmoid/logistic differential >*/
{
	float a, b;
	a = expf(-t);
	b = 1+a;
	return a / (b*b);
}

typedef float (*kernel)(float a);
/* nonlinear kernel interface */
/*^*/

kernel kernel_c2f(char *c)
/*< operator selector >*/
{
	if(strcmp(c,"gaussian")==0)	return kernel_gaussian;
	else if(strcmp(c,"sigmoid")==0)	return kernel_sigmoid;
	else if(strcmp(c,"dsigmoid")==0)	return kernel_dsigmoid;
	else if(strcmp(c,"func4")==0)	return kernel_func4;
	else return NULL;
}


