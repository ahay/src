/* Similarity operator */

/*
  Copyright (C) 2012 University of Texas at Austin
  
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


float similor_gaussian(float a, float b, float *p)
/*< Gaussian similarity >*/
{
	float t;
	t = (a-b)/p[0];
	return expf(-t*t/2.0);
}

typedef float (*similor)(float a, float b, float *p);
/* similarity operator interface */
/*^*/

similor similor_c2f(char *c)
/*< operator selector >*/
{
	if(strcmp(c,"gaussian")==0)	return similor_gaussian;
	else return NULL;
}


