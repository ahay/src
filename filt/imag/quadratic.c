/* Solving quadratic equations. */
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

#include <float.h>
/*^*/

#include <rsf.h>

#include "quadratic.h"

#ifndef _quadratic_h

#define EPS FLT_EPSILON
#define HUG FLT_MAX
/*^*/

#endif

float quadratic_solve (float a, float b, float c) 
/*< solves a x^2 + 2 b x + c == 0 for smallest positive x >*/
{
    float d;
    
    if (fabsf(c) < EPS && 
	((b > EPS && a  < -EPS) || (b < - EPS && a > EPS))) 
	return (-2.*b/a);
    
    d = b*b - a*c;
    if (d < 0.) return HUG;
    
    d = sqrtf(d)+fabsf(b);
    
    if (b*c <= 0.) {
	if (d > EPS && fabsf(c) > 0.) return (fabsf(c)/d);
    } else if (b*a <= 0.) {
	if (fabsf(a) > EPS && d > 0.) return (d/fabsf(a));
    } 
    
    return HUG;
}

/* 	$Id$	 */
