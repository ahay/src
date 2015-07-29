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

#include "quadratic.h"
#include "c99.h"
#include "_defs.h"

float sf_quadratic_solve (float a, float b, float c) 
/*< solves a x^2 + 2 b x + c == 0 for smallest positive x >*/
{
    float d;
    
    if (fabsf(c) < SF_EPS && 
	((b >   SF_EPS && a < -SF_EPS) || 
	 (b < - SF_EPS && a >  SF_EPS))) 
	return (-2.*b/a);
    
    d = b*b - a*c;
    if (d < 0.) return SF_HUGE;
    
    d = sqrtf(d)+fabsf(b);
    
    if (b*c <= 0.) {
	if (d > SF_EPS && fabsf(c) > 0.) return (fabsf(c)/d);
    } else if (b*a <= 0.) {
	if (fabsf(a) > SF_EPS && d > 0.) return (d/fabsf(a));
    } 
    
    return SF_HUGE;
}

/* 	$Id: quadratic.c 7107 2011-04-10 02:04:14Z ivlad $	 */
