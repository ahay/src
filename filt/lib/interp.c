/* Basic interpolation functions */
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

#include "interp.h"
#include "error.h"

#ifndef _sf_interp_h

typedef void (*sf_interpolator)(float,int,float*);
/* generic interpolation interface */
/*^*/

#endif

void sf_bin_int (float x, int n, float* w) 
/*< nearest neighbor >*/
{
  int i;

  w[0] = 1.;
  
  for (i = 1; i < n; i++)
      w[i] = 0.;
}

void sf_lin_int (float x, int n, float* w) 
/*< linear >*/
{
    int i;
    
    if (n == 1) {
	w[0] = 1.;
    } else {
	w[1] = x;
	w[0] = 1. - x;
	for (i = 2; i < n; i++)
	    w[i] = 0.;
    }
}

void sf_lg_int (float x, int n, float* w) 
/*< Lagrangian >*/
{
    int i, j, nc;
    float f, xi;

    nc = (n-1)*0.5;
    for (i=0; i < n; i++) {
	f = 1.;
	xi = x + nc - i;
	for (j=0; j < n; j++) {
	    if (i != j) f *= (1. + xi / (i - j));
	}
	w[i] = f;
    }
}

void sf_taylor (float x, int n, float* w) 
/*< Taylor >*/
{
    int i;
    float xi;

    xi = 1.;
    for (i=0; i < n; i++) {
	if (i > 0) xi *= (x + 1. -i)/i;
	w[i] = xi;
    }
}


/* 	$Id$	 */
