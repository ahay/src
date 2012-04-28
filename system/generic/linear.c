/* Linear 1-D interpolation */
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

#include <rsf.h> 

#include "linear.h"

static int n;
static float *x, *a, *b;

void linear_init(int n1 /* trace length */) 
/*< initialize >*/
{
    n = n1;
    x = sf_floatalloc(n);
    a = sf_floatalloc(n);
    b = sf_floatalloc(n);
}

void linear_close (void) 
/*< free allocated storage >*/
{
    free (x);
    free (a);
    free (b);
}

void linear_coeffs(float** table)
/*< feel coefficients table >*/
{
    int k;
    float xk, fk, xp, fp;
    
    for (k=n-1; k >= 0; k--) {
	x[k] = xk = table[k][0];
	a[k] = fk = table[k][1];
	
	if (k < n-1) {
	    b[k] = (fp-fk)/(xp-xk);
	} else {
	    b[k] = 0.;
	}
 
	xp = xk;
	fp = fk;
    }
}

float linear_eval(float y)
/*< evaluate a cubic spline >*/
{
    float dh=0., s;
    int k;

    /* find the interval for x */
    for (k=n-1; k >=0; k--) {
	dh = y - x[k];
	if (dh >= 0.) break;
    }
    if (k == n-1) {
	s = a[k];
    } else if (k < 0) {
	s = a[0];
    } else {
	s = a[k] + dh*b[k];
    }

    return s;
}
