/* Real Toeplitz solver */
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

#include "toeplitz.h"

static double dprod (int j, const double* a, const double* b) 
/* dot product */
{
    int i;
    double c;

    c = 0.;
    for (i=1; i <= j; i++) {
	c += a[j-i]*b[i];
    }
    return c;
}

void toeplitz (int n           /* matrix size */,
	       const double *r /* [n] top row of the matrix */, 
	       double *f       /* [n] inverted in place */,
	       double *a       /* [n] work array */)
/*< apply the solver >*/
{    
    int i,j;
    double e,c,w, bot, v;
    
    a[0] = 1.;
    v=r[0];
    f[0] /= v;
    
    for (j=1; j < n; j++) {
	e = dprod(j,a,r);
	c = -e/v;

	v += c*e;
       
	for (i=1; i <= j/2; i++) {
	    bot  = a[j-i] + c*a[i];
	    a[i] += c*a[j-i];
	    a[j-i] = bot;
	}
	a[j] = c;
       
	w = dprod(j,f,r);
	c = (f[j]-w)/v;
       
	for (i=0; i < j; i++) {
	    f[i] += c*a[j-i];
	}
	f[j] = c;
    }
}

/* 	$Id$	 */

