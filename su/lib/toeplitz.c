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

void stoepf (int n, float r[], float g[], float f[], float a[])
/*<
****************************************************************************
Solve a symmetric Toeplitz linear system of equations Rf=g for f
(float version) 
******************************************************************************
Input:
n		dimension of system
r		array[n] of top row of Toeplitz matrix
g		array[n] of right-hand-side column vector

Output:
f		array[n] of solution (left-hand-side) column vector
a		array[n] of solution to Ra=v (Claerbout, FGDP, p. 57)
******************************************************************************
Notes:
This routine does NOT solve the case when the main diagonal is zero, it
just silently returns.

The left column of the Toeplitz matrix is assumed to be equal to the top
row (as specified in r); i.e., the Toeplitz matrix is assumed symmetric.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
****************************************************************************>*/
{
    int i,j;
    float v,e,c,w,bot;

    if (r[0] == 0.0) return;

    a[0] = 1.0;
    v = r[0];
    f[0] = g[0]/r[0];

    for (j=1; j<n; j++) {
		
	/* solve Ra=v as in Claerbout, FGDP, p. 57 */
	a[j] = 0.0;
	f[j] = 0.0;
	for (i=0,e=0.0; i<j; i++)
	    e += a[i]*r[j-i];
	c = e/v;
	v -= c*e;
	for (i=0; i<=j/2; i++) {
	    bot = a[j-i]-c*a[i];
	    a[i] -= c*a[j-i];
	    a[j-i] = bot;
	}

	/* use a and v above to get f[i], i = 0,1,2,...,j */
	for (i=0,w=0.0; i<j; i++)
	    w += f[i]*r[j-i];
	c = (w-g[j])/v;
	for (i=0; i<=j; i++)
	    f[i] -= c*a[j-i];
    }
}


