/* Complex Toeplitz solver */
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
/*^*/

#include "ctoeplitz.h"
/*^*/

static int n;
static float complex *a;

static float complex cdprod (int j, 
			     const float complex* a, const float complex* b) 
/* complex dot product */
{
    int i;
    float complex c;
    c = 0.;
    for (i=1; i <= j; i++) {
	c += a[j-i]*conjf(b[i]);
    }
    return c;
}

void ctoeplitz_init (int n_in /* matrix size */)
/*< initialize >*/
{
    n = n_in;
    a = sf_complexalloc (n);
    a[0] = 1.;
}

void ctoeplitz_solve (const float complex *r /* top row of the matrix */, 
		      float complex *f       /* inverted in place */)
/*< apply the solver >*/
{    
    int i,j;
    float complex e,c,w, bot;
    float v;
    
    v=crealf(r[0]);
    f[0] /= v;
    
    for (j=1; j < n; j++) {
	e = cdprod(j,a,r);
	c = -e/v;

	v += crealf(c)*crealf(e) + cimagf(c)*cimagf(e);
       
	for (i=1; i <= j/2; i++) {
	    bot  = a[j-i] + c*conjf(a[i]);
	    a[i] = a[i] + c*conjf(a[j-i]);
	    a[j-i] = bot;
	}
	a[j] = c;
       
	w = cdprod(j,f,r);
	c = (f[j]-w)/v;
       
	for (i=0; i < j; i++) {
	    f[i] += c*conjf(a[j-i]);
	}
	f[j] = c;
    }
}

void ctoeplitz_close()
/*< free allocated storage >*/
{
    free (a);
}

/* 	$Id$	 */

