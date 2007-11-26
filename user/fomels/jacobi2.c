/* Jacobi-like method for general complex matrix */
/*
  Copyright (C) 2007 University of Texas at Austin
  
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
#include <float.h>

#include <rsf.h>

static sf_complex *ak;

void jacobi2_init(int n /* matrix size */)
/*< initialize >*/
{
    ak = sf_complexalloc(n);
}

void jacobi2_close(void)
/*< free allocated storage >*/
{
    free(ak);
}

float jacobi2(sf_complex** a /* matrix to rotate */, 
	      int n, int k, int l   /* indeces */) 
/*< Jacobi-like transformation >*/
/* The algorithm is modified from 

   J.P. Killingbeck, A.Grosjean, and G. Jolicard, 2004, A simple
   method for complex eigenvalues: Journal of Physics A: Mathematical
   and General, v. 37, L567-L572. */
{
    int i;
    sf_complex t, s, akl;

    if (k==l || cabsf(a[k][l]) < FLT_EPSILON) return 0.0f;
    

    t = 0.5*(a[l][l]-a[k][k]);
    s = csqrtf(t*t+a[l][k]*a[k][l]);

    if (cabsf(t) < FLT_EPSILON && cabsf(s) < FLT_EPSILON) return 0.0f;

    if (cabsf(t+s) > cabsf(t-s)) {
	t = a[k][l]/(t+s);
    } else {
	t = a[k][l]/(t-s);
    }

    akl = a[k][l] + t*(a[k][k]-a[l][l])-t*t*a[l][k];
    for (i=0; i < n; i++) {
	ak[i] = a[i][l] + t*a[i][k];
    }
    for (i=0; i < n; i++) {
	a[k][i] -= t*a[l][i];
    }
    for (i=0; i < n; i++) {
	a[i][l] = ak[i];
    }
    a[k][l] = akl;

    return cabsf(t);
}
