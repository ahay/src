/* Complex tridiagonal solver. */
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

#include <rsf.h>
/*^*/

#include "ctridiagonal.h"

#ifndef _ctridiagonal_h

typedef struct CTris *ctris;
/* abstract data type */
/*^*/

#endif

struct CTris {
    int n;
    float complex *d, *o;  
};

ctris ctridiagonal_init (int n)
/*< Initialize a solver object for n by n tridiagonal matrix >*/
{
    ctris slv;
    
    slv = (ctris) sf_alloc (1, sizeof(*slv));
    
    slv->n = n;
    slv->d = sf_complexalloc (n);
    slv->o = sf_complexalloc (n-1);

    return slv;
}

void ctridiagonal_define (ctris slv, float complex* diag, float complex* offd)
/*< Fill the matrix with diagonal and off-diagonal elements.
...
The matrix is symmetric but not necessarily self-adjoint.
 >*/
{
    int k;
    float complex t;

    slv->d[0] = diag[0];
    for (k = 1; k < slv->n; k++) {
	t = offd[k-1]; 
	slv->o[k-1] = t / slv->d[k-1];
	slv->d[k] = diag[k] - t * slv->o[k-1];
    }
}

void ctridiagonal_const_define (ctris slv, 
				float complex diag, float complex offd)
/*< Fill the matrix with constant diagonal and off-diagonal elements.
...
The matrix is symmetric but not necessarily self-adjoint.
>*/
{
    int k;
    
    slv->d[0] = diag;
    for (k = 1; k < slv->n; k++) {
	slv->o[k-1] = offd / slv->d[k-1];
	slv->d[k] = diag - offd * slv->o[k-1];
    }
}

void ctridiagonal_solve (ctris slv, float complex* b)
/*< Invert in-place. The right-hand side b is replaced with the solution. >*/ 
{
    int k, n;

    n = slv->n;

    for (k = 1; k < n; k++) {
	b[k] -= slv->o[k-1] * b[k-1];
    }
    b[n-1] /= slv->d[n-1];
    for (k = n-2; k >= 0; k--) {
	b[k] = b[k] / slv->d[k] - slv->o[k] * b[k+1];
    }
}

void ctridiagonal_close (ctris slv)
/*< Free allocated storage >*/
{
    free (slv->d); 
    free (slv->o); 
    free (slv);
}

/* 	$Id$	 */
