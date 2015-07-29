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
    sf_complex *d, *o;  
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

void ctridiagonal_define (ctris slv, sf_complex* diag, sf_complex* offd)
/*< Fill the matrix with diagonal and off-diagonal elements.
...
The matrix is symmetric but not necessarily self-adjoint.
 >*/
{
    int k;
    sf_complex t;

    slv->d[0] = diag[0];
    for (k = 1; k < slv->n; k++) {
	t = offd[k-1]; 
#ifdef SF_HAS_COMPLEX_H
	slv->o[k-1] = t / slv->d[k-1];
	slv->d[k] = diag[k] - t * slv->o[k-1];
#else
	slv->o[k-1] = sf_cdiv(t,slv->d[k-1]);
	slv->d[k] = sf_csub(diag[k],sf_cmul(t,slv->o[k-1]));
#endif
    }
}

void ctridiagonal_const_define (ctris slv, 
				sf_complex diag, sf_complex offd)
/*< Fill the matrix with constant diagonal and off-diagonal elements.
...
The matrix is symmetric but not necessarily self-adjoint.
>*/
{
    int k;
    
    slv->d[0] = diag;
    for (k = 1; k < slv->n; k++) {
#ifdef SF_HAS_COMPLEX_H
	slv->o[k-1] = offd / slv->d[k-1];
	slv->d[k] = diag - offd * slv->o[k-1];
#else
	slv->o[k-1] = sf_cdiv(offd,slv->d[k-1]);
	slv->d[k] = sf_csub(diag,sf_cmul(offd,slv->o[k-1]));
#endif
    }
}

void ctridiagonal_solve (ctris slv, sf_complex* b)
/*< Invert in-place. The right-hand side b is replaced with the solution. >*/ 
{
    int k, n;

    n = slv->n;

    for (k = 1; k < n; k++) {
#ifdef SF_HAS_COMPLEX_H
	b[k] -= slv->o[k-1] * b[k-1];
#else
	b[k] = sf_csub(b[k],sf_cmul(slv->o[k-1],b[k-1]));
#endif
    }
#ifdef SF_HAS_COMPLEX_H
    b[n-1] /= slv->d[n-1];
#else
    b[n-1] = sf_cdiv(b[n-1],slv->d[n-1]);
#endif
    for (k = n-2; k >= 0; k--) {
#ifdef SF_HAS_COMPLEX_H
	b[k] = b[k] / slv->d[k] - slv->o[k] * b[k+1];
#else
	b[k] = sf_csub(sf_cdiv(b[k],slv->d[k]),
		       sf_cmul(slv->o[k],b[k+1]));
#endif
    }
}

void ctridiagonal_close (ctris slv)
/*< Free allocated storage >*/
{
    free (slv->d); 
    free (slv->o); 
    free (slv);
}

/* 	$Id: ctridiagonal.c 7107 2011-04-10 02:04:14Z ivlad $	 */
