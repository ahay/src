/* Tridiagonal matrix solver */
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

#include "alloc.h"

#include "tridiagonal.h"

#ifndef _sf_tridiagonal_h

#include "_bool.h"
/*^*/

typedef struct sf_Tris *sf_tris;
/* abstract data type */
/*^*/

#endif

struct sf_Tris {
    int n;
    float *d[2], *o[2], *x[2];  
};

sf_tris sf_tridiagonal_init (int n /* matrix size */)
/*< initialize >*/
{
    sf_tris slv;
    
    slv = (sf_tris) sf_alloc (1, sizeof(*slv));
    
    slv->n = n;
    slv->d[0] = sf_floatalloc (n);
    slv->d[1] = sf_floatalloc (n);

    slv->o[0] = sf_floatalloc (n-1);
    slv->o[1] = sf_floatalloc (n-1);

    slv->x[0] = sf_floatalloc (n);
    slv->x[1] = sf_floatalloc (n);

    return slv;
}

void sf_tridiagonal_define (sf_tris slv /* solver object */, 
			    float* diag /* diagonal */, 
			    float* offd /* off-diagonal */)
/*< fill the matrix >*/
{
    int k;
    float t;

    slv->d[0][0] = diag[0];
    for (k = 1; k < slv->n; k++) {
	t = offd[k-1]; 
	slv->o[0][k-1] = t / slv->d[0][k-1];
	slv->d[0][k] = diag[k] - t * slv->o[0][k-1];
    }
    slv->d[1][slv->n-1] = diag[slv->n-1];
    for (k = slv->n-2; k >= 0; k--) {
	t = offd[k]; 
	slv->o[1][k] = t / slv->d[1][k+1];
	slv->d[1][k] = diag[k] - t * slv->o[1][k];
    }
}

void sf_tridiagonal_const_define (sf_tris slv /* solver object */, 
				  float diag  /* diagonal */, 
				  float offd  /* off-diagonal */,
				  bool damp   /* damping */)
/*< fill the matrix for the Toeplitz case >*/
{
    int k;
    
    slv->d[0][0] = damp? diag+offd: diag;
    for (k = 1; k < slv->n; k++) {
	slv->o[0][k-1] = offd / slv->d[0][k-1];
	slv->d[0][k] = diag - offd * slv->o[0][k-1];
    }
    if (damp) slv->d[0][slv->n-1] += offd;
    slv->d[1][slv->n-1] = damp? diag+offd: diag;
    for (k = slv->n-2; k >= 0; k--) {
	slv->o[1][k] = offd / slv->d[1][k+1];
	slv->d[1][k] = diag - offd * slv->o[1][k];
    }
    if (damp) slv->d[1][0] += offd;
}

void sf_tridiagonal_solve (sf_tris slv /* solver object */, 
			   float* b /* in - right-hand side, out - solution */)
/*< invert the matrix >*/
{
    int k;

    slv->x[0][0] = b[0];
    for (k = 1; k < slv->n; k++) {
	slv->x[0][k] = b[k] - slv->o[0][k-1] * slv->x[0][k-1];
    }
    slv->x[1][slv->n-1] = b[slv->n-1];
    for (k = slv->n-2; k >= 0; k--) {
	slv->x[1][k] = b[k] - slv->o[1][k] * slv->x[1][k+1];
    }
    b[slv->n-1] = slv->x[0][slv->n-1] / slv->d[0][slv->n-1];
    for (k = slv->n-2; k >= slv->n/2; k--) {
	b[k] = slv->x[0][k] / slv->d[0][k] - slv->o[0][k] * b[k+1];
    }
    b[0] = slv->x[1][0] / slv->d[1][0];
    for (k = 1; k < slv->n/2; k++) {
	b[k] = slv->x[1][k] / slv->d[1][k] - slv->o[1][k-1] * b[k-1];
    }
}

void sf_tridiagonal_close (sf_tris slv)
/*< free allocated storage >*/
{
    free (slv->d[0]); free (slv->d[1]); 
    free (slv->o[0]); free (slv->o[1]); 
    free (slv->x[0]); free (slv->x[1]); 
    free (slv);
}

/* 	$Id: tridiagonal.c 7107 2011-04-10 02:04:14Z ivlad $	 */
