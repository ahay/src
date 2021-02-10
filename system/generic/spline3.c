/* Cubic spline interpolation */
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

#include "spline3.h"

static int n;
static float *h, *a, *b, *c, *d, **coeff, *x, x0, dx, *fp;
static sf_tris slv;

void spline3_init(int n1     /* trace length */,
		 float *fp1 /* end-point derivatives */) 
/*< initialize >*/
{
    n = n1;
    h = sf_floatalloc(n-1);
    x = sf_floatalloc(n-1);
    a = sf_floatalloc(n);
    coeff = sf_floatalloc2(n,4);
    /* rename for convenience */
    b = coeff[1];
    c = coeff[2];
    d = coeff[3];

    fp = fp1;

    slv = sf_tridiagonal_init (NULL==fp? n-2: n);
}

void spline3_init1(int n1, float o1, float d1)
/*< initialize for regular trace interpolation >*/ 
{
    n = n1;
    x0 = o1;
    dx = d1;

    coeff = sf_floatalloc2(n-1,4);
    /* rename for convenience */
    b = coeff[1];
    c = coeff[2];
    d = coeff[3];
    
    fp = NULL;

    slv = sf_tridiagonal_init (n-2);
    sf_tridiagonal_const_define (slv,4.*d1,d1,false);
}

void spline3_close (void) 
/*< free allocated storage >*/
{
    free (h);
    free (a);
    free (*coeff);
    free (coeff);
    sf_tridiagonal_close (slv);
}

void spline3_close1 (void) 
/*< free allocated storage for regular trace >*/
{
    free (*coeff);
    free (coeff);
    sf_tridiagonal_close (slv);
}

void spline_coeffs(float** table)
/*< feel spline coefficients table >*/
{
    int k;
    float xk, fk;
    
    for (k=0; k < n-1; k++) {
	x[k] = xk = table[k][0];                /* table  knots */
	h[k] = table[k+1][0] - xk;       /* interval length */
	coeff[0][k] = fk = table[k][1];
	b[k] = (table[k+1][1]-fk)/h[k];  /* divided difference */
    }
    for (k=0; k < n-2; k++) {
	a[k+1] = 2.*(h[k+1] + h[k]);         /* diagonal */
	c[k+1] = b[k+1] - b[k];            /* right-hand side */
    }

    /* solve the tridiagonal system */
    if (NULL == fp) {
	c[0] = 0.;
	c[n-1] = 0.;
	sf_tridiagonal_define (slv,a+1,h+1);
	sf_tridiagonal_solve(slv,c+1);
    } else {
	a[0] = 2*h[0];
	a[n-1] = 2*h[n-2];
	c[0] = b[0] - fp[0];
	c[n-1] = fp[1] - b[n-2];
	sf_tridiagonal_define (slv,a,h);
	sf_tridiagonal_solve(slv,c);
    }

    for (k=0; k < n-1; k++) {
	d[k] = (c[k+1]-c[k])/h[k];
	b[k] -= (c[k+1]+2.*c[k])*h[k];
	c[k] *= 3.;
    }
}

void spline_coeffs1(float* table1)
/*< feel spline coefficients table for regular trace >*/
{
    int k;
    float fk;
    
    for (k=0; k < n-1; k++) {
	coeff[0][k] = fk = table1[k];
	b[k] = (table1[k+1]-fk)/dx;  /* divided difference */
    }
    for (k=0; k < n-2; k++) {
	c[k+1] = b[k+1] - b[k];            /* right-hand side */
    }
    c[0] = 0;

    /* solve the tridiagonal system */
    sf_tridiagonal_solve(slv,c+1);

    for (k=0; k < n-1; k++) {
	if (k < n-2) {
	    d[k] = (c[k+1]-c[k])/dx;
	    b[k] -= (c[k+1]+2.*c[k])*dx;
	} else {
	    d[k] = -c[k]/dx;
	    b[k] -= 2.*c[k]*dx;
	}
	c[k] *= 3.;
    }
}

float spline_eval(float y)
/*< evaluate a cubic spline >*/
{
    float dh=0., s;
    int i, k;

    /* find the interval for x */
    for (k=n-2; k >=0; k--) {
	dh = y - x[k];
	if (dh >= 0.) break;
    }
    if (k < 0) k = 0;
    
    /* evaluate cubic by Horner's rule */
    s = coeff[3][k];
    for (i=2; i >=0; i--) {
	s = s*dh + coeff[i][k];
    }
    return s;
}

float spline_eval1(float y)
/*< evaluate a cubic spline for regular trace >*/
{
    float dh=0., s;
    int i, k;

    k = (y-x0)/dx;
    if (k < 0)   k=0;
    if (k > n-2) k=n-2;
    dh = y - (x0+k*dx);
    
    /* evaluate cubic by Horner's rule */
    s = coeff[3][k];
    for (i=2; i >=0; i--) {
	s = s*dh + coeff[i][k];
    }
    return s;
}

/* 	$Id$	 */
