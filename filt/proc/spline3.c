#include <rsf.h> 

#include "spline3.h"
#include "tridiagonal.h"

static int n;
static float *h, *a, *b, *c, *d, **coeff, *x, x0, dx;
static tris slv;

void spine3_init(int n1) {
    n = n1;
    h = sf_floatalloc(n-1);
    x = sf_floatalloc(n-1);
    a = sf_floatalloc(n-2);
    coeff = sf_floatalloc2(n-1,4);
    /* rename for convenience */
    b = coeff[1];
    c = coeff[2];
    d = coeff[3];

    slv = tridiagonal_init (n-2);
}

void spine3_init1(int n1, float o1, float d1) {
    n = n1;
    x0 = o1;
    dx = d1;

    coeff = sf_floatalloc2(n-1,4);
    /* rename for convenience */
    b = coeff[1];
    c = coeff[2];
    d = coeff[3];

    slv = tridiagonal_init (n-2);
    tridiagonal_const_define (slv,4.*d1,d1);
}

void spline3_close (void) {
    free (h);
    free (a);
    free (*coeff);
    free (coeff);
    tridiagonal_close (slv);
}

void spline3_close1 (void) {
    free (*coeff);
    free (coeff);
    tridiagonal_close (slv);
}

void spline_coeffs(float** table)
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
	a[k] = 2.*(h[k+1] + h[k]);         /* diagonal */
	c[k+1] = b[k+1] - b[k];            /* right-hand side */
    }
    c[0] = 0;

    /* solve the tridiagonal system */
    tridiagonal_define (slv,a,h);
    tridiagonal_solve(slv,c+1);

    for (k=0; k < n-1; k++) {
	if (k < n-2) {
	    d[k] = (c[k+1]-c[k])/h[k];
	    b[k] -= (c[k+1]+2.*c[k])*h[k];
	} else {
	    d[k] = -c[k]/h[k];
	    b[k] -= 2.*c[k]*h[k];
	}
	c[k] *= 3.;
    }
}

void spline_coeffs1(float* table1)
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
    tridiagonal_solve(slv,c+1);

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

/* Function: spline_eval
   ---------------------
   Evaluate a cubic spline
   y - where to evaluate
*/
float spline_eval(float y)
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

/* 	$Id: spline3.c,v 1.4 2003/10/01 22:45:56 fomels Exp $	 */
