#include <rsf.h> 

#include "spline3.h"
#include "tridiagonal.h"

static int n;
static float *h, *a, *b, *c, *d, **coeff, *x;
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

void spline3_close (void) {
    free (h);
    free (a);
    free (*coeff);
    free (coeff);
    tridiagonal_close (slv);
}

/* Function: spline_coeffs
   -----------------------
   Compute spline coefficients for interpolating natural cubic spline
   n - number of knots
   x[n] - knots
   f[n] - function values
   coeff[4][n] - coefficients
*/
void spline_coeffs(const float** table)
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
