/* getcoefn.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Subroutine */ int getcoefn_(integer *n1, float *d1, float *w1, integer *n2, 
	float *d2, float *w2, float *coef, integer *ier)
{
    /* Initialized data */

    static float tol = 1e-10f;

    /* System generated locals */
    integer coef_dim1, coef_offset, i__1, i__2;

    /* Builtin functions */
    double cos(double);

    /* Local variables */
    static integer i1, i2;
    static float t1, t2, fred;

/* ========================================================== */

/*  Sets up Fourier Sine/Cosine Transform Multiplier Array for */
/*  discrete Helmholtz Operator with Dirichlet conditions on */
/*  top and bottom of rectangular grid and Neumann conditions */
/*  on sides. */

/*  The operator is */

/*    u |---> u + w1^2 d^2 u/dx1^2 + w2^2 d^2 u/dx2^2 */

/*  Designed for use with vfftpak vectorized FFTs. */

/* ========================================================== */

/*      include 'system.h' */
    /* Parameter adjustments */
    coef_dim1 = *n1;
    coef_offset = 1 + coef_dim1;
    coef -= coef_offset;

    /* Function Body */
    if (*ier != 0) {
	return 0;
    }
    if (dabs(*d1) < tol) {
      *ier=123;
	return 0;
    }
/*         write(ipdmp,*)' Error: GETCOEFN' */
/*         write(ipdmp,*)' d1 = ',d1 */
/*         write(ipdmp,*)' less than tolerance currently hardwired' */
/*         write(ipdmp,*)' into this routine: tol = ',tol */
/*         ier=123 */
/*         return */
/*      end if */
    if (dabs(*d2) < tol) {
      *ier=123;
	return 0;
    }
/*         write(ipdmp,*)' Error: GETCOEFN' */
/*         write(ipdmp,*)' d2 = ',d2 */
/*         write(ipdmp,*)' less than tolerance currently hardwired' */
/*         write(ipdmp,*)' into this routine: tol = ',tol */
/*         ier=123 */
/*         return */
/*      end if */
    fred = 3.1415927f;
    i__1 = *n1;
    for (i1 = 1; i1 <= i__1; ++i1) {
	t1 = *w1 * 2 * (1.f - cos(fred * i1 / *n1)) / (*d1 * *d1);
	i__2 = *n2;
	for (i2 = 1; i2 <= i__2; ++i2) {
	    coef[i1 + i2 * coef_dim1] = t1;
	}
    }
/* note that sum here starts at 2 because the first cosine transform */
/* coefficient is 0 (dc component!). note that this convention also */
/* neatly takes care of the possible difficulty with n2=1. */
    i__1 = *n2;
    for (i2 = 2; i2 <= i__1; ++i2) {
	t2 = *w2 * 2 * (1.f - cos(fred * (i2 - 1) / (*n2 - 1))) / (*d2 * *d2);
	i__2 = *n1;
	for (i1 = 1; i1 <= i__2; ++i1) {
	    coef[i1 + i2 * coef_dim1] += 1.0f + t2;
	}
    }
    return 0;
} /* getcoefn_ */

