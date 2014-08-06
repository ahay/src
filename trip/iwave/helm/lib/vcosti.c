/* vcosti.f -- translated by f2c (version 20100827).
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

/* Table of constant values */

static float c_b2 = 1.f;

/* areful! Anything free comes with no guarantee. */
/* Subroutine */ int vcosti_(integer *n, float *wsave)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sin(double), cos(double);

    /* Local variables */
    static integer k, kc;
    static float fk, dt, pi;
    static integer nm1, np1, ns2;
    extern double pimach_(float *);
    extern /* Subroutine */ int vrffti_(integer *, float *);

/* ***BEGIN PROLOGUE  VCOSTI */
/* ***DATE WRITTEN   860701   (YYMMDD) */
/* ***REVISION DATE  900509   (YYMMDD) */
/* ***CATEGORY NO.  J1A3 */
/* ***KEYWORDS  FAST FOURIER TRANSFORM, COSINE TRANSFORM, MULTIPLE */
/*             SEQUENCES */
/* ***AUTHOR  BOISVERT, R. F. (NIST) */
/* ***PURPOSE  Initialize for VCOST. */
/* ***DESCRIPTION */

/*  Subroutine VCOSTI initializes the array WSAVE which is used in */
/*  subroutine VCOST.  The prime factorization of N together with */
/*  a tabulation of the trigonometric functions are computed and */
/*  stored in WSAVE. */

/*  Input Parameter */

/*  N       the length of the sequence to be transformed.  The method */
/*          is most efficient when N-1 is a product of small primes. */

/*  Output Parameter */

/*  WSAVE   a work array which must be dimensioned at least 3*N+15. */
/*          Different WSAVE arrays are required for different values */
/*          of N.  The contents of WSAVE must not be changed between */
/*          calls of VCOST. */

/*  ----------------------------------------------------------------- */

/*  VCOSTI is a straightforward extension of the subprogram COSTI to */
/*  handle M simultaneous sequences.  COSTI was originally developed */
/*  by P. N. Swarztrauber of NCAR. */

/* ***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel */
/*               Computations, (G. Rodrigue, ed.), Academic Press, 1982, */
/*               pp. 51-83. */
/* ***ROUTINES CALLED  VRFFTI */
/* ***END PROLOGUE  VCOSTI */
/* ***FIRST EXECUTABLE STATEMENT  VCOSTI */
    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    pi = pimach_(&c_b2);
    if (*n <= 3) {
	return 0;
    }
    nm1 = *n - 1;
    np1 = *n + 1;
    ns2 = *n / 2;
    dt = pi / (float) nm1;
    fk = 0.f;
    i__1 = ns2;
    for (k = 2; k <= i__1; ++k) {
	fk += 1.f;
	wsave[k] = sin(fk * dt) * 2.f;
/* L101: */
    }
    fk = 0.f;
    i__1 = ns2;
    for (k = 2; k <= i__1; ++k) {
	kc = np1 - k;
	fk += 1.f;
	wsave[kc] = cos(fk * dt) * 2.f;
/* L102: */
    }
    vrffti_(&nm1, &wsave[*n + 1]);
    return 0;
} /* vcosti_ */

