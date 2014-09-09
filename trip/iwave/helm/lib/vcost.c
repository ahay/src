/* vcost.f -- translated by f2c (version 20100827).
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

/* areful! Anything free comes with no guarantee. */
/* Subroutine */ int vcost_(integer *m, integer *n, float *x, float *xt, 
	integer *mdimx, float *wsave)
{
    /* System generated locals */
    integer x_dim1, x_offset, xt_dim1, xt_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static integer i__, j, k;
    static float t1, t2;
    static integer kc;
    static float xi;
    static integer nm1, np1;
    static float x1h;
    static integer ns2;
    static float tx2, x1p3;
    static integer modn;
    static float scale, factor;
    extern /* Subroutine */ int vrfftf_(integer *, integer *, float *, float *, 
	    integer *, float *);

/* ***BEGIN PROLOGUE  VCOST */
/* ***DATE WRITTEN   860701   (YYMMDD) */
/* ***REVISION DATE  900509   (YYMMDD) */
/* ***CATEGORY NO.  J1A3 */
/* ***KEYWORDS  FAST FOURIER TRANSFORM, COSINE TRANSFORM, MULTIPLE */
/*             SEQUENCES */
/* ***AUTHOR  BOISVERT, R. F. (NIST) */
/* ***PURPOSE  Cosine transform of one or more real, even sequences. */
/* ***DESCRIPTION */

/*  Subroutine VCOST computes the discrete Fourier cosine transform */
/*  of M even sequences X(J,I), J=1,...,M.  The transform is defined */
/*  below at output parameter X. */

/*  The array WSAVE which is used by subroutine VCOST must be */
/*  initialized by calling subroutine VCOSTI(N,WSAVE). */

/*  Input Parameters */

/*  M       the number of sequences to be transformed. */

/*  N       the length of the sequence to be transformed.  N must be */
/*          greater than 1.  The method is most efficient when N-1 is */
/*          is a product of small primes. */

/*  X       an array of size at least X(MDIMX,N) which contains the */
/*          the sequences to be transformed.  The sequences are stored */
/*          in the ROWS of X.  Thus, the Jth sequence is stored in */
/*          X(J,I), I=1,..,N. */

/*  XT      a work array of size at least XT(MDIMX,N-1). */

/*  MDIMX   the first dimension of the array X exactly as it appears in */
/*          the calling program. */

/*  WSAVE   a work array which must be dimensioned at least 3*N+15 */
/*          in the program that calls VCOST.  The WSAVE array must be */
/*          initialized by calling subroutine VCOSTI(N,WSAVE), and a */
/*          different WSAVE array must be used for each different */
/*          value of N.  This initialization does not have to be */
/*          repeated so long as N remains unchanged.  Thus subsequent */
/*          transforms can be obtained faster than the first. */

/*  Output Parameters */

/*  X       For I=1,...,N and J=1,...,M */

/*             X(J,I) = ( X(J,1)+(-1)**(I-1)*X(J,N) */

/*               + the sum from K=2 to K=N-1 */

/*                 2*X(J,K)*COS((K-1)*(I-1)*PI/(N-1)) )/SQRT(2*(N-1)) */

/*  WSAVE   contains initialization calculations which must not be */
/*          destroyed between calls of VCOST. */

/*  ----------------------------------------------------------------- */

/*  NOTE  -  A call of VCOST followed immediately by another call */
/*           of VCOST will return the original sequences X.  Thus, */
/*           VCOST is the correctly normalized inverse of itself. */

/*  ----------------------------------------------------------------- */

/*  VCOST is a straightforward extension of the subprogram COST to */
/*  handle M simultaneous sequences.  The scaling of the sequences */
/*  computed by VCOST is different than that of COST.  COST was */
/*  originally developed by P. N. Swarztrauber of NCAR. */

/* ***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel */
/*               Computations, (G. Rodrigue, ed.), Academic Press, 1982, */
/*               pp. 51-83. */
/* ***ROUTINES CALLED  VRFFTF */
/* ***END PROLOGUE  VCOST */
/* ***FIRST EXECUTABLE STATEMENT  VCOST */
    /* Parameter adjustments */
    xt_dim1 = *mdimx;
    xt_offset = 1 + xt_dim1;
    xt -= xt_offset;
    x_dim1 = *mdimx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --wsave;

    /* Function Body */
    if (*m <= 0) {
	goto L900;
    }
    if (*n <= 1) {
	goto L900;
    }
    if (*n > 3) {
	goto L400;
    }
    if (*n == 3) {
	goto L300;
    }

/*  CASE  N = 2 */

    scale = sqrt(.5f);
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	x1h = scale * (x[j + x_dim1] + x[j + (x_dim1 << 1)]);
	x[j + (x_dim1 << 1)] = scale * (x[j + x_dim1] - x[j + (x_dim1 << 1)]);
	x[j + x_dim1] = x1h;
/* L210: */
    }
    goto L900;

/*  CASE  N = 3 */

L300:
    scale = .5f;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	x1p3 = x[j + x_dim1] + x[j + x_dim1 * 3];
	tx2 = x[j + (x_dim1 << 1)] + x[j + (x_dim1 << 1)];
	x[j + (x_dim1 << 1)] = scale * (x[j + x_dim1] - x[j + x_dim1 * 3]);
	x[j + x_dim1] = scale * (x1p3 + tx2);
	x[j + x_dim1 * 3] = scale * (x1p3 - tx2);
/* L310: */
    }
    goto L900;

/*  CASE  N .GT. 3 */

/*     ... PREPROCESSING */

L400:
    nm1 = *n - 1;
    np1 = *n + 1;
    ns2 = *n / 2;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	xt[j + xt_dim1] = x[j + x_dim1] - x[j + *n * x_dim1];
	x[j + x_dim1] += x[j + *n * x_dim1];
/* L410: */
    }
    i__1 = ns2;
    for (k = 2; k <= i__1; ++k) {
	kc = np1 - k;
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    t1 = x[j + k * x_dim1] + x[j + kc * x_dim1];
	    t2 = x[j + k * x_dim1] - x[j + kc * x_dim1];
	    xt[j + xt_dim1] += wsave[kc] * t2;
	    t2 = wsave[k] * t2;
	    x[j + k * x_dim1] = t1 - t2;
	    x[j + kc * x_dim1] = t1 + t2;
/* L420: */
	}
    }
    modn = *n % 2;
    if (modn != 0) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    x[j + (ns2 + 1) * x_dim1] += x[j + (ns2 + 1) * x_dim1];
/* L430: */
	}
    }
    i__2 = *m;
    for (j = 1; j <= i__2; ++j) {
	x[j + *n * x_dim1] = xt[j + xt_dim1];
/* L435: */
    }

/*     ... REAL PERIODIC TRANSFORM */

    vrfftf_(m, &nm1, &x[x_offset], &xt[xt_offset], mdimx, &wsave[np1]);

/*     ... POSTPROCESSING */

    factor = 1.f / sqrt((float) nm1);
    i__2 = *m;
    for (j = 1; j <= i__2; ++j) {
	xt[j + xt_dim1] = x[j + (x_dim1 << 1)];
	x[j + (x_dim1 << 1)] = factor * x[j + *n * x_dim1];
/* L440: */
    }
    i__2 = *n;
    for (i__ = 4; i__ <= i__2; i__ += 2) {
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    xi = x[j + i__ * x_dim1];
	    x[j + i__ * x_dim1] = x[j + (i__ - 2) * x_dim1] - x[j + (i__ - 1) 
		    * x_dim1];
	    x[j + (i__ - 1) * x_dim1] = xt[j + xt_dim1];
	    xt[j + xt_dim1] = xi;
/* L450: */
	}
    }
    if (modn != 0) {
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    x[j + *n * x_dim1] = xt[j + xt_dim1];
/* L460: */
	}
    }

/*     ... NORMALIZATION */

    scale = sqrt(.5f);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    x[j + i__ * x_dim1] = scale * x[j + i__ * x_dim1];
/* L490: */
	}
    }

/*  EXIT */

L900:
    return 0;
} /* vcost_ */

