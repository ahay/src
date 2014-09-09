/* vsint.f -- translated by f2c (version 20100827).
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
/* Subroutine */ int vsint_(integer *m, integer *n, float *x, float *xt, 
	integer *mdimx, float *wsave)
{
    /* System generated locals */
    integer x_dim1, x_offset, xt_dim1, xt_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static integer i__, j, k;
    static float t1, t2;
    static integer kc, nf;
    static float xh;
    static integer np1, ns2, modn;
    static float scale, sqrth;
    extern /* Subroutine */ int vrfftf_(integer *, integer *, float *, float *, 
	    integer *, float *);

/* ***BEGIN PROLOGUE  VSINT */
/* ***DATE WRITTEN   860701   (YYMMDD) */
/* ***REVISION DATE  900509   (YYMMDD) */
/* ***CATEGORY NO.  J1A3 */
/* ***KEYWORDS  FAST FOURIER TRANSFORM, SINE TRANSFORM, MULTIPLE */
/*             SEQUENCES */
/* ***AUTHOR  BOISVERT, R. F., (NIST) */
/* ***PURPOSE  Sine transform of one or more float, odd sequences. */
/* ***DESCRIPTION */

/*  Subroutine VSINT computes the discrete Fourier sine transform */
/*  of M odd sequences X(J,I), J=1,...,M.  The transform is defined */
/*  below at output parameter X. */

/*  The array WSAVE which is used by subroutine VSINT must be */
/*  initialized by calling subroutine VSINTI(N,WSAVE). */

/*  Input Parameters */

/*  M       the number of sequences to be transformed. */

/*  N       the length of the sequence to be transformed.  The method */
/*          is most efficient when N+1 is the product of small primes. */

/*  X       an array of size at least X(MDIMX,N+1) which contains the */
/*          the sequences to be transformed.  The sequences are stored */
/*          in the ROWS of X.  Thus, the Jth sequence is stored in */
/*          X(J,I), I=1,..,N.  The extra column of X is used as work */
/*          storage. */

/*  XT      a work array of size at least XT(MDIMX,N+1). */

/*  MDIMX   the first dimension of the array X exactly as it appears in */
/*          the calling program. */

/*  WSAVE   a work array with dimension at least INT(2.5*N+15) */
/*          in the program that calls VSINT.  The WSAVE array must be */
/*          initialized by calling subroutine VSINTI(N,WSAVE), and a */
/*          different WSAVE array must be used for each different */
/*          value of N.  This initialization does not have to be */
/*          repeated so long as N remains unchanged. */

/*  Output Parameters */

/*  X       for I=1,...,N and J=1,...,M */

/*               X(J,I)= the sum from k=1 to k=N */

/*                    2*X(J,K)*SIN(K*I*PI/(N+1))/SQRT(2*(N+1)) */

/*  WSAVE   contains initialization calculations which must not be */
/*          destroyed between calls of VSINT. */

/*  ----------------------------------------------------------------- */

/*  NOTE  -  A call of VSINT followed immediately by another call */
/*           of VSINT will return the original sequences X.  Thus, */
/*           VSINT is the correctly normalized inverse of itself. */

/*  ----------------------------------------------------------------- */

/*  VSINT is a straightforward extension of the subprogram SINT to */
/*  handle M simultaneous sequences.  The scaling of the sequences */
/*  computed by VSINT is different than that of SINT.  SINT was */
/*  originally developed by P. N. Swarztrauber of NCAR. */

/* ***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel */
/*               Computations, (G. Rodrigue, ed.), Academic Press, 1982, */
/*               pp. 51-83. */
/* ***ROUTINES CALLED  VRFFTF */
/* ***END PROLOGUE  VSINT */
/* ***FIRST EXECUTABLE STATEMENT  SINT */
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
    if (*n > 2) {
	goto L300;
    }

/*  CASE   N = 2 */

    sqrth = sqrt(.5f);
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	xh = sqrth * (x[j + x_dim1] + x[j + (x_dim1 << 1)]);
	x[j + (x_dim1 << 1)] = sqrth * (x[j + x_dim1] - x[j + (x_dim1 << 1)]);
	x[j + x_dim1] = xh;
/* L201: */
    }
    goto L900;

/*  CASE   N .GT. 2 */

/*     ... PREPROCESSING */

L300:
    np1 = *n + 1;
    ns2 = *n / 2;
    i__1 = *m;
    for (j = 1; j <= i__1; ++j) {
	xt[j + xt_dim1] = 0.f;
/* L301: */
    }
    i__1 = ns2;
    for (k = 1; k <= i__1; ++k) {
	kc = np1 - k;
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    t1 = x[j + k * x_dim1] - x[j + kc * x_dim1];
	    t2 = wsave[k] * (x[j + k * x_dim1] + x[j + kc * x_dim1]);
	    xt[j + (k + 1) * xt_dim1] = t1 + t2;
	    xt[j + (kc + 1) * xt_dim1] = t2 - t1;
/* L310: */
	}
    }
    modn = *n % 2;
    if (modn != 0) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    xt[j + (ns2 + 2) * xt_dim1] = x[j + (ns2 + 1) * x_dim1] * 4.f;
/* L320: */
	}
    }

/*     ... REAL PERIODIC TRANSFORM */

    nf = ns2 + 1;
    vrfftf_(m, &np1, &xt[xt_offset], &x[x_offset], mdimx, &wsave[nf]);

/*     ... POSTPROCESSING */

    i__2 = *m;
    for (j = 1; j <= i__2; ++j) {
	x[j + x_dim1] = xt[j + xt_dim1] * .5f;
/* L330: */
    }
    i__2 = *n;
    for (i__ = 3; i__ <= i__2; i__ += 2) {
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    x[j + (i__ - 1) * x_dim1] = -xt[j + i__ * xt_dim1];
/* L340: */
	}
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    x[j + i__ * x_dim1] = x[j + (i__ - 2) * x_dim1] + xt[j + (i__ - 1)
		     * xt_dim1];
/* L345: */
	}
/* L350: */
    }
    if (modn == 0) {
	i__2 = *m;
	for (j = 1; j <= i__2; ++j) {
	    x[j + *n * x_dim1] = -xt[j + (*n + 1) * xt_dim1];
/* L360: */
	}
    }

/*     ... NORMALIZATION */

    scale = sqrt(.5f);
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
	i__1 = *m;
	for (j = 1; j <= i__1; ++j) {
	    x[j + i__ * x_dim1] = scale * x[j + i__ * x_dim1];
/* L370: */
	}
    }

/*  EXIT */

L900:
    return 0;
} /* vsint_ */

