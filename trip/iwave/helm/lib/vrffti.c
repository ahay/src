/* vrffti.f -- translated by f2c (version 20100827).
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
/* Subroutine */ int vrffti_(integer *n, float *wsave)
{
    extern /* Subroutine */ int vrfti1_(integer *, float *, float *);

/* ***BEGIN PROLOGUE  VRFFTI */
/* ***DATE WRITTEN   860701   (YYMMDD) */
/* ***REVISION DATE  900509   (YYMMDD) */
/* ***CATEGORY NO.  J1A1 */
/* ***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM, */
/*             MULTIPLE SEQUENCES */
/* ***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST) */
/* ***PURPOSE  Initialization for VRFFTF and VRFFTB. */
/* ***DESCRIPTION */

/*  Subroutine VRFFTI initializes the array WSAVE which is used in */
/*  both VRFFTF and VRFFTB.  The prime factorization of N together with */
/*  a tabulation of certain trigonometric functions are computed and */
/*  stored in the array WSAVE. */

/*  Input Parameter */

/*  N       the length of the sequence to be transformed.  There is no */
/*          restriction on N. */

/*  Output Parameter */

/*  WSAVE   a work array which must be dimensioned at least N+15. */
/*          The same work array can be used for both VRFFTF and VRFFTB */
/*          as long as N remains unchanged.  Different WSAVE arrays */
/*          are required for different values of N.  The contents of */
/*          WSAVE must not be changed between calls of VRFFTF or VRFFTB. */


/*              * * * * * * * * * * * * * * * * * * * * * */
/*              *                                       * */
/*              *         PROGRAM SPECIFICATIONS        * */
/*              *                                       * */
/*              * * * * * * * * * * * * * * * * * * * * * */


/*     DIMENSION OF    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15) */
/*     ARGUMENTS */

/*     LATEST          AUGUST 1, 1985 */
/*     REVISION */

/*     SUBPROGRAMS     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3, */
/*     REQUIRED        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2, */
/*                     VRADB3, VRADB4, VRADB5, VRADBG, PIMACH */

/*     SPECIAL         NONE */
/*     CONDITIONS */

/*     COMMON          NONE */
/*     BLOCKS */

/*     I/O             NONE */

/*     PRECISION       SINGLE */

/*     SPECIALIST      ROLAND SWEET */

/*     LANGUAGE        FORTRAN */

/*     HISTORY         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE */
/*                     NATIONAL BUREAU OF STANDARDS (BOULDER). */

/*     ALGORITHM       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION */
/*                     OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM. */

/*     PORTABILITY     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77. */
/*                     THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN */
/*                     THE FUNCTION PIMACH. */

/*     REQUIRED        COS,SIN */
/*     RESIDENT */
/*     ROUTINES */


/* ***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel */
/*               Computations, (G. Rodrigue, ed.), Academic Press, 1982, */
/*               pp. 51-83. */
/* ***ROUTINES CALLED  VRFTI1 */
/* ***END PROLOGUE  VRFFTI */

/*     VRFFTPK, VERSION 1, AUGUST 1985 */

/* ***FIRST EXECUTABLE STATEMENT  VRFFTI */
    /* Parameter adjustments */
    --wsave;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
    vrfti1_(n, &wsave[1], &wsave[*n + 1]);
    return 0;
} /* vrffti_ */

