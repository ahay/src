/* vrfftf.f -- translated by f2c (version 20100827).
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
/* Subroutine */ int vrfftf_(integer *m, integer *n, float *r__, float *rt, 
	integer *mdimr, float *wsave)
{
    /* System generated locals */
    integer r_dim1, r_offset, rt_dim1, rt_offset;

    /* Local variables */
    extern /* Subroutine */ int vrftf1_(integer *, integer *, float *, float *, 
	    integer *, float *, float *);

/* ***BEGIN PROLOGUE  VRFFTF */
/* ***DATE WRITTEN   850801   (YYMMDD) */
/* ***REVISION DATE  900509   (YYMMDD) */
/* ***CATEGORY NO.  J1A1 */
/* ***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM, */
/*             FOURIER ANALYSIS, FORWARD TRANSFORM, MULTIPLE SEQUENCES */
/* ***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST) */
/* ***PURPOSE  Forward float periodic transform, M sequences. */
/* ***DESCRIPTION */

/*  Subroutine VRFFTF computes the Fourier coefficients (forward */
/*  transform) of a number of float periodic sequences.  Specifically, */
/*  for each sequence the subroutine claculates the independent */
/*  Fourier coefficients described below at output parameter R. */

/*  The array WSAVE which is used by subroutine VRFFTF must be */
/*  initialized by calling subroutine VRFFTI(N,WSAVE). */


/*  Input Parameters */

/*  M       the number of sequences to be transformed. */

/*  N       the length of the sequences to be transformed.  The method */
/*          is most efficient when N is a product of small primes, */
/*          however n may be any positive integer. */

/*  R       afloat two-dimensional array of size MDIMX x N containing the */
/*          the sequences to be transformed.  The sequences are stored */
/*          in the ROWS of R.  Thus, the I-th sequence to be transformed, */
/*          X(I,J), J=0,1,...,N-1, is stored as */

/*               R(I,J) = X(I,J-1) , J=1, 2, . . . , N. */

/*  RT      a float two-dimensional work array of size MDIMX x N. */

/*  MDIMR   the row (or first) dimension of the arrays R and RT exactly */
/*          as they appear in the calling program.  This parameter is */
/*          used to specify the variable dimension of these arrays. */

/*  WSAVE   a float one-dimensional work array which must be dimensioned */
/*          at least N+15.  The WSAVE array must be initialized by */
/*          calling subroutine VRFFTI.  A different WSAVE array must be */
/*          used for each different value of N.  This initialization does */
/*          not have to be repeated so long as N remains unchanged.  The */
/*          same WSAVE array may be used by VRFFTF and VRFFTB. */

/*  Output Parameters */

/*  R       contains the Fourier coefficients F(K) for each of the M */
/*          input sequences.  Specifically, row I of R, R(I,J), */
/*          J=1,2,..,N, contains the independent Fourier coefficients */
/*          F(I,K), for the I-th input sequence stored as */

/*             R(I,1) = REAL( F(I,0) ), */
/*                    = SQRT(1/N)*SUM(J=0,N-1)[ X(I,J) ], */

/*             R(I,2*K) = REAL( F(I,K) ) */
/*                      = SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*COS(2J*K*PI/N)] */

/*             R(I,2*K+1) = IMAG( F(I,K) ) */
/*                        =-SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*SIN(2J*K*PI/N)] */

/*                   for K = 1, 2, . . . , M-1, */

/*              and, when N is even, */

/*              R(I,N) = REAL( F(I,N/2) ). */
/*                     = SQRT(1/N)*SUM(J=0,N-1)[ (-1)**J*X(I,J) ]. */

/*  WSAVE   contains results which must not be destroyed between calls */
/*          to VRFFTF or VRFFTB. */

/*  ----------------------------------------------------------------- */

/*  NOTE  -  A call of VRFFTF followed immediately by a call of */
/*           of VRFFTB will return the original sequences R.  Thus, */
/*           VRFFTB is the correctly normalized inverse of VRFFTF. */

/*  ----------------------------------------------------------------- */

/*  VRFFTF is a straightforward extension of the subprogram RFFTF to */
/*  handle M simultaneous sequences.  RFFTF was originally developed */
/*  by P. N. Swarztrauber of NCAR. */


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
/* ***ROUTINES CALLED  VRFTF1 */
/* ***END PROLOGUE  VRFFTF */

/*     VRFFTPK, VERSION 1, AUGUST 1985 */

/* ***FIRST EXECUTABLE STATEMENT  VRFFTF */
    /* Parameter adjustments */
    --wsave;
    rt_dim1 = *mdimr;
    rt_offset = 1 + rt_dim1;
    rt -= rt_offset;
    r_dim1 = *mdimr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;

    /* Function Body */
    if (*n == 1) {
	return 0;
    }
    vrftf1_(m, n, &r__[r_offset], &rt[rt_offset], mdimr, &wsave[1], &wsave[*n 
	    + 1]);
    return 0;
} /* vrfftf_ */

