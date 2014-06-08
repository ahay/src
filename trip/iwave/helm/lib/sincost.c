/* sincost.f -- translated by f2c (version 20100827).
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

static integer c__1 = 1;

/* Subroutine */ int sincost_(integer *m, integer *n, float *x, float *xt, float 
	*wsave)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m1; /* n1; */
    extern /* Subroutine */ int f2c_scopy(integer *, float *, integer *, float *, 
	    integer *), vcost_(integer *, integer *, float *, float *, integer *
	    , float *), vsint_(integer *, integer *, float *, float *, integer *,
	     float *), vcosti_(integer *, float *), vsinti_(integer *, float *);

/* ========================================================================== */

/* W. W. Symes, 2/91 */

/* Interface with vfftpk routines to implement 2d fast sine/cosine transform */
/* of an input array x(1:m,1:n), stored as a linear array x(*). */
/* Uses workspace arrays xt(*) and wsave(*). */

/* The COLUMNS of x are sine-transformed; the ROWS are */
/* cosine-transformed. */

/* VERY IMPORTANT: it is assumed that the columns of the */
/* input array x satisfy Dirichlet conditions at the bottom, i.e. */
/* that x(m,:)=0. Thus the nonredundant part of x is, say, */
/* the first (m-1) rows. Therefore x transpose is passed to */
/* the vfftpak routine VSINT as an nx(m-1) array, with */
/* mdimx = m. This automatically takes care of the need in */
/* VSINT for an extra column of workspace. As a result, the last */
/* column of x transpose (haveing been used as workspace by VSINT) must be */
/* zeroed out again. */

/* ALSO VERY IMPORTANT: the vfftpk routines VSINT and VCOST are */
/* efficient if the ROW length of the input (i.e. the number of */
/* columns) is N+1, resp. N-1, where N is a product of small */
/* primes. Thus if x has dimensions m,n, m and n-1 should be */
/* products of small primes. The obvious way to arrange this, */
/* if say x lives somewhere else as a 2^k x 2^j matrix, is to */
/* pad x with an extra copy of its last column (this is consistent */
/* with the Neumann condition in the column direction). */

/* ========================================================================= */

    /* Parameter adjustments */
    --wsave;
    --xt;
    --x;

    /* Function Body */
    m1 = *m - 1;
/*    n1 = *n - 1; */
/* transform the rows if n > 1: */
    if (*n > 1) {
	vcosti_(n, &wsave[1]);
	vcost_(&m1, n, &x[1], &xt[1], m, &wsave[1]);
    }
/* transpose x into xt */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f2c_scopy(n, &x[i__], m, &xt[(i__ - 1) * *n + 1], &c__1);
    }
/* transform the columns */
    vsinti_(&m1, &wsave[1]);
    vsint_(n, &m1, &xt[1], &x[1], n, &wsave[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xt[m1 * *n + i__] = 0.f;
    }
/* transpose back */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	f2c_scopy(m, &xt[i__], n, &x[(i__ - 1) * *m + 1], &c__1);
    }
    return 0;
} /* sincost_ */

