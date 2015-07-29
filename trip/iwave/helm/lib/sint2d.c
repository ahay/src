/* sint2d.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int sint2d_(integer *m, integer *n, float *x, float *xt, float *
	wsave)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m1, n1;
    extern /* Subroutine */ int f2c_scopy(integer *, float *, integer *, float *, 
	    integer *), vsint_(integer *, integer *, float *, float *, integer *
	    , float *), vsinti_(integer *, float *);


/* Interface with vfftpk routines to implement 2d fast sine transform */
/* of an input array x(1:m,1:n), stored as a linear array x(*). */
/* Uses workspace arrays xt(*) and wsave(*). */

/* VERY IMPORTANT: it is assumed that the columns of the */
/* input array x satisfy Dirichlet conditions on the bottom and */
/* on the right both ends, i.e. that x(:,n)=x(m,:)=0. */
/* Thus the nonredundant part of x is, the (m-1)x(n-1) principal */
/* submatrix. Therefore x is passed to */
/* the vfftpak routine VSINT as an (m-1)x(n-1) array, with */
/* mdimx = m. This automatically takes care of the need in */
/* VSINT for an extra column of workspace. As a result, the last */
/* column of x (haveing been used as workspace by VSINT) must be */
/* zeroed out again. */




    /* Parameter adjustments */
    --wsave;
    --xt;
    --x;

    /* Function Body */
    m1 = *m - 1;
    n1 = *n - 1;

/* transform the rows */

    if (*n > 1) {
	vsinti_(&n1, &wsave[1]);
	vsint_(m, &n1, &x[1], &xt[1], m, &wsave[1]);
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    x[n1 * *m + i__] = 0.f;
	}
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
} /* sint2d_ */

