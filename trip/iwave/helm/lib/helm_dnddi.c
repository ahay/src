/* helm_dnddi.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

/* modified WWS 10.12.11 to use pow instead of libf2c version - thus this package
   is completely independent of libf2c, and need not link to it.
*/

/* include C math library headers to avoid having to use libf2c version of pow */
#include <math.h>

#include "f2c.h"

/* Subroutine */ int helm_(integer *nsx, integer *ntx, float *dt, float *dx, 
	float *st, float *sx, float *p, float *d__, float *x, float *y, float *work, 
	integer *lenwork, integer *ier)
{

    /* System generated locals */
    integer i__1, i__2;
    double d__1, d__2;

    /* Builtin functions */
    /*    double pow_dd(double *, double *); */
    /* replace with pow from C math library - OK as double=double */

    /* Local variables */
    static integer ptr_to_x__, i__, j, ptr_to_xt__, nd, jx, jy, ptr_to_coef__,
	/* ntx1, */ nsx2, nfft, ptr_to_wsave__;
    extern /* Subroutine */ int sint2d_(integer *, integer *, float *, float *, 
	    float *), getcoef_(integer *, float *, float *, integer *, float *, 
	    float *, float *, integer *);
    static integer stridex;


/* 081106: Modified to eliminate Fortran i/o. name changed WWS */
/* ------------------------------------------------------------ */

/* Incore powers of 2d Helmholtz operator, defined with */
/* Dirichlet conditions on top and side boundaries, Neumann */
/* condition on the bottom. The top boundary is d depth units */
/* below the top of the grid. */

/* Implements */

/*     y = (I - D^t*diag(s)*D)^p x */

/* Here s >= 0 (scale) is float vector of scales in the two spatial directions, */
/* p is a float scalar power, and D and D^t are the gradient and divergence */
/* operators respectively. The scaled Discrete Laplacian (D^t*diag(s)*D) */
/* with the above-mentioned boundary conditions is implemented via DFT. */

/* METHOD: */
/* First the columns are cut off at the NEXT row of data below depth d, */
/* then evenly extended, and an extra row of zeroes is added at the bottom. */
/* The last row is used as */
/* workspace in SINT2D; setting it to zero is consistent with the */
/* Dirichlet condition. Note that the input columns are viewed as */
/* the nonredundant part of input sequences satisfying the Dirichlet */
/* condition at the top, so the returned array will not necessarily */
/* vanish on the top row even though the input array does - the */
/* Dirichlet condition is implicitly imposed on the (virtual) row */
/* above the top. The top row should be small, however. */
/* In the same vein the last column is regarded as redundant, */
/* and filled with zeroes. This is necessary for the sake of */
/* efficiency, but will obviously result in error unless the */
/* assumption is correct. In order for the vfftpk sine transform of */
/* the rows to work efficiently, the input length must be one */
/* less than a power of two (or product of small primes). Since */
/* typically the input length is actually a power of two, ignoring */
/* the last column is essential. */

/* Then the extended array is passed to the routine SINT2D, where */
/* the rows are cosine-transformed, the columns sine-transformed. */
/* Then an appropriate multiplier array is applied. Since SINT2D */
/* implements an idempotent operator, the scaled extended transformed */
/* array is once more passed back to SINT2D, then the image is */
/* extracted from the subarray aligned with the original data. */

/* ------------------------------------------------------------ */

/* ARGUMENTS */
/* input array */
/* output array */
/* work array */
/* input trace step */
/* input sample step */
/* weight in trace direction */
/* weight in sample direction */
/* power */
/* cutoff depth */
/* INTERNAL VARIABLES */
/* length of work array */
/* number of input traces */
/* number of input samples */
/* error flag */

/* ------------------------------------------------------ */

    /* Parameter adjustments */
    --work;
    --y;
    --x;

    /* Function Body */
    if (*ier != 0) {
	return 0;
    }
/* check to make sure that scale factor s is >= 0 */
    if (*sx < 0.f || *st < 0.f) {
	*ier = 44;
	return 0;
    }
/*     useful aux. numbers 
    ntx1 = *ntx - 1; */
/* Computing MAX */
    i__1 = (integer) (*d__ / *dt) + 1;
    nd = max(i__1,0);
    nsx2 = (*nsx - nd) << 1;
    nfft = *ntx * nsx2;
/* nfft < 2*ntx*nsx */
/* ptr_to_wsave  < 1+6*ntx*nsx */
/* lenwork > 6*ntx*nsx+3*max(ntx,2*nsx)+21 */
/* set up pointers into work vector */
    ptr_to_x__ = 1;
    ptr_to_xt__ = nfft + 1;
    ptr_to_coef__ = (nfft << 1) + 1;
    ptr_to_wsave__ = nfft * 3 + 1;
    if (ptr_to_wsave__ + (max(*ntx,nsx2) * 3 + 20) > *lenwork) {
	*ier = 767;
	return 0;
    } else {
	i__1 = ptr_to_wsave__ + (max(*ntx,nsx2) * 3 + 20);
	for (i__ = 1; i__ <= i__1; ++i__) {
	    work[i__] = 0.f;
	}
    }
    getcoef_(&nsx2, dt, st, ntx, dx, sx, &work[ptr_to_coef__], ier);
    if (*ier != 0) {
	return 0;
    }
/* copy the data onto the workvector, extending the columns evenly. */

/* input       output */

/*  j            j, nsx2-j    j = nd+1:nsx */
/*  nsx2         (zeroes) */

    stridex = *nsx;
    i__1 = *ntx;
    for (j = 1; j <= i__1; ++j) {
	jx = (j - 1) * stridex;
	jy = (j - 1) * nsx2 + ptr_to_x__ - 1;
	i__2 = *nsx - nd;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    work[jy + i__] = x[jx + i__ + nd];
	    work[jy + nsx2 - i__] = x[jx + i__ + nd];
	}
	work[jy + nsx2] = 0.f;
    }

/* transform, Fourier multiplier, inverse transform */

    sint2d_(&nsx2, ntx, &work[ptr_to_x__], &work[ptr_to_xt__], &work[
	    ptr_to_wsave__]);
    i__1 = nfft;
    for (j = 1; j <= i__1; ++j) {
      d__1 = (double) work[ptr_to_coef__ - 1 + j];
      d__2 = (double) (*p);
      /*	work[ptr_to_x__ - 1 + j] *= pow_dd(&d__1, &d__2); */
      work[ptr_to_x__ - 1 + j] *= pow(d__1, d__2);
    }
    sint2d_(&nsx2, ntx, &work[ptr_to_x__], &work[ptr_to_xt__], &work[
	    ptr_to_wsave__]);
/* zero out last column */
    if (*ntx > 1) {
	i__1 = nsx2;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    work[ptr_to_x__ - 1 + (*ntx - 1) * nsx2 + i__] = 0.f;
	}
    }
/* copy output data into output buffer */
    i__1 = *ntx;
    for (j = 1; j <= i__1; ++j) {
	jx = (j - 1) * stridex;
	jy = (j - 1) * nsx2 + ptr_to_x__ - 1;
	i__2 = nd;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[jx + i__] = 0.f;
	}
	i__2 = *nsx - nd;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[jx + i__ + nd] = work[jy + i__];
	}
	work[jy + nsx2] = 0.f;
    }
    return 0;
} /* helm_ */

