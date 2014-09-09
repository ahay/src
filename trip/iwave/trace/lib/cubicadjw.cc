#include "cubic.h"
#include "utils.h" /* added 06.09.12 WWS */

int cubicadj_getworksize(int ntout, int ntin) { return ntin + 7*ntout; }

/* ORIGINAL FORTRAN DOCUMENTATION */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*     Cubic Spline Adjoint */
/*     by Ann Campbell, Summer 1991 */

/*     function:  This program generates the adjoint to the cubic */
/*                spline and uses it to create an output vector */
/*                with a specific initial value, sample rate, and */
/*                number of data points. */


/* 	the parameters are: */
/* 	x1i:	initial x value for input data */
/* 	hi:	sample rate for input data */
/* 	vi:	function values for input data */
/* 	ni:	number of input data points */
/* 	x1o:	initial x value for output data */
/* 	ho:	sample rate for output data */
/* 	vo:	output vector */
/* 	no:	number of output data points */
/* 	iend:	type of end condition to be used */
/* 		=1	linear ends, s(1)=s(n)=0 */
/* 		=2	parabolic ends,s(1)=s(2),s(n)=s(n-1) */
/* 		=3	cubic ends,s(1),s(n) are extrapolated */
/* 	work:	vector of workspace */
/* 	wl:	length of work vector */
/* 	ierr:	error flag, should equal 0 if no error */
/*       i:      do loop counter */
/*       first:  first component used in gaussian */
/*               elimination */
/*       last:   last component used in gaussian */
/*               elimination */
/*       j,k:    used in intermediate calculations */
/*       interval:the interval of the output vector that contains */
/*               a specific component of the input vector */
/*       trial:  do loop counter for the gaussian elimination */
/*       u:      the current input x value */
/*       dx:     the difference between the x value of a component */
/*               of the input vector and the x value of the output */
/*               vector that is closest to it */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* cubicadjw.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/
/* not needed, nor is link to f2c, since no Fortran intrinsics used 
#include "f2c.h"
*/

/* consts added 27.02.09 WWS */

/* Subroutine */
int cubicadj_(float const *x1i, 
	      float const *hi, 
	      float const *vi, 
	      int const *ni, 
	      float const *x1o, 
	      float const *ho, 
	      float *vo, 
	      int const *no, 
	      int const *iend, 
	      float *work,
	      int const *wl) {
  /* System generated locals */
  int i__1;
  float r__1;

  /* Local variables */
  static int interval, i__, j, k;
  static int ierr=0;
  static float u, dx;
  static int last;
  static float temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8;
  static int trial, first;

  /* 	now check if an appropriate end condition has been provided */

  /* Parameter adjustments */
  --vi;
  --vo;
  --work;

  /* Function Body */

  /*    if 
	(1) number of samples are same for input and output
	(2) initial abscissae are same to within 0.01% of input step
	(3) total abscissae interval lengths (nt*dt) are same to within
	0.01% of input step

	then regard interpolation as copy, to working precision, and
	implement as a copy.

	- WWS 06.09.12
  */
  if ((ni == no) &&
      (iwave_abs(((float)(*ni))*(*hi)-((float)(*no))*(*ho)) < 0.01 * iwave_abs(*hi)) &&
      (iwave_abs(*x1i-*x1o) < 0.01 * iwave_abs(*hi))) {
    for (i__=0;i__<*ni;i__++) vo[i__]=vi[i__];
    return 0;
  }


  if (*iend < 1 || *iend > 3) {
    ierr = 11;
    return ierr;
  }

  /* 	check if there are enough input and output data points for */
  /* 	the program to work properly */

  if (*no < 4) {
    ierr = 12;
    return ierr;
  }
  if (*ni < 3) {
    ierr = 12;
    return ierr;
  }
  /* 	then test to see if the length of the workspace vector */
  /* 	is sufficient to calculate the adjoint. */
  /*       if it is not, set an error code and return to the calling */
  /*       program. */
  /* 	*note:  to work correctly, this program requires the input */
  /* 	vector to have length at least 3 and the output vector to */
  /* 	have length at least 4 */

  if (*wl < *ni + *no * 7 - 8) {
    ierr = 10;
    return ierr;
  }

  /* 	now, set everything to zero */

  i__1 = *wl;
  for (i__ = 1; i__ <= i__1; ++i__) {
    work[i__] = 0.f;
    /* L1: */
  }
  first = 0;
  last = 0;
  interval = 0;
  trial = 0;
  u = 0.f;
  dx = 0.f;

  /*     to calculate the adjoint, the first step is to determine */
  /*     what interval of the output values that each of the input */
  /*     data points lies in. */
  /*     the intervals are stored in work(1 . .ni) */

  interval = 1;
  u = *x1i;
  i__1 = *ni;
  for (i__ = 1; i__ <= i__1; ++i__) {
    if (u >= *x1o + (*no - 1) * *ho) {
      work[i__] = (float) (*no - 1);
      u += *hi;
      goto L10;
    }
    if (u >= *x1o + (interval - 1) * *ho) {
      if (u <= *x1o + interval * *ho) {
	work[i__] = (float) interval;
	u += *hi;
	goto L10;
      }
    }
    interval = 1;
    j = *no + 1;
  L5:
    k = (interval + j) / 2;
    if (u < *x1o + (k - 1) * *ho) {
      j = k;
    } else {
      interval = k;
    }
    if (j > interval + 1) {
      goto L5;
    } else {
      work[i__] = (float) interval;
      u += *hi;
      goto L10;
    }
  L10:
    ;
  }

  /*     each component of the final output vector is the sum of */
  /*     five numbers, all stored in the work vector. */
  /*     three of the five sets of numbers that sum to the final */
  /*     vector are calculated in the next part of the program, using */
  /*     the values calculated above. */
  /*     the first set (a) is stored in work(ni+1..ni+no), the second */
  /*     set (b) is stored in (1+ni+no..ni+2*no), the third set(c) is */
  /*     stored in (1+ni+2*no..ni+3*no). */
  /*     the components of each of these sets initially equals the sum of */
  /*     the elements in the input vector in a specific interval of */
  /*     the output vector multiplied by dx itself(c), dx squared(b), */
  /*     or dx cubed(a). */

  i__1 = *ni;
  for (i__ = 1; i__ <= i__1; ++i__) {
      j = (int) work[i__];
    dx = *x1i + (i__ - 1) * *hi - (*x1o + (j - 1) * *ho);
    /* Computing 3rd power */
    r__1 = dx;
    work[j + *ni] += vi[i__] * (r__1 * (r__1 * r__1));
    /* Computing 2nd power */
    r__1 = dx;
    work[j + *ni + *no] += vi[i__] * (r__1 * r__1);
    work[j + *ni + (*no << 1)] += vi[i__] * dx;
    /* L20: */
  }
  temp1 = work[*ni + *no - 1];
  temp2 = work[*ni + *no * 3 - 1];
  for (i__ = *no - 1; i__ >= 2; --i__) {
    work[i__ + *ni] = work[i__ + *ni - 1] / (*ho * 6) - work[i__ + *ni] / 
      (*ho * 6);
    work[i__ + *ni + *no] /= 2;
    work[i__ + *ni + (*no << 1)] = *ho / 6 * work[*ni + (*no << 1) + i__ 
						  - 1] + *ho * 2 / 6 * work[*ni + (*no << 1) + i__];
    /* L30: */
  }
  work[*ni + 1] *= -1 / (*ho * 6);
  work[*ni + *no] = temp1 / (*ho * 6);
  work[*ni + *no + 1] /= 2;
  work[*ni + (*no << 1)] = 0.f;
  work[*ni + (*no << 1) + 1] = *ho * 2 / 6 * work[*ni + (*no << 1) + 1];
  work[*ni + *no * 3] = *ho / 6 * temp2;

  /*     now I make changes based on the end conditions which */
  /*     transform the length of the sets of numbers to no-2 */
  /*     [(start+2...start+no-1) instead of (start+1...start+no)] */

  /* 	for iend=1:  no changes */
  /* 	for iend=2, */

  if (*iend == 2) {
    work[*ni + 2] = work[*ni + 1] + work[*ni + 2];
    work[*ni + *no - 1] += work[*ni + *no];
    work[*ni + *no + 2] = work[*ni + *no + 1] + work[*ni + *no + 2];
    work[*ni + (*no << 1) - 1] += work[*ni + (*no << 1)];
    work[*ni + (*no << 1) + 2] = work[*ni + (*no << 1) + 1] + work[*ni + (
									  *no << 1) + 2];
    work[*ni + *no * 3 - 1] += work[*ni + *no * 3];
  }
  if (*iend == 3) {
    work[*ni + 2] += work[*ni + 1] * 2;
    work[*ni + 3] -= work[*ni + 1];
    work[*ni + *no - 2] -= work[*ni + *no];
    work[*ni + *no - 1] += work[*ni + *no] * 2;
    work[*ni + *no + 2] += work[*ni + *no + 1] * 2;
    work[*ni + *no + 3] -= work[*ni + *no + 1];
    work[*ni + (*no << 1) - 2] -= work[*ni + (*no << 1)];
    work[*ni + (*no << 1) - 1] += work[*ni + (*no << 1)] * 2;
    work[*ni + (*no << 1) + 2] += work[*ni + (*no << 1) + 1] * 2;
    work[*ni + (*no << 1) + 3] -= work[*ni + (*no << 1) + 1];
    work[*ni + *no * 3 - 2] -= work[*ni + *no * 3];
    work[*ni + *no * 3 - 1] += work[*ni + *no * 3] * 2;
  }

  /*     for the three sets of numbers that I have been dealing with, */
  /*     [work(ni+1..ni+no),work(ni+no+1..ni+2*no) and work(ni+2*no+1.. */
  /*     ni+3*no)], the next step is for each of them to be multiplied */
  /*     by the inverse of a certain matrix.  Since I know the matrix, */
  /*     rather than its inverse, gaussian elimination is used. */

  /*     the components of the tridiagonal matrix are stored in */
  /*     work(ni+3*no+1..ni+4*no-2),work(ni+4*no-1..ni+5*no-4), */
  /*     and work(ni+5*no-3..ni+6*no-6). */

  /*     the 3 trials represent the multiplication of the inverse */
  /*     on the matrix on the three sets of numbers.  In each trial, */
  /*     the set of numbers is stored in work(ni+6*no-5...ni+7*no-8). */

  for (trial = 1; trial <= 3; ++trial) {
    i__1 = *no - 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
      work[*ni + *no * 3 + i__] = *ho;
      work[*ni + (*no << 2) - 2 + i__] = *ho * 4;
      work[*ni + *no * 5 - 4 + i__] = *ho;
      if (trial == 1) {
	work[*ni + *no * 6 - 6 + i__] = work[i__ + 1 + *ni];
      } else if (trial == 2) {
	work[*ni + *no * 6 - 6 + i__] = work[i__ + 1 + *ni + *no];
      } else if (trial == 3) {
	work[*ni + *no * 6 - 6 + i__] = work[i__ + 1 + *ni + (*no << 
							      1)];
      }
      /* L40: */
    }
    first = 2;
    last = *no - 2;

    /*     again, I need to make adjustments for end conditions */

    /*     iend=1  no changes */

    if (*iend == 2) {
      work[*ni + 1 + (*no << 2) - 2] += *ho;
      work[*ni + *no * 5 - 4] += *ho;
    } else if (*iend == 3) {
      work[*ni + (*no << 2) - 1] = *ho * 6;
      work[*ni + *no * 3 + 2] = 0.f;
      work[*ni + *no * 6 - 7] = 0.f;
      work[*ni + *no * 5 - 4] = *ho * 6;
    }

    /*     this is the gaussian elimination */

    i__1 = last;
    for (i__ = first; i__ <= i__1; ++i__) {
      work[*ni + *no * 3 + i__] /= work[*ni + (*no << 2) + i__ - 3];
      work[*ni + (*no << 2) - 2 + i__] -= work[*ni + *no * 3 + i__] * 
	work[*ni + *no * 5 + i__ - 5];
      work[*ni + *no * 6 - 6 + i__] -= work[*ni + *no * 3 + i__] * work[
									*ni + *no * 6 + i__ - 7];
      /* L60: */
    }
    work[last + *ni + *no * 6 - 6] /= work[last + *ni + (*no << 2) - 2];
    i__1 = first - 1;
    for (j = last - 1; j >= i__1; --j) {
      work[j + *ni + *no * 6 - 6] = (work[j + *ni + *no * 6 - 6] - work[
									j + *ni + *no * 5 - 4] * work[j + *ni + *no * 6 - 5]) / 
	work[j + *ni + (*no << 2) - 2];
      /* L70: */
    }
    i__1 = last;
    for (i__ = first - 1; i__ <= i__1; ++i__) {
      if (trial == 1) {
	work[*ni + i__] = work[*ni + *no * 6 - 6 + i__];
      } else if (trial == 2) {
	work[*ni + *no + i__] = work[*ni + *no * 6 - 6 + i__];
      } else if (trial == 3) {
	work[*ni + (*no << 1) + i__] = work[*ni + *no * 6 - 6 + i__];
      }
      /* L80: */
    }
    /* L100: */
  }
  /*     to get the final values for the three sets of numbers, */
  /*     each of the sets must be multiplied by the same matrix. */
  /*     this is accomplished in the following steps. */

  temp3 = work[*ni + *no - 2];
  temp4 = work[*ni + (*no << 1) - 2];
  temp5 = work[*ni + *no * 3 - 2];
  temp6 = work[*ni + *no - 3];
  temp7 = work[*ni + (*no << 1) - 3];
  temp8 = work[*ni + *no * 3 - 3];
  for (i__ = *no - 2; i__ >= 3; --i__) {
    work[*ni + i__] = 6 / *ho * work[*ni + i__ - 2] - 12 / *ho * work[*ni 
								      + i__ - 1] + 6 / *ho * work[*ni + i__];
    work[*ni + *no + i__] = 6 / *ho * work[*ni + *no + i__ - 2] - 12 / *
      ho * work[*ni + *no + i__ - 1] + 6 / *ho * work[*ni + *no + 
						      i__];
    work[*ni + (*no << 1) + i__] = 6 / *ho * work[*ni + (*no << 1) + i__ 
						  - 2] - 12 / *ho * work[*ni + (*no << 1) + i__ - 1] + 6 / *ho *
      work[*ni + (*no << 1) + i__];
    /* L110: */
  }
  work[*ni + 2] = -12 / *ho * work[*ni + 1] + 6 / *ho * work[*ni + 2];
  work[*ni + 1] = 6 / *ho * work[*ni + 1];
  work[*ni + *no - 1] = 6 / *ho * temp6 - 12 / *ho * temp3;
  work[*ni + *no] = 6 / *ho * temp3;
  work[*ni + *no + 2] = -12 / *ho * work[*ni + *no + 1] + 6 / *ho * work[*
									 ni + *no + 2];
  work[*ni + *no + 1] = 6 / *ho * work[*ni + *no + 1];
  work[*ni + (*no << 1) - 1] = 6 / *ho * temp7 - 12 / *ho * temp4;
  work[*ni + (*no << 1)] = 6 / *ho * temp4;
  work[*ni + (*no << 1) + 2] = -12 / *ho * work[*ni + (*no << 1) + 1] + 6 / 
    *ho * work[*ni + (*no << 1) + 2];
  work[*ni + (*no << 1) + 1] = 6 / *ho * work[*ni + (*no << 1) + 1];
  work[*ni + *no * 3 - 1] = 6 / *ho * temp8 - 12 / *ho * temp5;
  work[*ni + *no * 3] = 6 / *ho * temp5;

  /*     now, the object is to determine the remaining two sets of */
  /*     numbers that sum with the three already calculated to */
  /*     produce the output vector. */
  /*     for efficiency, I will 'recycle' the part of the workspace */
  /*     used in calculating the previous three sets that is */
  /*     no longer needed. */
  /*     Thus, I need to zero out this part of the work vector. */

  i__1 = *ni + *no * 7 - 8;
  for (i__ = *ni + *no * 3 + 1; i__ <= i__1; ++i__) {
    work[i__] = 0.f;
    /* L115: */
  }

  /*     the two sets of numbers to be calculated are stored in */
  /*     work(ni+3*no+1..ni+4*no)(d) and work(ni+4*no..ni+5*no)(e). */
  /*     the components of each set should equal the sum of the */
  /*     elements in the input vector in a specific interval of the */
  /*     output vector (d) or this amount multiplied by dx(e). */

  i__1 = *ni;
  for (i__ = 1; i__ <= i__1; ++i__) {
      j = (int) work[i__];
    work[*ni + *no * 3 + j] += vi[i__];
    dx = *x1i + (i__ - 1) * *hi - (*x1o + (j - 1) * *ho);
    work[*ni + (*no << 2) + j] += vi[i__] * dx;
    /* L120: */
  }

  /*     the last set (e) is also multiplied by a matrix, as reflected */
  /*     in the following calculations */

  for (i__ = *no - 1; i__ >= 2; --i__) {
    work[*ni + (*no << 2) + i__] = work[*ni + (*no << 2) + i__ - 1] / *ho 
      - work[*ni + (*no << 2) + i__] / *ho;
    /* L130: */
  }
  work[*ni + (*no << 2) + 1] = -1 / *ho * work[*ni + (*no << 2) + 1];
  work[*ni + *no * 5] = temp2 / *ho;

  /*     now all five sets of numbers can be summed to produce my */
  /*     output vector. */

  i__1 = *no;
  for (i__ = 1; i__ <= i__1; ++i__) {
    vo[i__] = work[*ni + *no * 3 + i__] + work[*ni + i__] + work[*ni + *
								 no + i__] - work[*ni + (*no << 1) + i__] + work[*ni + (*no << 
															2) + i__];

    /*     the next multiplication translates the transpose to being */
    /*     the adjoint */

    vo[i__] *= *hi / *ho;
    /* L140: */
  }
  return ierr;
} /* cubicadj_ */

