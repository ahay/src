/* vradf5.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int vradf5_(integer *mp, integer *ido, integer *l1, float *cc,
	 float *ch, integer *mdimc, float *wa1, float *wa2, float *wa3, float *wa4)
{
    /* System generated locals */
    integer cc_dim1, cc_dim2, cc_dim3, cc_offset, ch_dim1, ch_dim2, ch_offset,
	     i__1, i__2, i__3;

    /* Builtin functions */
    double cos(double), sin(double);

    /* Local variables */
    static integer i__, k, m, ic;
    static float arg, ti11, ti12, tr11, tr12;
    static integer idp2;
    extern double pimach_(float *);


/*     VRFFTPK, VERSION 1, AUGUST 1985 */

    /* Parameter adjustments */
    --wa4;
    --wa3;
    --wa2;
    --wa1;
    ch_dim1 = *mdimc;
    ch_dim2 = *ido;
    ch_offset = 1 + ch_dim1 * (1 + ch_dim2 * 6);
    ch -= ch_offset;
    cc_dim1 = *mdimc;
    cc_dim2 = *ido;
    cc_dim3 = *l1;
    cc_offset = 1 + cc_dim1 * (1 + cc_dim2 * (1 + cc_dim3));
    cc -= cc_offset;

    /* Function Body */
    arg = pimach_(&c_b2) * 2.f / 5.f;
    tr11 = cos(arg);
    ti11 = sin(arg);
    tr12 = cos(arg * 2.f);
    ti12 = sin(arg * 2.f);
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *mp;
	for (m = 1; m <= i__2; ++m) {
	    ch[m + ((k * 5 + 1) * ch_dim2 + 1) * ch_dim1] = cc[m + ((k + 
		    cc_dim3) * cc_dim2 + 1) * cc_dim1] + (cc[m + ((k + 
		    cc_dim3 * 5) * cc_dim2 + 1) * cc_dim1] + cc[m + ((k + (
		    cc_dim3 << 1)) * cc_dim2 + 1) * cc_dim1]) + (cc[m + ((k + 
		    (cc_dim3 << 2)) * cc_dim2 + 1) * cc_dim1] + cc[m + ((k + 
		    cc_dim3 * 3) * cc_dim2 + 1) * cc_dim1]);
	    ch[m + (*ido + (k * 5 + 2) * ch_dim2) * ch_dim1] = cc[m + ((k + 
		    cc_dim3) * cc_dim2 + 1) * cc_dim1] + tr11 * (cc[m + ((k + 
		    cc_dim3 * 5) * cc_dim2 + 1) * cc_dim1] + cc[m + ((k + (
		    cc_dim3 << 1)) * cc_dim2 + 1) * cc_dim1]) + tr12 * (cc[m 
		    + ((k + (cc_dim3 << 2)) * cc_dim2 + 1) * cc_dim1] + cc[m 
		    + ((k + cc_dim3 * 3) * cc_dim2 + 1) * cc_dim1]);
	    ch[m + ((k * 5 + 3) * ch_dim2 + 1) * ch_dim1] = ti11 * (cc[m + ((
		    k + cc_dim3 * 5) * cc_dim2 + 1) * cc_dim1] - cc[m + ((k + 
		    (cc_dim3 << 1)) * cc_dim2 + 1) * cc_dim1]) + ti12 * (cc[m 
		    + ((k + (cc_dim3 << 2)) * cc_dim2 + 1) * cc_dim1] - cc[m 
		    + ((k + cc_dim3 * 3) * cc_dim2 + 1) * cc_dim1]);
	    ch[m + (*ido + (k * 5 + 4) * ch_dim2) * ch_dim1] = cc[m + ((k + 
		    cc_dim3) * cc_dim2 + 1) * cc_dim1] + tr12 * (cc[m + ((k + 
		    cc_dim3 * 5) * cc_dim2 + 1) * cc_dim1] + cc[m + ((k + (
		    cc_dim3 << 1)) * cc_dim2 + 1) * cc_dim1]) + tr11 * (cc[m 
		    + ((k + (cc_dim3 << 2)) * cc_dim2 + 1) * cc_dim1] + cc[m 
		    + ((k + cc_dim3 * 3) * cc_dim2 + 1) * cc_dim1]);
	    ch[m + ((k * 5 + 5) * ch_dim2 + 1) * ch_dim1] = ti12 * (cc[m + ((
		    k + cc_dim3 * 5) * cc_dim2 + 1) * cc_dim1] - cc[m + ((k + 
		    (cc_dim3 << 1)) * cc_dim2 + 1) * cc_dim1]) - ti11 * (cc[m 
		    + ((k + (cc_dim3 << 2)) * cc_dim2 + 1) * cc_dim1] - cc[m 
		    + ((k + cc_dim3 * 3) * cc_dim2 + 1) * cc_dim1]);
/* L1001: */
	}
/* L101: */
    }
    if (*ido == 1) {
	return 0;
    }
    idp2 = *ido + 2;
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch[m + (i__ - 1 + (k * 5 + 1) * ch_dim2) * ch_dim1] = cc[m + (
			i__ - 1 + (k + cc_dim3) * cc_dim2) * cc_dim1] + (wa1[
			i__ - 2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) * 
			cc_dim2) * cc_dim1] + wa1[i__ - 1] * cc[m + (i__ + (k 
			+ (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa4[i__ - 
			2] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * cc_dim2) * 
			cc_dim1] + wa4[i__ - 1] * cc[m + (i__ + (k + cc_dim3 *
			 5) * cc_dim2) * cc_dim1])) + (wa2[i__ - 2] * cc[m + (
			i__ - 1 + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] + 
			wa2[i__ - 1] * cc[m + (i__ + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] + (wa3[i__ - 2] * cc[m + (i__ - 1 
			+ (k + (cc_dim3 << 2)) * cc_dim2) * cc_dim1] + wa3[
			i__ - 1] * cc[m + (i__ + (k + (cc_dim3 << 2)) * 
			cc_dim2) * cc_dim1]));
		ch[m + (i__ + (k * 5 + 1) * ch_dim2) * ch_dim1] = cc[m + (i__ 
			+ (k + cc_dim3) * cc_dim2) * cc_dim1] + (wa1[i__ - 2] 
			* cc[m + (i__ + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa4[i__ - 2] *
			 cc[m + (i__ + (k + cc_dim3 * 5) * cc_dim2) * cc_dim1]
			 - wa4[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) 
			* cc_dim2) * cc_dim1])) + (wa2[i__ - 2] * cc[m + (i__ 
			+ (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] - wa2[i__ - 
			1] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * cc_dim2) * 
			cc_dim1] + (wa3[i__ - 2] * cc[m + (i__ + (k + (
			cc_dim3 << 2)) * cc_dim2) * cc_dim1] - wa3[i__ - 1] * 
			cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1]));
		ch[m + (i__ - 1 + (k * 5 + 3) * ch_dim2) * ch_dim1] = cc[m + (
			i__ - 1 + (k + cc_dim3) * cc_dim2) * cc_dim1] + tr11 *
			 (wa1[i__ - 2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)
			) * cc_dim2) * cc_dim1] + wa1[i__ - 1] * cc[m + (i__ 
			+ (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + wa4[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * 
			cc_dim2) * cc_dim1] + wa4[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 5) * cc_dim2) * cc_dim1]) + tr12 * (wa2[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] + wa2[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 3) * cc_dim2) * cc_dim1] + wa3[i__ - 2] * 
			cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] + wa3[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 2)) * cc_dim2) * cc_dim1]) + ti11 * (wa1[i__ - 2] *
			 cc[m + (i__ + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 1)) * cc_dim2) * cc_dim1] - (wa4[i__ - 2] *
			 cc[m + (i__ + (k + cc_dim3 * 5) * cc_dim2) * cc_dim1]
			 - wa4[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) 
			* cc_dim2) * cc_dim1])) + ti12 * (wa2[i__ - 2] * cc[m 
			+ (i__ + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] - 
			wa2[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] - (wa3[i__ - 2] * cc[m + (i__ + (
			k + (cc_dim3 << 2)) * cc_dim2) * cc_dim1] - wa3[i__ - 
			1] * cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2)
			 * cc_dim1]));
		ch[m + (ic - 1 + (k * 5 + 2) * ch_dim2) * ch_dim1] = cc[m + (
			i__ - 1 + (k + cc_dim3) * cc_dim2) * cc_dim1] + tr11 *
			 (wa1[i__ - 2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)
			) * cc_dim2) * cc_dim1] + wa1[i__ - 1] * cc[m + (i__ 
			+ (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + wa4[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * 
			cc_dim2) * cc_dim1] + wa4[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 5) * cc_dim2) * cc_dim1]) + tr12 * (wa2[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] + wa2[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 3) * cc_dim2) * cc_dim1] + wa3[i__ - 2] * 
			cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] + wa3[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 2)) * cc_dim2) * cc_dim1]) - (ti11 * (wa1[i__ - 2] 
			* cc[m + (i__ + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 1)) * cc_dim2) * cc_dim1] - (wa4[i__ - 2] *
			 cc[m + (i__ + (k + cc_dim3 * 5) * cc_dim2) * cc_dim1]
			 - wa4[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) 
			* cc_dim2) * cc_dim1])) + ti12 * (wa2[i__ - 2] * cc[m 
			+ (i__ + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] - 
			wa2[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] - (wa3[i__ - 2] * cc[m + (i__ + (
			k + (cc_dim3 << 2)) * cc_dim2) * cc_dim1] - wa3[i__ - 
			1] * cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2)
			 * cc_dim1])));
		ch[m + (i__ + (k * 5 + 3) * ch_dim2) * ch_dim1] = cc[m + (i__ 
			+ (k + cc_dim3) * cc_dim2) * cc_dim1] + tr11 * (wa1[
			i__ - 2] * cc[m + (i__ + (k + (cc_dim3 << 1)) * 
			cc_dim2) * cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 
			+ (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa4[
			i__ - 2] * cc[m + (i__ + (k + cc_dim3 * 5) * cc_dim2) 
			* cc_dim1] - wa4[i__ - 1] * cc[m + (i__ - 1 + (k + 
			cc_dim3 * 5) * cc_dim2) * cc_dim1])) + tr12 * (wa2[
			i__ - 2] * cc[m + (i__ + (k + cc_dim3 * 3) * cc_dim2) 
			* cc_dim1] - wa2[i__ - 1] * cc[m + (i__ - 1 + (k + 
			cc_dim3 * 3) * cc_dim2) * cc_dim1] + (wa3[i__ - 2] * 
			cc[m + (i__ + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] - wa3[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 2)) * cc_dim2) * cc_dim1])) + (ti11 * (wa4[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * 
			cc_dim2) * cc_dim1] + wa4[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 5) * cc_dim2) * cc_dim1] - (wa1[i__ - 2] *
			 cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] + wa1[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 1)) * cc_dim2) * cc_dim1])) + ti12 * (wa3[i__ - 2] 
			* cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] + wa3[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 2)) * cc_dim2) * cc_dim1] - (wa2[i__ - 2] * cc[m + 
			(i__ - 1 + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] + 
			wa2[i__ - 1] * cc[m + (i__ + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1])));
		ch[m + (ic + (k * 5 + 2) * ch_dim2) * ch_dim1] = ti11 * (wa4[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * 
			cc_dim2) * cc_dim1] + wa4[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 5) * cc_dim2) * cc_dim1] - (wa1[i__ - 2] *
			 cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] + wa1[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 1)) * cc_dim2) * cc_dim1])) + ti12 * (wa3[i__ - 2] 
			* cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] + wa3[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 2)) * cc_dim2) * cc_dim1] - (wa2[i__ - 2] * cc[m + 
			(i__ - 1 + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] + 
			wa2[i__ - 1] * cc[m + (i__ + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1])) - (cc[m + (i__ + (k + cc_dim3) *
			 cc_dim2) * cc_dim1] + tr11 * (wa1[i__ - 2] * cc[m + (
			i__ + (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] - 
			wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) 
			* cc_dim2) * cc_dim1] + (wa4[i__ - 2] * cc[m + (i__ + 
			(k + cc_dim3 * 5) * cc_dim2) * cc_dim1] - wa4[i__ - 1]
			 * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * cc_dim2) * 
			cc_dim1])) + tr12 * (wa2[i__ - 2] * cc[m + (i__ + (k 
			+ cc_dim3 * 3) * cc_dim2) * cc_dim1] - wa2[i__ - 1] * 
			cc[m + (i__ - 1 + (k + cc_dim3 * 3) * cc_dim2) * 
			cc_dim1] + (wa3[i__ - 2] * cc[m + (i__ + (k + (
			cc_dim3 << 2)) * cc_dim2) * cc_dim1] - wa3[i__ - 1] * 
			cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1])));
		ch[m + (i__ - 1 + (k * 5 + 5) * ch_dim2) * ch_dim1] = cc[m + (
			i__ - 1 + (k + cc_dim3) * cc_dim2) * cc_dim1] + tr12 *
			 (wa1[i__ - 2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)
			) * cc_dim2) * cc_dim1] + wa1[i__ - 1] * cc[m + (i__ 
			+ (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa4[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * 
			cc_dim2) * cc_dim1] + wa4[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 5) * cc_dim2) * cc_dim1])) + tr11 * (wa2[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] + wa2[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 3) * cc_dim2) * cc_dim1] + (wa3[i__ - 2] *
			 cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] + wa3[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 2)) * cc_dim2) * cc_dim1])) + (ti12 * (wa1[i__ - 2]
			 * cc[m + (i__ + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 1)) * cc_dim2) * cc_dim1] - (wa4[i__ - 2] *
			 cc[m + (i__ + (k + cc_dim3 * 5) * cc_dim2) * cc_dim1]
			 - wa4[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) 
			* cc_dim2) * cc_dim1])) - ti11 * (wa2[i__ - 2] * cc[m 
			+ (i__ + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] - 
			wa2[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] - (wa3[i__ - 2] * cc[m + (i__ + (
			k + (cc_dim3 << 2)) * cc_dim2) * cc_dim1] - wa3[i__ - 
			1] * cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2)
			 * cc_dim1])));
		ch[m + (ic - 1 + (k * 5 + 4) * ch_dim2) * ch_dim1] = cc[m + (
			i__ - 1 + (k + cc_dim3) * cc_dim2) * cc_dim1] + tr12 *
			 (wa1[i__ - 2] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)
			) * cc_dim2) * cc_dim1] + wa1[i__ - 1] * cc[m + (i__ 
			+ (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa4[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * 
			cc_dim2) * cc_dim1] + wa4[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 5) * cc_dim2) * cc_dim1])) + tr11 * (wa2[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] + wa2[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 3) * cc_dim2) * cc_dim1] + (wa3[i__ - 2] *
			 cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] + wa3[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 2)) * cc_dim2) * cc_dim1])) - (ti12 * (wa1[i__ - 2]
			 * cc[m + (i__ + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 1)) * cc_dim2) * cc_dim1] - (wa4[i__ - 2] *
			 cc[m + (i__ + (k + cc_dim3 * 5) * cc_dim2) * cc_dim1]
			 - wa4[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) 
			* cc_dim2) * cc_dim1])) - ti11 * (wa2[i__ - 2] * cc[m 
			+ (i__ + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] - 
			wa2[i__ - 1] * cc[m + (i__ - 1 + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1] - (wa3[i__ - 2] * cc[m + (i__ + (
			k + (cc_dim3 << 2)) * cc_dim2) * cc_dim1] - wa3[i__ - 
			1] * cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2)
			 * cc_dim1])));
		ch[m + (i__ + (k * 5 + 5) * ch_dim2) * ch_dim1] = cc[m + (i__ 
			+ (k + cc_dim3) * cc_dim2) * cc_dim1] + tr12 * (wa1[
			i__ - 2] * cc[m + (i__ + (k + (cc_dim3 << 1)) * 
			cc_dim2) * cc_dim1] - wa1[i__ - 1] * cc[m + (i__ - 1 
			+ (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] + (wa4[
			i__ - 2] * cc[m + (i__ + (k + cc_dim3 * 5) * cc_dim2) 
			* cc_dim1] - wa4[i__ - 1] * cc[m + (i__ - 1 + (k + 
			cc_dim3 * 5) * cc_dim2) * cc_dim1])) + tr11 * (wa2[
			i__ - 2] * cc[m + (i__ + (k + cc_dim3 * 3) * cc_dim2) 
			* cc_dim1] - wa2[i__ - 1] * cc[m + (i__ - 1 + (k + 
			cc_dim3 * 3) * cc_dim2) * cc_dim1] + (wa3[i__ - 2] * 
			cc[m + (i__ + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] - wa3[i__ - 1] * cc[m + (i__ - 1 + (k + (
			cc_dim3 << 2)) * cc_dim2) * cc_dim1])) + (ti12 * (wa4[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * 
			cc_dim2) * cc_dim1] + wa4[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 5) * cc_dim2) * cc_dim1] - (wa1[i__ - 2] *
			 cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] + wa1[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 1)) * cc_dim2) * cc_dim1])) - ti11 * (wa3[i__ - 2] 
			* cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] + wa3[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 2)) * cc_dim2) * cc_dim1] - (wa2[i__ - 2] * cc[m + 
			(i__ - 1 + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] + 
			wa2[i__ - 1] * cc[m + (i__ + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1])));
		ch[m + (ic + (k * 5 + 4) * ch_dim2) * ch_dim1] = ti12 * (wa4[
			i__ - 2] * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * 
			cc_dim2) * cc_dim1] + wa4[i__ - 1] * cc[m + (i__ + (k 
			+ cc_dim3 * 5) * cc_dim2) * cc_dim1] - (wa1[i__ - 2] *
			 cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) * cc_dim2) * 
			cc_dim1] + wa1[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 1)) * cc_dim2) * cc_dim1])) - ti11 * (wa3[i__ - 2] 
			* cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1] + wa3[i__ - 1] * cc[m + (i__ + (k + (cc_dim3 
			<< 2)) * cc_dim2) * cc_dim1] - (wa2[i__ - 2] * cc[m + 
			(i__ - 1 + (k + cc_dim3 * 3) * cc_dim2) * cc_dim1] + 
			wa2[i__ - 1] * cc[m + (i__ + (k + cc_dim3 * 3) * 
			cc_dim2) * cc_dim1])) - (cc[m + (i__ + (k + cc_dim3) *
			 cc_dim2) * cc_dim1] + tr12 * (wa1[i__ - 2] * cc[m + (
			i__ + (k + (cc_dim3 << 1)) * cc_dim2) * cc_dim1] - 
			wa1[i__ - 1] * cc[m + (i__ - 1 + (k + (cc_dim3 << 1)) 
			* cc_dim2) * cc_dim1] + (wa4[i__ - 2] * cc[m + (i__ + 
			(k + cc_dim3 * 5) * cc_dim2) * cc_dim1] - wa4[i__ - 1]
			 * cc[m + (i__ - 1 + (k + cc_dim3 * 5) * cc_dim2) * 
			cc_dim1])) + tr11 * (wa2[i__ - 2] * cc[m + (i__ + (k 
			+ cc_dim3 * 3) * cc_dim2) * cc_dim1] - wa2[i__ - 1] * 
			cc[m + (i__ - 1 + (k + cc_dim3 * 3) * cc_dim2) * 
			cc_dim1] + (wa3[i__ - 2] * cc[m + (i__ + (k + (
			cc_dim3 << 2)) * cc_dim2) * cc_dim1] - wa3[i__ - 1] * 
			cc[m + (i__ - 1 + (k + (cc_dim3 << 2)) * cc_dim2) * 
			cc_dim1])));
/* L1002: */
	    }
/* L102: */
	}
/* L103: */
    }
    return 0;
} /* vradf5_ */

