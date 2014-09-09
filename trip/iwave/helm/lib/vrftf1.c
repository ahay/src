/* vrftf1.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int vrftf1_(integer *m, integer *n, float *c__, float *ch, 
	integer *mdimc, float *wa, float *fac)
{
    /* System generated locals */
    integer ch_dim1, ch_offset, c_dim1, c_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static integer i__, j, k1, l1, l2, na, kh, nf, ip, iw, ix2, ix3, ix4, ido,
	     idl1;
    static float scale;
    extern /* Subroutine */ int vradf2_(integer *, integer *, integer *, float 
	    *, float *, integer *, float *), vradf3_(integer *, integer *, 
	    integer *, float *, float *, integer *, float *, float *), vradf4_(
	    integer *, integer *, integer *, float *, float *, integer *, float *
	    , float *, float *), vradf5_(integer *, integer *, integer *, float *
	    , float *, integer *, float *, float *, float *, float *), vradfg_(
	    integer *, integer *, integer *, integer *, integer *, float *, 
	    float *, float *, float *, float *, integer *, float *);


/*     VRFFTPK, VERSION 1, AUGUST 1985 */

    /* Parameter adjustments */
    --wa;
    ch_dim1 = *mdimc;
    ch_offset = 1 + ch_dim1;
    ch -= ch_offset;
    c_dim1 = *mdimc;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    --fac;

    /* Function Body */
    nf = fac[2];
    na = 1;
    l2 = *n;
    iw = *n;
    i__1 = nf;
    for (k1 = 1; k1 <= i__1; ++k1) {
	kh = nf - k1;
	ip = fac[kh + 3];
	l1 = l2 / ip;
	ido = *n / l2;
	idl1 = ido * l1;
	iw -= (ip - 1) * ido;
	na = 1 - na;
	if (ip != 4) {
	    goto L102;
	}
	ix2 = iw + ido;
	ix3 = ix2 + ido;
	if (na != 0) {
	    goto L101;
	}
	vradf4_(m, &ido, &l1, &c__[c_offset], &ch[ch_offset], mdimc, &wa[iw], 
		&wa[ix2], &wa[ix3]);
	goto L110;
L101:
	vradf4_(m, &ido, &l1, &ch[ch_offset], &c__[c_offset], mdimc, &wa[iw], 
		&wa[ix2], &wa[ix3]);
	goto L110;
L102:
	if (ip != 2) {
	    goto L104;
	}
	if (na != 0) {
	    goto L103;
	}
	vradf2_(m, &ido, &l1, &c__[c_offset], &ch[ch_offset], mdimc, &wa[iw]);
	goto L110;
L103:
	vradf2_(m, &ido, &l1, &ch[ch_offset], &c__[c_offset], mdimc, &wa[iw]);
	goto L110;
L104:
	if (ip != 3) {
	    goto L106;
	}
	ix2 = iw + ido;
	if (na != 0) {
	    goto L105;
	}
	vradf3_(m, &ido, &l1, &c__[c_offset], &ch[ch_offset], mdimc, &wa[iw], 
		&wa[ix2]);
	goto L110;
L105:
	vradf3_(m, &ido, &l1, &ch[ch_offset], &c__[c_offset], mdimc, &wa[iw], 
		&wa[ix2]);
	goto L110;
L106:
	if (ip != 5) {
	    goto L108;
	}
	ix2 = iw + ido;
	ix3 = ix2 + ido;
	ix4 = ix3 + ido;
	if (na != 0) {
	    goto L107;
	}
	vradf5_(m, &ido, &l1, &c__[c_offset], &ch[ch_offset], mdimc, &wa[iw], 
		&wa[ix2], &wa[ix3], &wa[ix4]);
	goto L110;
L107:
	vradf5_(m, &ido, &l1, &ch[ch_offset], &c__[c_offset], mdimc, &wa[iw], 
		&wa[ix2], &wa[ix3], &wa[ix4]);
	goto L110;
L108:
	if (ido == 1) {
	    na = 1 - na;
	}
	if (na != 0) {
	    goto L109;
	}
	vradfg_(m, &ido, &ip, &l1, &idl1, &c__[c_offset], &c__[c_offset], &
		c__[c_offset], &ch[ch_offset], &ch[ch_offset], mdimc, &wa[iw])
		;
	na = 1;
	goto L110;
L109:
	vradfg_(m, &ido, &ip, &l1, &idl1, &ch[ch_offset], &ch[ch_offset], &ch[
		ch_offset], &c__[c_offset], &c__[c_offset], mdimc, &wa[iw]);
	na = 0;
L110:
	l2 = l1;
/* L111: */
    }
    scale = sqrt(1.f / *n);
    if (na == 1) {
	goto L113;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    c__[i__ + j * c_dim1] = scale * ch[i__ + j * ch_dim1];
/* L112: */
	}
    }
    return 0;
L113:
    i__2 = *n;
    for (j = 1; j <= i__2; ++j) {
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    c__[i__ + j * c_dim1] = scale * c__[i__ + j * c_dim1];
/* L114: */
	}
    }
    return 0;
} /* vrftf1_ */

