/* vradfg.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int vradfg_(integer *mp, integer *ido, integer *ip, integer *
	l1, integer *idl1, float *cc, float *c1, float *c2, float *ch, float *ch2, 
	integer *mdimc, float *wa)
{
    /* System generated locals */
    integer ch_dim1, ch_dim2, ch_dim3, ch_offset, cc_dim1, cc_dim2, cc_dim3, 
	    cc_offset, c1_dim1, c1_dim2, c1_dim3, c1_offset, c2_dim1, c2_dim2,
	     c2_offset, ch2_dim1, ch2_dim2, ch2_offset, i__1, i__2, i__3, 
	    i__4;

    /* Builtin functions */
    double cos(double), sin(double);

    /* Local variables */
    static integer i__, j, k, l, m, j2, ic, jc, lc, ik, is;
    static float dc2, ai1, ai2, ar1, ar2, ds2;
    static integer nbd;
    static float dcp, arg, dsp, tpi, ar1h, ar2h;
    static integer idp2, ipp2, idij, ipph;
    extern double pimach_(float *);


/*     VRFFTPK, VERSION 1, AUGUST 1985 */

    /* Parameter adjustments */
    --wa;
    ch2_dim1 = *mdimc;
    ch2_dim2 = *idl1;
    ch2_offset = 1 + ch2_dim1 * (1 + ch2_dim2);
    ch2 -= ch2_offset;
    ch_dim1 = *mdimc;
    ch_dim2 = *ido;
    ch_dim3 = *l1;
    ch_offset = 1 + ch_dim1 * (1 + ch_dim2 * (1 + ch_dim3));
    ch -= ch_offset;
    c2_dim1 = *mdimc;
    c2_dim2 = *idl1;
    c2_offset = 1 + c2_dim1 * (1 + c2_dim2);
    c2 -= c2_offset;
    c1_dim1 = *mdimc;
    c1_dim2 = *ido;
    c1_dim3 = *l1;
    c1_offset = 1 + c1_dim1 * (1 + c1_dim2 * (1 + c1_dim3));
    c1 -= c1_offset;
    cc_dim1 = *mdimc;
    cc_dim2 = *ido;
    cc_dim3 = *ip;
    cc_offset = 1 + cc_dim1 * (1 + cc_dim2 * (1 + cc_dim3));
    cc -= cc_offset;

    /* Function Body */
    tpi = pimach_(&c_b2) * 2.f;
    arg = tpi / (float) (*ip);
    dcp = cos(arg);
    dsp = sin(arg);
    ipph = (*ip + 1) / 2;
    ipp2 = *ip + 2;
    idp2 = *ido + 2;
    nbd = (*ido - 1) / 2;
    if (*ido == 1) {
	goto L119;
    }
    i__1 = *idl1;
    for (ik = 1; ik <= i__1; ++ik) {
	i__2 = *mp;
	for (m = 1; m <= i__2; ++m) {
	    ch2[m + (ik + ch2_dim2) * ch2_dim1] = c2[m + (ik + c2_dim2) * 
		    c2_dim1];
/* L1001: */
	}
/* L101: */
    }
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch[m + ((k + j * ch_dim3) * ch_dim2 + 1) * ch_dim1] = c1[m + (
			(k + j * c1_dim3) * c1_dim2 + 1) * c1_dim1];
/* L1002: */
	    }
/* L102: */
	}
/* L103: */
    }
    if (nbd > *l1) {
	goto L107;
    }
    is = -(*ido);
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	is += *ido;
	idij = is;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    idij += 2;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    ch[m + (i__ - 1 + (k + j * ch_dim3) * ch_dim2) * ch_dim1] 
			    = wa[idij - 1] * c1[m + (i__ - 1 + (k + j * 
			    c1_dim3) * c1_dim2) * c1_dim1] + wa[idij] * c1[m 
			    + (i__ + (k + j * c1_dim3) * c1_dim2) * c1_dim1];
		    ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * ch_dim1] = 
			    wa[idij - 1] * c1[m + (i__ + (k + j * c1_dim3) * 
			    c1_dim2) * c1_dim1] - wa[idij] * c1[m + (i__ - 1 
			    + (k + j * c1_dim3) * c1_dim2) * c1_dim1];
/* L1004: */
		}
/* L104: */
	    }
/* L105: */
	}
/* L106: */
    }
    goto L111;
L107:
    is = -(*ido);
    i__1 = *ip;
    for (j = 2; j <= i__1; ++j) {
	is += *ido;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    idij = is;
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		idij += 2;
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    ch[m + (i__ - 1 + (k + j * ch_dim3) * ch_dim2) * ch_dim1] 
			    = wa[idij - 1] * c1[m + (i__ - 1 + (k + j * 
			    c1_dim3) * c1_dim2) * c1_dim1] + wa[idij] * c1[m 
			    + (i__ + (k + j * c1_dim3) * c1_dim2) * c1_dim1];
		    ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * ch_dim1] = 
			    wa[idij - 1] * c1[m + (i__ + (k + j * c1_dim3) * 
			    c1_dim2) * c1_dim1] - wa[idij] * c1[m + (i__ - 1 
			    + (k + j * c1_dim3) * c1_dim2) * c1_dim1];
/* L1008: */
		}
/* L108: */
	    }
/* L109: */
	}
/* L110: */
    }
L111:
    if (nbd < *l1) {
	goto L115;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    c1[m + (i__ - 1 + (k + j * c1_dim3) * c1_dim2) * c1_dim1] 
			    = ch[m + (i__ - 1 + (k + j * ch_dim3) * ch_dim2) *
			     ch_dim1] + ch[m + (i__ - 1 + (k + jc * ch_dim3) *
			     ch_dim2) * ch_dim1];
		    c1[m + (i__ - 1 + (k + jc * c1_dim3) * c1_dim2) * c1_dim1]
			     = ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * 
			    ch_dim1] - ch[m + (i__ + (k + jc * ch_dim3) * 
			    ch_dim2) * ch_dim1];
		    c1[m + (i__ + (k + j * c1_dim3) * c1_dim2) * c1_dim1] = 
			    ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * 
			    ch_dim1] + ch[m + (i__ + (k + jc * ch_dim3) * 
			    ch_dim2) * ch_dim1];
		    c1[m + (i__ + (k + jc * c1_dim3) * c1_dim2) * c1_dim1] = 
			    ch[m + (i__ - 1 + (k + jc * ch_dim3) * ch_dim2) * 
			    ch_dim1] - ch[m + (i__ - 1 + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1];
/* L1012: */
		}
/* L112: */
	    }
/* L113: */
	}
/* L114: */
    }
    goto L121;
L115:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    c1[m + (i__ - 1 + (k + j * c1_dim3) * c1_dim2) * c1_dim1] 
			    = ch[m + (i__ - 1 + (k + j * ch_dim3) * ch_dim2) *
			     ch_dim1] + ch[m + (i__ - 1 + (k + jc * ch_dim3) *
			     ch_dim2) * ch_dim1];
		    c1[m + (i__ - 1 + (k + jc * c1_dim3) * c1_dim2) * c1_dim1]
			     = ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * 
			    ch_dim1] - ch[m + (i__ + (k + jc * ch_dim3) * 
			    ch_dim2) * ch_dim1];
		    c1[m + (i__ + (k + j * c1_dim3) * c1_dim2) * c1_dim1] = 
			    ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * 
			    ch_dim1] + ch[m + (i__ + (k + jc * ch_dim3) * 
			    ch_dim2) * ch_dim1];
		    c1[m + (i__ + (k + jc * c1_dim3) * c1_dim2) * c1_dim1] = 
			    ch[m + (i__ - 1 + (k + jc * ch_dim3) * ch_dim2) * 
			    ch_dim1] - ch[m + (i__ - 1 + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1];
/* L1016: */
		}
/* L116: */
	    }
/* L117: */
	}
/* L118: */
    }
    goto L121;
L119:
    i__1 = *idl1;
    for (ik = 1; ik <= i__1; ++ik) {
	i__2 = *mp;
	for (m = 1; m <= i__2; ++m) {
	    c2[m + (ik + c2_dim2) * c2_dim1] = ch2[m + (ik + ch2_dim2) * 
		    ch2_dim1];
/* L1020: */
	}
/* L120: */
    }
L121:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		c1[m + ((k + j * c1_dim3) * c1_dim2 + 1) * c1_dim1] = ch[m + (
			(k + j * ch_dim3) * ch_dim2 + 1) * ch_dim1] + ch[m + (
			(k + jc * ch_dim3) * ch_dim2 + 1) * ch_dim1];
		c1[m + ((k + jc * c1_dim3) * c1_dim2 + 1) * c1_dim1] = ch[m + 
			((k + jc * ch_dim3) * ch_dim2 + 1) * ch_dim1] - ch[m 
			+ ((k + j * ch_dim3) * ch_dim2 + 1) * ch_dim1];
/* L1022: */
	    }
/* L122: */
	}
/* L123: */
    }

    ar1 = 1.f;
    ai1 = 0.f;
    i__1 = ipph;
    for (l = 2; l <= i__1; ++l) {
	lc = ipp2 - l;
	ar1h = dcp * ar1 - dsp * ai1;
	ai1 = dcp * ai1 + dsp * ar1;
	ar1 = ar1h;
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch2[m + (ik + l * ch2_dim2) * ch2_dim1] = c2[m + (ik + 
			c2_dim2) * c2_dim1] + ar1 * c2[m + (ik + (c2_dim2 << 
			1)) * c2_dim1];
		ch2[m + (ik + lc * ch2_dim2) * ch2_dim1] = ai1 * c2[m + (ik + 
			*ip * c2_dim2) * c2_dim1];
/* L1024: */
	    }
/* L124: */
	}
	dc2 = ar1;
	ds2 = ai1;
	ar2 = ar1;
	ai2 = ai1;
	i__2 = ipph;
	for (j = 3; j <= i__2; ++j) {
	    jc = ipp2 - j;
	    ar2h = dc2 * ar2 - ds2 * ai2;
	    ai2 = dc2 * ai2 + ds2 * ar2;
	    ar2 = ar2h;
	    i__3 = *idl1;
	    for (ik = 1; ik <= i__3; ++ik) {
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    ch2[m + (ik + l * ch2_dim2) * ch2_dim1] += ar2 * c2[m + (
			    ik + j * c2_dim2) * c2_dim1];
		    ch2[m + (ik + lc * ch2_dim2) * ch2_dim1] += ai2 * c2[m + (
			    ik + jc * c2_dim2) * c2_dim1];
/* L1025: */
		}
/* L125: */
	    }
/* L126: */
	}
/* L127: */
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *idl1;
	for (ik = 1; ik <= i__2; ++ik) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		ch2[m + (ik + ch2_dim2) * ch2_dim1] += c2[m + (ik + j * 
			c2_dim2) * c2_dim1];
/* L1028: */
	    }
/* L128: */
	}
/* L129: */
    }

    if (*ido < *l1) {
	goto L132;
    }
    i__1 = *l1;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *ido;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		cc[m + (i__ + (k * cc_dim3 + 1) * cc_dim2) * cc_dim1] = ch[m 
			+ (i__ + (k + ch_dim3) * ch_dim2) * ch_dim1];
/* L1030: */
	    }
/* L130: */
	}
/* L131: */
    }
    goto L135;
L132:
    i__1 = *ido;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		cc[m + (i__ + (k * cc_dim3 + 1) * cc_dim2) * cc_dim1] = ch[m 
			+ (i__ + (k + ch_dim3) * ch_dim2) * ch_dim1];
/* L1033: */
	    }
/* L133: */
	}
/* L134: */
    }
L135:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *mp;
	    for (m = 1; m <= i__3; ++m) {
		cc[m + (*ido + (j2 - 2 + k * cc_dim3) * cc_dim2) * cc_dim1] = 
			ch[m + ((k + j * ch_dim3) * ch_dim2 + 1) * ch_dim1];
		cc[m + ((j2 - 1 + k * cc_dim3) * cc_dim2 + 1) * cc_dim1] = ch[
			m + ((k + jc * ch_dim3) * ch_dim2 + 1) * ch_dim1];
/* L1036: */
	    }
/* L136: */
	}
/* L137: */
    }
    if (*ido == 1) {
	return 0;
    }
    if (nbd < *l1) {
	goto L141;
    }
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *l1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = *ido;
	    for (i__ = 3; i__ <= i__3; i__ += 2) {
		ic = idp2 - i__;
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    cc[m + (i__ - 1 + (j2 - 1 + k * cc_dim3) * cc_dim2) * 
			    cc_dim1] = ch[m + (i__ - 1 + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1] + ch[m + (i__ - 1 + (k + jc * 
			    ch_dim3) * ch_dim2) * ch_dim1];
		    cc[m + (ic - 1 + (j2 - 2 + k * cc_dim3) * cc_dim2) * 
			    cc_dim1] = ch[m + (i__ - 1 + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1] - ch[m + (i__ - 1 + (k + jc * 
			    ch_dim3) * ch_dim2) * ch_dim1];
		    cc[m + (i__ + (j2 - 1 + k * cc_dim3) * cc_dim2) * cc_dim1]
			     = ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * 
			    ch_dim1] + ch[m + (i__ + (k + jc * ch_dim3) * 
			    ch_dim2) * ch_dim1];
		    cc[m + (ic + (j2 - 2 + k * cc_dim3) * cc_dim2) * cc_dim1] 
			    = ch[m + (i__ + (k + jc * ch_dim3) * ch_dim2) * 
			    ch_dim1] - ch[m + (i__ + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1];
/* L1038: */
		}
/* L138: */
	    }
/* L139: */
	}
/* L140: */
    }
    return 0;
L141:
    i__1 = ipph;
    for (j = 2; j <= i__1; ++j) {
	jc = ipp2 - j;
	j2 = j + j;
	i__2 = *ido;
	for (i__ = 3; i__ <= i__2; i__ += 2) {
	    ic = idp2 - i__;
	    i__3 = *l1;
	    for (k = 1; k <= i__3; ++k) {
		i__4 = *mp;
		for (m = 1; m <= i__4; ++m) {
		    cc[m + (i__ - 1 + (j2 - 1 + k * cc_dim3) * cc_dim2) * 
			    cc_dim1] = ch[m + (i__ - 1 + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1] + ch[m + (i__ - 1 + (k + jc * 
			    ch_dim3) * ch_dim2) * ch_dim1];
		    cc[m + (ic - 1 + (j2 - 2 + k * cc_dim3) * cc_dim2) * 
			    cc_dim1] = ch[m + (i__ - 1 + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1] - ch[m + (i__ - 1 + (k + jc * 
			    ch_dim3) * ch_dim2) * ch_dim1];
		    cc[m + (i__ + (j2 - 1 + k * cc_dim3) * cc_dim2) * cc_dim1]
			     = ch[m + (i__ + (k + j * ch_dim3) * ch_dim2) * 
			    ch_dim1] + ch[m + (i__ + (k + jc * ch_dim3) * 
			    ch_dim2) * ch_dim1];
		    cc[m + (ic + (j2 - 2 + k * cc_dim3) * cc_dim2) * cc_dim1] 
			    = ch[m + (i__ + (k + jc * ch_dim3) * ch_dim2) * 
			    ch_dim1] - ch[m + (i__ + (k + j * ch_dim3) * 
			    ch_dim2) * ch_dim1];
/* L1042: */
		}
/* L142: */
	    }
/* L143: */
	}
/* L144: */
    }
    return 0;
} /* vradfg_ */

