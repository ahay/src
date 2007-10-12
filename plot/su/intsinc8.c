/* Copyright (c) Colorado School of Mines, 2007.*/
/* All rights reserved.                       */

#include <rsf.h>

/* these are used by both ints8c and ints8r */
#define LTABLE 8
#define NTABLE 513

#include "mksinc.h"
#include "inttable8.h"

void ints8c (int nxin         /* number of x values at which y(x) is input */, 
	     float dxin       /* x sampling interval for input y(x) */, 
	     float fxin       /* x value of first sample input */, 
	     sf_complex *yin  /* [nxin] array of input y(x) values:  yin[0] = y(fxin), etc. */,
	     sf_complex yinl  /* value used to extrapolate yin values to left of yin[0] */, 
	     sf_complex yinr  /* value used to extrapolate yin values to right of yin[nxin-1] */, 
	     int nxout        /* number of x values a which y(x) is output */, 
	     float *xout      /* [nxout] array of x values at which y(x) is output */, 
	     sf_complex *yout /* [nxout] array of output y(x):  yout[0] = y(xout[0]), etc. */)
/*< Interpolation of a uniformly-sampled complex function y(x) via a
table of 8-coefficient sinc approximations; maximum error for frequiencies
less than 0.6 nyquist is less than one percent.

Because extrapolation of the input function y(x) is defined by the
left and right values yinl and yinr, the xout values are not restricted
to lie within the range of sample locations defined by nxin, dxin, and
fxin. >*/
/******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
        static float table[NTABLE][LTABLE];
	static int tabled=0;
	int jtable;
	float frac;

	/* tabulate sinc interpolation coefficients if not already tabulated */
	if (!tabled) {
		for (jtable=1; jtable<NTABLE-1; jtable++) {
			frac = (float)jtable/(float)(NTABLE-1);
			mksinc(frac,LTABLE,&table[jtable][0]);
		}
		for (jtable=0; jtable<LTABLE; jtable++) {
			table[0][jtable] = 0.0;
			table[NTABLE-1][jtable] = 0.0;
		}
		table[0][LTABLE/2-1] = 1.0;
		table[NTABLE-1][LTABLE/2] = 1.0;
		tabled = 1;
	}

	/* interpolate using tabulated coefficients */
	intt8c(NTABLE,table,nxin,dxin,fxin,yin,yinl,yinr,nxout,xout,yout);
}

void ints8r (int nxin /* number of x values at which y(x) is input */, 
	     float dxin /* x sampling interval for input y(x) */, 
	     float fxin /* x value of first sample input */, 
	     float *yin /* [nxin] array of input y(x) values:  yin[0] = y(fxin), etc. */, 
	     float yinl  /* value used to extrapolate yin values to left of yin[0] */, 
	     float yinr  /* alue used to extrapolate yin values to right of yin[nxin-1] */, 
	     int nxout   /* number of x values a which y(x) is output */, 
	     float *xout /* [nxout] array of x values at which y(x) is output */, 
	     float *yout /* [nxout] array of output y(x):  yout[0] = y(xout[0]), etc. */)
/*< Interpolation of a uniformly-sampled real function y(x) via a
table of 8-coefficient sinc approximations; maximum error for frequiencies
less than 0.6 nyquist is less than one percent.

Because extrapolation of the input function y(x) is defined by the
left and right values yinl and yinr, the xout values are not restricted
to lie within the range of sample locations defined by nxin, dxin, and
fxin. >*/
/******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	static float table[NTABLE][LTABLE];
	static int tabled=0;
	int jtable;
	float frac;

	/* tabulate sinc interpolation coefficients if not already tabulated */
	if (!tabled) {
		for (jtable=1; jtable<NTABLE-1; jtable++) {
			frac = (float)jtable/(float)(NTABLE-1);
			mksinc(frac,LTABLE,&table[jtable][0]);
		}
		for (jtable=0; jtable<LTABLE; jtable++) {
			table[0][jtable] = 0.0;
			table[NTABLE-1][jtable] = 0.0;
		}
		table[0][LTABLE/2-1] = 1.0;
		table[NTABLE-1][LTABLE/2] = 1.0;
		tabled = 1;
	}

	/* interpolate using tabulated coefficients */
	intt8r(NTABLE,table,nxin,dxin,fxin,yin,yinl,yinr,nxout,xout,yout);
}
