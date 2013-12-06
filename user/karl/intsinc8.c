/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
INTSINC8 - Functions to interpolate uniformly-sampled data via 8-coeff. sinc
		approximations:

ints8c	interpolation of a uniformly-sampled complex function y(x) via an
	 8-coefficient sinc approximation.
ints8r	Interpolation of a uniformly-sampled real function y(x) via a
		table of 8-coefficient sinc approximations

******************************************************************************
Function Prototypes:
void ints8c (int nxin, float dxin, float fxin, complex yin[], 
	sf_complex yinl, sf_complex yinr, int nxout, float xout[], sf_complex yout[]);
void ints8r (int nxin, float dxin, float fxin, float yin[], 
	float yinl, float yinr, int nxout, float xout[], float yout[]);

******************************************************************************
Input:
nxin		number of x values at which y(x) is input
dxin		x sampling interval for input y(x)
fxin		x value of first sample input
yin		array[nxin] of input y(x) values:  yin[0] = y(fxin), etc.
yinl		value used to extrapolate yin values to left of yin[0]
yinr		value used to extrapolate yin values to right of yin[nxin-1]
nxout		number of x values a which y(x) is output
xout		array[nxout] of x values at which y(x) is output

Output:
yout		array[nxout] of output y(x):  yout[0] = y(xout[0]), etc.

******************************************************************************
Notes:
Because extrapolation of the input function y(x) is defined by the
left and right values yinl and yinr, the xout values are not restricted
to lie within the range of sample locations defined by nxin, dxin, and
fxin.

The maximum error for frequiencies less than 0.6 nyquist is less than
one percent.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
/**************** end self doc ********************************/

#include <rsf.h>
#include "intsinc8.h"
#include "mksinc.h"
#include "inttable8.h"


/* these are used by both ints8c and ints8r */
#define LTABLE 8
#define NTABLE 513

void ints8r (int nxin, float dxin, float fxin, float yin[], 
	float yinl, float yinr, int nxout, float xout[], float yout[])
/*< real interp using table of 8 point filters >*/
/*****************************************************************************
Interpolation of a uniformly-sampled real function y(x) via a
table of 8-coefficient sinc approximations; maximum error for frequiencies
less than 0.6 nyquist is less than one percent.
******************************************************************************
Input:
nxin		number of x values at which y(x) is input
dxin		x sampling interval for input y(x)
fxin		x value of first sample input
yin		array[nxin] of input y(x) values:  yin[0] = y(fxin), etc.
yinl		value used to extrapolate yin values to left of yin[0]
yinr		value used to extrapolate yin values to right of yin[nxin-1]
nxout		number of x values a which y(x) is output
xout		array[nxout] of x values at which y(x) is output

Output:
yout		array[nxout] of output y(x):  yout[0] = y(xout[0]), etc.
******************************************************************************
Notes:
Because extrapolation of the input function y(x) is defined by the
left and right values yinl and yinr, the xout values are not restricted
to lie within the range of sample locations defined by nxin, dxin, and
fxin.
******************************************************************************
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
