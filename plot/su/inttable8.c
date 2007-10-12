/* Copyright (c) Colorado School of Mines, 2007.*/
/* All rights reserved.                       */

/******************************************************************************
AUTHOR:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
  
#include <rsf.h>

void intt8c (int ntable       /* number of tabulated interpolation operators; ntable>=2 */, 
	     float table[][8] /* array of tabulated 8-point interpolation operators */,
	     int nxin         /* number of x values at which y(x) is input */, 
	     float dxin       /* x sampling interval for input y(x) */, 
	     float fxin       /* x value of first sample input */, 
	     sf_complex *yin  /* array of input y(x) values:  yin[0] = y(fxin), etc. */,
	     sf_complex yinl  /* value used to extrapolate yin values to left of yin[0] */, 
	     sf_complex yinr  /* value used to extrapolate yin values to right of yin[nxin-1] */, 
	     int nxout        /* number of x values a which y(x) is output */, 
	     float *xout      /*number of x values a which y(x) is output */, 
	     sf_complex *yout /* array of output y(x) values:  yout[0] = y(xout[0]), etc. */)
/*< interpolation of a uniformly-sampled complex function y(x)
  via a table of 8-coefficient interpolators 

  The table of interpolation operators must be as follows:

Let d be the distance, expressed as a fraction of dxin, from a particular
xout value to the sampled location xin just to the left of xout.  Then,
for d = 0.0,

table[0][0:7] = 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0

are the weights applied to the 8 input samples nearest xout.
Likewise, for d = 1.0,

table[ntable-1][0:7] = 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0

are the weights applied to the 8 input samples nearest xout.  In general,
for d = (float)itable/(float)(ntable-1), table[itable][0:7] are the
weights applied to the 8 input samples nearest xout.  If the actual sample
distance d does not exactly equal one of the values for which interpolators
are tabulated, then the interpolator corresponding to the nearest value of
d is used.

Because extrapolation of the input function y(x) is defined by the left
and right values yinl and yinr, the xout values are not restricted to lie
within the range of sample locations defined by nxin, dxin, and fxin. >*/

{
	int ioutb,ixout,ixoutn,kyin,ktable,itable,nxinm8;
	float xoutb,xoutf,xouts,xoutn,fntablem1,frac,
		yinlr,yinli,yinrr,yinri,yinir,yinii,sumr,sumi,
		*table00,*pyin,*ptable;
	sf_complex *yin0;

	/* compute constants */
	ioutb = -3-8;
	xoutf = fxin;
	xouts = 1.0/dxin;
	xoutb = 8.0-xoutf*xouts;
	fntablem1 = (float)(ntable-1);
	nxinm8 = nxin-8;
	yin0 = &yin[0];
	table00 = &table[0][0];
	yinlr = crealf(yinl);  yinli = cimagf(yinl);
	yinrr = crealf(yinr);  yinri = cimagf(yinr);

	/* loop over output samples */
	for (ixout=0; ixout<nxout; ixout++) {

		/* determine pointers into table and yin */
		xoutn = xoutb+xout[ixout]*xouts;
		ixoutn = (int)xoutn;
		kyin = ioutb+ixoutn;
		pyin = (float*)(yin0+kyin);
		frac = xoutn-(float)ixoutn;
		ktable = frac>=0.0?frac*fntablem1+0.5:(frac+1.0)*fntablem1-0.5;
		ptable = table00+ktable*8;
		
		/* if totally within input array, use fast method */
		if (kyin>=0 && kyin<=nxinm8) {
		    yout[ixout] = 
			sf_cmplx(pyin[0]*ptable[0]+
				 pyin[2]*ptable[1]+
				 pyin[4]*ptable[2]+
				 pyin[6]*ptable[3]+
				 pyin[8]*ptable[4]+
				 pyin[10]*ptable[5]+
				 pyin[12]*ptable[6]+
				 pyin[14]*ptable[7],
				 pyin[1]*ptable[0]+
				 pyin[3]*ptable[1]+
				 pyin[5]*ptable[2]+
				 pyin[7]*ptable[3]+
				 pyin[9]*ptable[4]+
				 pyin[11]*ptable[5]+
				 pyin[13]*ptable[6]+
				 pyin[15]*ptable[7]);
		/* else handle end effects with care */
		} else {

			sumr = sumi = 0.0;
			for (itable=0; itable<8; itable++,kyin++) {
				if (kyin<0) {
					yinir = yinlr;
					yinii = yinli;
				} else if (kyin>=nxin) {
					yinir = yinrr;
					yinii = yinri;
				} else {
				    yinir = crealf(yin[kyin]);
				    yinii = cimagf(yin[kyin]);
				}
				sumr += yinir*(*ptable);
				sumi += yinii*(*ptable++);
			}
			yout[ixout] = sf_cmplx(sumr,sumi);
		}
	}
}

void intt8r (int ntable       /* number of tabulated interpolation operators; ntable>=2 */, 
	     float table[][8] /* array of tabulated 8-point interpolation operators */,
	     int nxin         /* number of x values at which y(x) is input */, 
	     float dxin       /* x sampling interval for input y(x) */, 
	     float fxin       /* x value of first sample input */, 
	     float *yin       /* array of input y(x) values:  yin[0] = y(fxin), etc. */,
	     float yinl       /* value used to extrapolate yin values to left of yin[0] */, 
	     float yinr       /* value used to extrapolate yin values to right of yin[nxin-1] */, 
	     int nxout        /* number of x values a which y(x) is output */, 
	     float *xout      /*number of x values a which y(x) is output */, 
	     float *yout      /* array of output y(x) values:  yout[0] = y(xout[0]), etc. */)
/*< interpolation of a uniformly-sampled complex function y(x)
  via a table of 8-coefficient interpolators 

  The table of interpolation operators must be as follows:

Let d be the distance, expressed as a fraction of dxin, from a particular
xout value to the sampled location xin just to the left of xout.  Then,
for d = 0.0,

table[0][0:7] = 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0

are the weights applied to the 8 input samples nearest xout.
Likewise, for d = 1.0,

table[ntable-1][0:7] = 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0

are the weights applied to the 8 input samples nearest xout.  In general,
for d = (float)itable/(float)(ntable-1), table[itable][0:7] are the
weights applied to the 8 input samples nearest xout.  If the actual sample
distance d does not exactly equal one of the values for which interpolators
are tabulated, then the interpolator corresponding to the nearest value of
d is used.

Because extrapolation of the input function y(x) is defined by the left
and right values yinl and yinr, the xout values are not restricted to lie
within the range of sample locations defined by nxin, dxin, and fxin. >*/
{
	int ioutb,nxinm8,ixout,ixoutn,kyin,ktable,itable;
	float xoutb,xoutf,xouts,xoutn,frac,fntablem1,yini,sum,
		*yin0,*table00,*pyin,*ptable;

	/* compute constants */
	ioutb = -3-8;
	xoutf = fxin;
	xouts = 1.0/dxin;
	xoutb = 8.0-xoutf*xouts;
	fntablem1 = (float)(ntable-1);
	nxinm8 = nxin-8;
	yin0 = &yin[0];
	table00 = &table[0][0];

	/* loop over output samples */
	for (ixout=0; ixout<nxout; ixout++) {

		/* determine pointers into table and yin */
		xoutn = xoutb+xout[ixout]*xouts;
		ixoutn = (int)xoutn;
		kyin = ioutb+ixoutn;
		pyin = yin0+kyin;
		frac = xoutn-(float)ixoutn;
		ktable = frac>=0.0?frac*fntablem1+0.5:(frac+1.0)*fntablem1-0.5;
		ptable = table00+ktable*8;
		
		/* if totally within input array, use fast method */
		if (kyin>=0 && kyin<=nxinm8) {
			yout[ixout] = 
				pyin[0]*ptable[0]+
				pyin[1]*ptable[1]+
				pyin[2]*ptable[2]+
				pyin[3]*ptable[3]+
				pyin[4]*ptable[4]+
				pyin[5]*ptable[5]+
				pyin[6]*ptable[6]+
				pyin[7]*ptable[7];
		
		/* else handle end effects with care */
		} else {
	
			/* sum over 8 tabulated coefficients */
			for (itable=0,sum=0.0; itable<8; itable++,kyin++) {
				if (kyin<0)
					yini = yinl;
				else if (kyin>=nxin)
					yini = yinr;
				else
					yini = yin[kyin];
				sum += yini*(*ptable++);
			}
			yout[ixout] = sum;
		}
	}
}
