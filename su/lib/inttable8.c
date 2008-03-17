/* Interpolation tables for 8-point sinc */
/*
  Copyright © 2007, Colorado School of Mines,
  All rights reserved.
  
  
  Redistribution and use in source and binary forms, with or 
  without modification, are permitted provided that the following 
  conditions are met:
  
  *  Redistributions of source code must retain the above copyright 
  notice, this list of conditions and the following disclaimer.
  *  Redistributions in binary form must reproduce the above 
  copyright notice, this list of conditions and the following 
  disclaimer in the documentation and/or other materials provided 
  with the distribution.
  *  Neither the name of the Colorado School of Mines nor the names of
  its contributors may be used to endorse or promote products 
  derived from this software without specific prior written permission.
  
  Warranty Disclaimer:
  THIS SOFTWARE IS PROVIDED BY THE COLORADO SCHOOL OF MINES AND CONTRIBUTORS 
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
  COLORADO SCHOOL OF MINES OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
  POSSIBILITY OF SUCH DAMAGE.
  
  
  Export Restriction Disclaimer:
  We believe that CWP/SU: Seismic Un*x is a low technology product that does
  not appear on the Department of Commerce CCL list of restricted exports.
  Accordingly, we believe that our product meets the qualifications of
  an ECCN (export control classification number) of EAR99 and we believe
  it fits the qualifications of NRR (no restrictions required), and
  is thus not subject to export restrictions of any variety.
  
  Approved Reference Format:
  In publications, please refer to SU as per the following example:
  Cohen, J. K. and Stockwell, Jr. J. W., (200_), CWP/SU: Seismic Un*x 
  Release No. __: an open source software  package for seismic 
  research and processing, 
  Center for Wave Phenomena, Colorado School of Mines.
  
  Articles about SU in peer-reviewed journals:
  Saeki, T., (1999), A guide to Seismic Un*x (SU)(2)---examples of data processing (part 1), data input and preparation of headers, Butsuri-Tansa (Geophysical Exploration), vol. 52, no. 5, 465-477.
  Stockwell, Jr. J. W. (1999), The CWP/SU: Seismic Un*x Package, Computers and Geosciences, May 1999.
  Stockwell, Jr. J. W. (1997), Free Software in Education: A case study of CWP/SU: Seismic Un*x, The Leading Edge, July 1997.
  Templeton, M. E., Gough, C.A., (1998), Web Seismic Un*x: Making seismic reflection processing more accessible, Computers and Geosciences.
  
  Acknowledgements:
  SU stands for CWP/SU:Seismic Un*x, a processing line developed at Colorado 
  School of Mines, partially based on Stanford Exploration Project (SEP) 
  software.
*/

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
	     int nxout        /* number of x values at which y(x) is output */, 
	     float *xout      /* x values a which y(x) is output */, 
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
	     int nxout        /* number of x values at which y(x) is output */, 
	     float *xout      /* x values at which y(x) is output */, 
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
