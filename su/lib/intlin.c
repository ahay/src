/* Linear interpolation */
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

/*********************** self documentation **********************/
/*****************************************************************************
INTLIN - evaluate y(x) via linear interpolation of y(x[0]), y(x[1]), ...

intlin		evaluate y(x) via linear interpolation of y(x[0]), y(x[1]), ...

******************************************************************************
Function Prototype:
void intlin (int nin, float xin[], float yin[], float yinl, float yinr,
	int nout, float xout[], float yout[]);

******************************************************************************
Input:
nin		length of xin and yin arrays
xin		array[nin] of monotonically increasing or decreasing x values
yin		array[nin] of input y(x) values
yinl		value used to extraplate y(x) to left of input yin values
yinr		value used to extraplate y(x) to right of input yin values
nout		length of xout and yout arrays
xout		array[nout] of x values at which to evaluate y(x)

Output:
yout		array[nout] of linearly interpolated y(x) values

******************************************************************************
Notes:
xin values must be monotonically increasing or decreasing.

Extrapolation of the function y(x) for xout values outside the range
spanned by the xin values in performed as follows:

	For monotonically increasing xin values,
		yout=yinl if xout<xin[0], and yout=yinr if xout>xin[nin-1].

	For monotonically decreasing xin values, 
		yout=yinl if xout>xin[0], and yout=yinr if xout<xin[nin-1].

If nin==1, then the monotonically increasing case is used.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
/**************** end self doc ********************************/

#include "intlin.h"
#include "xindex.h"

void intlin (int nin /* length of xin and yin arrays */, 
	     const float *xin /* array[nin] of monotonically increasing or decreasing x values */, 
	     const float *yin /* array[nin] of input y(x) values */, 
	     float yinl /* value used to extraplate y(x) to left of input yin values */, 
	     float yinr /* value used to extraplate y(x) to right of input yin values */, 
	     int nout /* length of xout and yout arrays */, 
	     const float *xout /* array[nout] of x values at which to evaluate y(x) */, 
	     float *yout /* array[nout] of linearly interpolated y(x) values */)
/*< evaluate y(x) via linear interpolation of y(x[0]), y(x[1]), ...

Notes:
xin values must be monotonically increasing or decreasing.

Extrapolation of the function y(x) for xout values outside the range
spanned by the xin values in performed as follows:

	For monotonically increasing xin values,
		yout=yinl if xout<xin[0], and yout=yinr if xout>xin[nin-1].

	For monotonically decreasing xin values, 
		yout=yinl if xout>xin[0], and yout=yinr if xout<xin[nin-1].

        If nin==1, then the monotonically increasing case is used. >*/
/******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	static int idx;
	int jout;
	float x;

	/* if input x values are monotonically increasing, then */
	if (xin[0]<=xin[nin-1]) {
		for (jout=0; jout<nout; jout++) {
			x = xout[jout];
			if (x<xin[0])
				yout[jout] = yinl;
			else if (x>xin[nin-1])
				yout[jout] = yinr;
			else if (x==xin[nin-1] || nin==1)
				yout[jout] = yin[nin-1];
			else {
				xindex(nin,xin,x,&idx);
				yout[jout] = yin[idx]+(x-xin[idx])
					*(yin[idx+1]-yin[idx])
					/(xin[idx+1]-xin[idx]);
			}
		}
	
	/* else, if input x values are monotonically decreasing, then */
	} else {
		for (jout=0; jout<nout; jout++) {
			x = xout[jout];
			if (x>xin[0])
				yout[jout] = yinl;
			else if (x<xin[nin-1])
				yout[jout] = yinr;
			else if (x==xin[nin-1] || nin==1)
				yout[jout] = yin[nin-1];
			else {
				xindex(nin,xin,x,&idx);
				yout[jout] = yin[idx]+(x-xin[idx])
					*(yin[idx+1]-yin[idx])
					/(xin[idx+1]-xin[idx]);
			}
		}
	}
}
