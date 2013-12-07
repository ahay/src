/* Convert yx to xy */
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
YXTOXY - Compute a regularly-sampled, monotonically increasing function x(y)
	from a regularly-sampled, monotonically increasing function y(x) by
	inverse linear interpolation.

yxtoxy		compute a regularly sampled function x(y) from a regularly
		sampled, monotonically increasing function y(x)

******************************************************************************
Function Prototype:
void yxtoxy (int nx, float dx, float fx, float y[], 
	int ny, float dy, float fy, float xylo, float xyhi, float x[]);

******************************************************************************
Input:
nx		number of samples of y(x)
dx		x sampling interval; dx>0.0 is required
fx		first x
y		array[nx] of y(x) values; y[0] < y[1] < ... < y[nx-1] required
ny		number of samples of x(y)
dy		y sampling interval; dy>0.0 is required
fy		first y
xylo		x value assigned to x(y) when y is less than smallest y(x)
xyhi		x value assigned to x(y) when y is greater than largest y(x)

Output:
x		array[ny] of x(y) values

******************************************************************************
 Notes:
 User must ensure that:
 (1) dx>0.0 && dy>0.0
 (2) y[0] < y[1] < ... < y[nx-1] */

/*******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
/**************** end self doc ********************************/
#include "yxtoxy.h"

void yxtoxy (int nx /* number of samples of y(x) */, 
	     float dx /* x sampling interval; dx>0.0 is required */, 
	     float fx /* first x */, 
	     const float *y /* array[nx] of y(x) values; y[0] < y[1] < ... < y[nx-1] required */, 
	     int ny /* number of samples of x(y) */, 
	     float dy /* y sampling interval; dy>0.0 is required */, 
	     float fy /* first y */, 
	     float xylo /* x value assigned to x(y) when y is less than smallest y(x) */, 
	     float xyhi /* x value assigned to x(y) when y is greater than largest y(x) */, 
	     float *x /* array[ny] of x(y) values */)
/*< Compute a regularly-sampled, monotonically increasing function x(y)
from a regularly-sampled, monotonically increasing function y(x) by
inverse linear interpolation.
***
Notes:
User must ensure that:
(1) dx>0.0 && dy>0.0
(2) y[0] < y[1] < ... < y[nx-1] >*/
/******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	int nxi,nyo,jxi1,jxi2,jyo;
	float dxi,fxi,dyo,fyo,fyi,yo,xi1,yi1,yi2; 

	nxi = nx; dxi = dx; fxi = fx;
	nyo = ny; dyo = dy; fyo = fy;
	fyi = y[0];

	/* loop over output y less than smallest input y */
	for (jyo=0,yo=fyo; jyo<nyo; jyo++,yo+=dyo) {
		if (yo>=fyi) break;
		x[jyo] = xylo;
	}

	/* loop over output y between smallest and largest input y */
	if (jyo==nyo-1 && yo==fyi) {
		x[jyo++] = fxi;
		yo += dyo;
	}
	jxi1 = 0;
	jxi2 = 1;
	xi1 = fxi;
	while (jxi2<nxi && jyo<nyo) {
		yi1 = y[jxi1];
		yi2 = y[jxi2];
		if (yi1<=yo && yo<=yi2) {
			x[jyo++] = xi1+dxi*(yo-yi1)/(yi2-yi1);
			yo += dyo;
		} else {
			jxi1++;
			jxi2++;
			xi1 += dxi;
		}
	}

	/* loop over output y greater than largest input y */
	while (jyo<nyo)
		x[jyo++] = xyhi;
}
