/* Find index of x */
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
XINDEX - determine index of x with respect to an array of x values

xindex		determine index of x with respect to an array of x values

******************************************************************************
Input:
nx		number of x values in array ax
ax		array[nx] of monotonically increasing or decreasing x values
x		the value for which index is to be determined
index		index determined previously (used to begin search)

Output:
index		for monotonically increasing ax values, the largest index
		for which ax[index]<=x, except index=0 if ax[0]>x;
		for monotonically decreasing ax values, the largest index
		for which ax[index]>=x, except index=0 if ax[0]<x

******************************************************************************
Notes:
This function is designed to be particularly efficient when called
repeatedly for slightly changing x values; in such cases, the index 
returned from one call should be used in the next.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 12/25/89
*****************************************************************************/
/**************** end self doc ********************************/

void xindex (int nx /* number of x values in array ax */,
	     const float *ax /* array[nx] of monotonically increasing or decreasing x values */, 
	     float x /* the value for which index is to be determined */, 
	     int *index /* index determined previously (used to begin search) */)
/*< Determine index of x with respect to an array of x values

Output:
index		for monotonically increasing ax values, the largest index
		for which ax[index]<=x, except index=0 if ax[0]>x;
		for monotonically decreasing ax values, the largest index
		for which ax[index]>=x, except index=0 if ax[0]<x
******************************************************************************
Notes:
This function is designed to be particularly efficient when called
repeatedly for slightly changing x values; in such cases, the index 
returned from one call should be used in the next. >*/
/******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 12/25/89
*****************************************************************************/
{
	int lower,upper,middle,step;

	/* initialize lower and upper indices and step */
	lower = *index;
	if (lower<0) lower = 0;
	if (lower>=nx) lower = nx-1;
	upper = lower+1;
	step = 1;

	/* if x values increasing */
	if (ax[nx-1]>ax[0]) {

		/* find indices such that ax[lower] <= x < ax[upper] */
		while (lower>0 && ax[lower]>x) {
			upper = lower;
			lower -= step;
			step += step;
		}
		if (lower<0) lower = 0;
		while (upper<nx && ax[upper]<=x) {
			lower = upper;
			upper += step;
			step += step;
		}
		if (upper>nx) upper = nx;

		/* find index via bisection */
		while ((middle=(lower+upper)>>1)!=lower) {
			if (x>=ax[middle])
				lower = middle;
			else
				upper = middle;
		}

	/* else, if not increasing */
	} else {

		/* find indices such that ax[lower] >= x > ax[upper] */
		while (lower>0 && ax[lower]<x) {
			upper = lower;
			lower -= step;
			step += step;
		}
		if (lower<0) lower = 0;
		while (upper<nx && ax[upper]>=x) {
			lower = upper;
			upper += step;
			step += step;
		}
		if (upper>nx) upper = nx;

		/* find index via bisection */
		while ((middle=(lower+upper)>>1)!=lower) {
			if (x<=ax[middle])
				lower = middle;
			else
				upper = middle;
		}
	}

	/* return lower index */
	*index = lower;
}
