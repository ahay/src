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

void stoepd (int n /* dimension of system */, 
	     double *r /* [n] array of top row of Toeplitz matrix */, 
	     double *g /* [n] array of right-hand-side column vector */, 
	     double *f /* [n] array of solution (left-hand-side) column vector */, 
	     double *a /* [n] array[n] of solution to Ra=v (Claerbout, FGDP, p. 57) */)
/*< Solve a symmetric Toeplitz linear system of equations Rf=g for f
(double version) 

This routine does NOT solve the case when the main diagonal is zero, it
just silently returns.

The left column of the Toeplitz matrix is assumed to be equal to the top
row (as specified in r); i.e., the Toeplitz matrix is assumed symmetric. >*/
/******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	int i,j;
	double v,e,c,w,bot;

	if (r[0] == 0.0) return;

	a[0] = 1.0;
	v = r[0];
	f[0] = g[0]/r[0];

	for (j=1; j<n; j++) {
		
		/* solve Ra=v as in Claerbout, FGDP, p. 57 */
		a[j] = 0.0;
		f[j] = 0.0;
		for (i=0,e=0.0; i<j; i++)
			e += a[i]*r[j-i];
		c = e/v;
		v -= c*e;
		for (i=0; i<=j/2; i++) {
			bot = a[j-i]-c*a[i];
			a[i] -= c*a[j-i];
			a[j-i] = bot;
		}

		/* use a and v above to get f[i], i = 0,1,2,...,j */
		for (i=0,w=0.0; i<j; i++)
			w += f[i]*r[j-i];
		c = (w-g[j])/v;
		for (i=0; i<=j; i++)
			f[i] -= c*a[j-i];
	}
}

void stoepf (int n /* dimension of system */, 
	     float *r /* [n] array of top row of Toeplitz matrix */, 
	     float *g /* [n] array of right-hand-side column vector */, 
	     float *f /* [n] array of solution (left-hand-side) column vector */, 
	     float *a /* [n] array[n] of solution to Ra=v (Claerbout, FGDP, p. 57) */)
/*< 
Solve a symmetric Toeplitz linear system of equations Rf=g for f
(float version) 

This routine does NOT solve the case when the main diagonal is zero, it
just silently returns.

The left column of the Toeplitz matrix is assumed to be equal to the top
row (as specified in r); i.e., the Toeplitz matrix is assumed symmetric. >*/
/******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
{
	int i,j;
	float v,e,c,w,bot;

	if (r[0] == 0.0) return;

	a[0] = 1.0;
	v = r[0];
	f[0] = g[0]/r[0];

	for (j=1; j<n; j++) {
		
		/* solve Ra=v as in Claerbout, FGDP, p. 57 */
		a[j] = 0.0;
		f[j] = 0.0;
		for (i=0,e=0.0; i<j; i++)
			e += a[i]*r[j-i];
		c = e/v;
		v -= c*e;
		for (i=0; i<=j/2; i++) {
			bot = a[j-i]-c*a[i];
			a[i] -= c*a[j-i];
			a[j-i] = bot;
		}

		/* use a and v above to get f[i], i = 0,1,2,...,j */
		for (i=0,w=0.0; i<j; i++)
			w += f[i]*r[j-i];
		c = (w-g[j])/v;
		for (i=0; i<=j; i++)
			f[i] -= c*a[j-i];
	}
}
