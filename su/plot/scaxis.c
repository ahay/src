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
#include <math.h>

#include <rsf.h>

void scaxis (float x1     /* first x value */, 
	     float x2     /* second x value */, 
	     int *nxnum   /* number of numbered values */, 
	     float *dxnum /* increment between numbered values (dxnum>0.0) */, 
	     float *fxnum /* first numbered value */)
/*< compute a readable scale for use in plotting axes

scaxis attempts to honor the user-specified nxnum.  However, nxnum
will be modified if necessary for readability.  Also, fxnum and nxnum
will be adjusted to compensate for roundoff error; in particular, 
fxnum will not be less than xmin-eps, and fxnum+(nxnum-1)*dxnum 
will not be greater than xmax+eps, where eps = 0.0001*(xmax-xmin).
xmin is the minimum of x1 and x2.  xmax is the maximum of x1 and x2. >*/
/******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 01/13/89
*****************************************************************************/
{
	int n,i,iloga;
	float d,f,rdint[4],eps,a,b,xmin,xmax;

	/* set readable intervals */
	rdint[0] = 1.0;  rdint[1] = 2.0;  rdint[2] = 5.0;  rdint[3] = 10.0;

	/* handle x1==x2 as a special case */
	if  (x1==x2) {
		*nxnum = 1;
		*dxnum = 1.0;
		*fxnum = x1;
		return;
	}

	/* determine minimum and maximum x */
	xmin = (x1<x2)?x1:x2;
	xmax = (x1>x2)?x1:x2;
	
	/* get desired number of numbered values */
	n = *nxnum;
	n = (2>n)?2:n;
	
	/* determine output parameters, adjusted for roundoff */
	a = (xmax-xmin)/(float)(n-1);
	iloga = (int)log10f(a);
	if (a<1.0) iloga = iloga - 1;
	b = a/powf(10.0,(double)iloga);
	for (i=0; i<3 && b>=sqrtf(rdint[i]*rdint[i+1]); i++);
	d = rdint[i]*powf(10.0,(float)iloga);
	f = ((int)(xmin/d))*d-d;
	eps = 0.0001*(xmax-xmin);
	while(f<(xmin-eps))
		 f = f+d;
	n = 1+(int)((xmax+eps-f)/d); 
        
	/* set output parameters before returning */
	*nxnum = n;
	*dxnum = d;
	*fxnum = f;
}
