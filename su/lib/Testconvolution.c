/* Copyright (c) Colorado School of Mines, 2011.
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
#include <stdio.h>

#include <rsf.h>

#include "convolution.h"

static void convolve_cwp1 (int lx, int ifx, float *x,
			   int ly, int ify, float *y, 
			   int lz, int ifz, float *z)
{
    int ilx=ifx+lx-1,ily=ify+ly-1,ilz=ifz+lz-1,i,j,jlow,jhigh;
    float sum;
	
    x -= ifx;  y -= ify;  z -= ifz;
    for (i=ifz; i<=ilz; ++i) {
	jlow = i-ily;  if (jlow<ifx) jlow = ifx;
	jhigh = i-ify;  if (jhigh>ilx) jhigh = ilx;
	for (j=jlow,sum=0.0; j<=jhigh; ++j)
	    sum += x[j]*y[i-j];
	z[i] = sum;
    }
}

#define IFC 0
#define ILC 60

static void pp1d (FILE *fp, char *title, int lx, int ifx, float x[])
/*****************************************************************************
Printer plot of a 1-dimensional array
******************************************************************************
Input:
fp		file pointer for output (e.g., stdout, stderr, etc.)
title		title of plot
lx		length of x
ifx		index of first x
x		array[lx] to be plotted
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/

{
    int ix,ibase,icx,ic;
    float xmin,xmax,xscale,xbase;

    /* print title */
    fprintf(fp,"\n");
    fprintf(fp,"%s\n",title);

    /* minimum and maximum x */
    xmin = x[0];
    xmax = x[0];
    for (ix=1; ix<lx; ix++) {
	xmin = SF_MIN(xmin,x[ix]);
	xmax = SF_MAX(xmax,x[ix]);
    }
    fprintf(fp,"minimum = %g\n",xmin);
    fprintf(fp,"maximum = %g\n",xmax);

    /* determine scale factors and shifts for converting x values to *s */
    if (xmin==xmax)
	xscale = 1.0;
    else
	xscale = (ILC-IFC)/(xmax-xmin);
    if (xmin<0.0 && xmax<0.0) {
	ibase = ILC;
	xbase = xmax;
    } else if (xmin<0.0 && xmax>=0.0) {
	ibase = IFC+(0.0-xmin)*xscale;
	xbase = 0.0;
    } else {
	ibase = IFC;
	xbase = xmin;
    }

    /* loop over x values */
    for (ix=0; ix<lx; ix++) {

	/* determine column corresponding to x value */
	icx = ibase+SF_NINT((x[ix]-xbase)*xscale);
	icx = SF_MAX(IFC,SF_MIN(ILC,icx));

	/* print the index, x value, and row of *s */
	fprintf(fp,"%4d %13.6e ",ifx+ix,x[ix]);
	for (ic=IFC; ic<SF_MIN(ibase,icx); ic++)
	    fprintf(fp," ");
	for (ic=SF_MIN(ibase,icx); ic<=SF_MAX(ibase,icx); ic++)
	    fprintf(fp,"*");
	fprintf(fp,"\n");
    }
}

int main(void)
{
    int it, lx,ly,lz,ifx,ify,ifz,i;
    float *x,*y,*z,*z1;
	
    /* loop over tests */
    for (it=0; it < 1000; it++) {
	sf_warning("test %d;", it);
	
	/* make parameters */
	lx = 1+100*genrand_real1();
	ly = 1+100*genrand_real1();
	lz = 1+100*genrand_real1();
	ifx = -100+200*genrand_real1();
	ify = -100+200*genrand_real1();
	ifz = -100+200*genrand_real1();
		
	/* allocate space */
	x = sf_floatalloc(lx);
	y = sf_floatalloc(ly);
	z = sf_floatalloc(lz);
	z1 = sf_floatalloc(lz);
			
	/* fill x and y arrays */
	for (i=0; i<lx; i++)
	    x[i] = i+1;
	for (i=0; i<ly; i++)
	    y[i] = i+1;

	/* convolve and check */
	convolve_cwp1(lx,ifx,x,ly,ify,y,lz,ifz,z1);
	convolve_cwp(lx,ifx,x,ly,ify,y,lz,ifz,z);
	for (i=0; i<lz; ++i)
	    if (z[i]!=z1[i]) break;
	if (i<lz) {
	    printf("%10d %10d %10d %10d %10d %10d\n",
		   lx,ifx,ly,ify,lz,ifz);
	    printf("z1[%d]=%g != z[%d]=%g\n",
		   i+ifz,z1[i],i+ifz,z[i]);
	    pp1d(stdout,"simple",lz,ifz,z1);
	    pp1d(stdout,"optimized",lz,ifz,z);
	    exit(-1);
	}
				
	/* free space */
	free(x);
	free(y);
	free(z);
	free(z1);
    }
    sf_warning(".");
}
