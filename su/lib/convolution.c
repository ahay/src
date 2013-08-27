/* Convolution routines. */

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

/* prototype of function used internally */
static void convolve_cwp_s(int,int,float*,int,int,float*,int,int,float*);

void convolve_cwp (int lx, int ifx, float *x,
		   int ly, int ify, float *y,
		   int lz, int ifz, float *z)
/*<
****************************************************************************
Compute z = x convolved with y; i.e.,

ifx+lx-1
z[i] =   sum    x[j]*y[i-j]  ;  i = ifz,...,ifz+lz-1
j=ifx
******************************************************************************
Input:
lx		length of x array
ifx		sample index of first x
x		array[lx] to be convolved with y
ly		length of y array
ify		sample index of first y
y		array[ly] with which x is to be convolved
lz		length of z array
ifz		sample index of first z

Output:
z		array[lz] containing x convolved with y
******************************************************************************
Notes:
The x samples are contained in x[0], x[1], ..., x[lx-1]; likewise for
the y and z samples.  The sample indices of the first x, y, and z values
determine the location of the origin for each array.  For example, if
z is to be a weighted average of the nearest 5 samples of y, one might
use 
...
x[0] = x[1] = x[2] = x[3] = x[4] = 1.0/5.0;
conv(5,-2,x,lx,0,y,ly,0,z);
...
In this example, the filter x is symmetric, with index of first sample = -2.

This function is optimized for architectures that can simultaneously perform
a multiply, add, and one load from memory; e.g., the IBM RISC System/6000.
Because, for each value of i, it accumulates the convolution sum z[i] in a
scalar, this function is not likely to be optimal for vector architectures.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 11/23/91
****************************************************************************>*/
#ifdef SIMPLE_CONV
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
#else
{
    int ilx=ifx+lx-1,ily=ify+ly-1,ilz=ifz+lz-1,
	i,j,ilow,ihigh,jlow,jhigh;
    float sa,sb,xa,xb,ya,yb,*t;

    /* if x is longer than y, swap x and y */
    if (lx>ly) {
	i = ifx;  ifx = ify;  ify = i;
	i = ilx;  ilx = ily;  ily = i;
	i = lx;  lx = ly;  ly = i;
	t = x;  x = y;  y = t;
    }
	
    /* handle short x with special code */
    if (lx>=1 && lx<=30) {
	convolve_cwp_s(lx,ifx,x,ly,ify,y,lz,ifz,z);
	return;
    }
	
    /* adjust pointers for indices of first samples */
    x -= ifx;
    y -= ify;
    z -= ifz;
		
    /* OFF LEFT:  i < ify+ifx */
	
    /* zero output for all i */
    ilow = ifz;
    ihigh = ify+ifx-1;  if (ihigh>ilz) ihigh = ilz;
    for (i=ilow; i<=ihigh; ++i)
	z[i] = 0.0;

    /* ROLLING ON:  ify+ifx <= i < ify+ilx */
	
    /* if necessary, do one i so that number of j in overlap is odd */
    if (i<ify+ilx && i<=ilz) {
	jlow = ifx;
	jhigh = i-ify;
	if ((jhigh-jlow)%2) {
	    sa = 0.0;
	    for (j=jlow; j<=jhigh; ++j)
		sa += x[j]*y[i-j];
	    z[i++] = sa;
	}
    }
	
    /* loop over pairs of i and j */
    ilow = i;
    ihigh = ilx+ify-1;  if (ihigh>ilz) ihigh = ilz;
    jlow = ifx;
    jhigh = ilow-ify;
    for (i=ilow; i<ihigh; i+=2,jhigh+=2) {
	sa = sb = 0.0;
	xb = x[jhigh+1];
	yb = 0.0;
	for (j=jhigh; j>=jlow; j-=2) {
	    sa += xb*yb;
	    ya = y[i-j];
	    sb += xb*ya;
	    xa = x[j];
	    sa += xa*ya;
	    yb = y[i+1-j];
	    sb += xa*yb;
	    xb = x[j-1];
	}
	z[i] = sa;
	z[i+1] = sb;
    }
	
    /* if number of i is odd */
    if (i==ihigh) {
	jlow = ifx;
	jhigh = i-ify;
	sa = 0.0;
	for (j=jlow; j<=jhigh; ++j)
	    sa += x[j]*y[i-j];
	z[i++] = sa;
    }
	
    /* MIDDLE:  ify+ilx <= i <= ily+ifx */
	
    /* determine limits for i and j */
    ilow = i;
    ihigh = ily+ifx;  if (ihigh>ilz) ihigh = ilz;
    jlow = ifx;
    jhigh = ilx;
	
    /* if number of j is even, do j in pairs with no leftover */
    if ((jhigh-jlow)%2) {
	for (i=ilow; i<ihigh; i+=2) {
	    sa = sb = 0.0;
	    yb = y[i+1-jlow];
	    xa = x[jlow];
	    for (j=jlow; j<jhigh; j+=2) {
		sb += xa*yb;
		ya = y[i-j];
		sa += xa*ya;
		xb = x[j+1];
		sb += xb*ya;
		yb = y[i-1-j];
		sa += xb*yb;
		xa = x[j+2];
	    }
	    z[i] = sa;
	    z[i+1] = sb;
	}
	
	/* else, number of j is odd, so do j in pairs with leftover */
    } else {
	for (i=ilow; i<ihigh; i+=2) {
	    sa = sb = 0.0;
	    yb = y[i+1-jlow];
	    xa = x[jlow];
	    for (j=jlow; j<jhigh; j+=2) {
		sb += xa*yb;
		ya = y[i-j];
		sa += xa*ya;
		xb = x[j+1];
		sb += xb*ya;
		yb = y[i-1-j];
		sa += xb*yb;
		xa = x[j+2];
	    }
	    z[i] = sa+x[jhigh]*y[i-jhigh];
	    z[i+1] = sb+x[jhigh]*y[i+1-jhigh];
	}
    }
	
    /* if number of i is odd */
    if (i==ihigh) {
	sa = 0.0;
	for (j=jlow; j<=jhigh; ++j)
	    sa += x[j]*y[i-j];
	z[i++] = sa;
    }

    /* ROLLING OFF:  ily+ifx < i <= ily+ilx */
	
    /* if necessary, do one i so that number of j in overlap is even */
    if (i<=ily+ilx && i<=ilz) {
	jlow = i-ily;
	jhigh = ilx;
	if (!((jhigh-jlow)%2)) {
	    sa = 0.0;
	    for (j=jlow; j<=jhigh; ++j)
		sa += x[j]*y[i-j];
	    z[i++] = sa;
	}
    }
	
    /* number of j is now even, so loop over both i and j in pairs */
    ilow = i;
    ihigh = ily+ilx;  if (ihigh>ilz) ihigh = ilz;
    jlow = ilow-ily;
    jhigh = ilx-2; /* Dave's new patch */
    for (i=ilow; i<ihigh; i+=2,jlow+=2) {
	sa = sb = 0.0;
	xa = x[jlow];
	yb = 0.0;
	for (j=jlow; j<jhigh; j+=2) {
	    sb += xa*yb;
	    ya = y[i-j];
	    sa += xa*ya;
	    xb = x[j+1];
	    sb += xb*ya;
	    yb = y[i-1-j];
	    sa += xb*yb;
	    xa = x[j+2];
	}
	sb += xa*yb;
	ya = y[i-j];
	sa += xa*ya;
	xb = x[j+1];
	sb += xb*ya;
	yb = y[i-1-j];
	sa += xb*yb;
	z[i] = sa;
	z[i+1] = sb;
    }
	
    /* if number of i is odd */
    if (i==ihigh) {
	jlow = i-ily;
	jhigh = ilx;
	sa = 0.0;
	for (j=jlow; j<=jhigh; ++j)
	    sa += x[j]*y[i-j];
	z[i++] = sa;
    }
	
    /* OFF RIGHT:  ily+ilx < i */
	
    /* zero output for all i */
    ilow = i;
    ihigh = ilz;
    for (i=ilow; i<=ihigh; ++i)
	z[i] = 0.0;
}

/* internal function optimized for short x */
static void convolve_cwp_s (int lx, int ifx, float *x,
			    int ly, int ify, float *y, 
			    int lz, int ifz, float *z)
{
    int ilx=ifx+lx-1,ily=ify+ly-1,ilz=ifz+lz-1,
	i,j,ilow,ihigh,jlow,jhigh;
    float x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,
	x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,
	x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,
	ya,yb,z0,z1,sum;
	
    x -= ifx;
    y -= ify;
    z -= ifz;
	
		
    /* OFF LEFT:  i < ifx+ify */
    ilow = ifz;
    ihigh = ify+ifx-1;  if (ihigh>ilz) ihigh = ilz;
    for (i=ilow; i<=ihigh; ++i)
	z[i] = 0.0;
	
    /* ROLLING ON:  ify+ifx <= i < ify+ilx */
    ilow = ify+ifx;  if (ilow<ifz) ilow = ifz;
    ihigh = ify+ilx-1;  if (ihigh>ilz) ihigh = ilz;
    jlow = ifx;
    jhigh = ilow-ify;
    for (i=ilow; i<=ihigh; ++i,++jhigh) {
	for (j=jlow,sum=0.0; j<=jhigh; ++j)
	    sum += x[j]*y[i-j];
	z[i] = sum;
    }
	
    /* MIDDLE:  ify+ilx <= i <= ily+ifx */
    ilow = ify+ilx;  if (ilow<ifz) ilow = ifz;
    ihigh = ily+ifx;  if (ihigh>ilz) ihigh = ilz;
    if (lx==1) {
	x0 = x[ifx];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==2) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==3) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==4) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==5) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==6) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==7) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==8) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==9) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==10) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==11) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==12) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	x11 = x[ifx+11];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
	    ya = y[i-ifx-11];  z0 += x11*ya;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==13) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	x11 = x[ifx+11];
	x12 = x[ifx+12];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
	    ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
	    yb = y[i-ifx-12];  z0 += x12*yb;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==14) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	x11 = x[ifx+11];
	x12 = x[ifx+12];
	x13 = x[ifx+13];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
	    ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
	    yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
	    ya = y[i-ifx-13];  z0 += x13*ya;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==15) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	x11 = x[ifx+11];
	x12 = x[ifx+12];
	x13 = x[ifx+13];
	x14 = x[ifx+14];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
	    ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
	    yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
	    ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
	    yb = y[i-ifx-14];  z0 += x14*yb;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==16) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	x11 = x[ifx+11];
	x12 = x[ifx+12];
	x13 = x[ifx+13];
	x14 = x[ifx+14];
	x15 = x[ifx+15];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
	    ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
	    yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
	    ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
	    yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
	    ya = y[i-ifx-15];  z0 += x15*ya;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==17) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	x11 = x[ifx+11];
	x12 = x[ifx+12];
	x13 = x[ifx+13];
	x14 = x[ifx+14];
	x15 = x[ifx+15];
	x16 = x[ifx+16];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
	    ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
	    yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
	    ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
	    yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
	    ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
	    yb = y[i-ifx-16];  z0 += x16*yb;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==18) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	x11 = x[ifx+11];
	x12 = x[ifx+12];
	x13 = x[ifx+13];
	x14 = x[ifx+14];
	x15 = x[ifx+15];
	x16 = x[ifx+16];
	x17 = x[ifx+17];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
	    ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
	    yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
	    ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
	    yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
	    ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
	    yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
	    ya = y[i-ifx-17];  z0 += x17*ya;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==19) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	x11 = x[ifx+11];
	x12 = x[ifx+12];
	x13 = x[ifx+13];
	x14 = x[ifx+14];
	x15 = x[ifx+15];
	x16 = x[ifx+16];
	x17 = x[ifx+17];
	x18 = x[ifx+18];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
	    ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
	    yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
	    ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
	    yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
	    ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
	    yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
	    ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
	    yb = y[i-ifx-18];  z0 += x18*yb;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==20) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	x11 = x[ifx+11];
	x12 = x[ifx+12];
	x13 = x[ifx+13];
	x14 = x[ifx+14];
	x15 = x[ifx+15];
	x16 = x[ifx+16];
	x17 = x[ifx+17];
	x18 = x[ifx+18];
	x19 = x[ifx+19];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
	    ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
	    yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
	    ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
	    yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
	    ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
	    yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
	    ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
	    yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
	    ya = y[i-ifx-19];  z0 += x19*ya;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==21) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	x11 = x[ifx+11];
	x12 = x[ifx+12];
	x13 = x[ifx+13];
	x14 = x[ifx+14];
	x15 = x[ifx+15];
	x16 = x[ifx+16];
	x17 = x[ifx+17];
	x18 = x[ifx+18];
	x19 = x[ifx+19];
	x20 = x[ifx+20];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
	    ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
	    yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
	    ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
	    yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
	    ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
	    yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
	    ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
	    yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
	    ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
	    yb = y[i-ifx-20];  z0 += x20*yb;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==22) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	x11 = x[ifx+11];
	x12 = x[ifx+12];
	x13 = x[ifx+13];
	x14 = x[ifx+14];
	x15 = x[ifx+15];
	x16 = x[ifx+16];
	x17 = x[ifx+17];
	x18 = x[ifx+18];
	x19 = x[ifx+19];
	x20 = x[ifx+20];
	x21 = x[ifx+21];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
	    ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
	    yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
	    ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
	    yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
	    ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
	    yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
	    ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
	    yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
	    ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
	    yb = y[i-ifx-20];  z0 += x20*yb;  z1 += x21*yb;
	    ya = y[i-ifx-21];  z0 += x21*ya;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==23) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	x11 = x[ifx+11];
	x12 = x[ifx+12];
	x13 = x[ifx+13];
	x14 = x[ifx+14];
	x15 = x[ifx+15];
	x16 = x[ifx+16];
	x17 = x[ifx+17];
	x18 = x[ifx+18];
	x19 = x[ifx+19];
	x20 = x[ifx+20];
	x21 = x[ifx+21];
	x22 = x[ifx+22];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
	    ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
	    yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
	    ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
	    yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
	    ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
	    yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
	    ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
	    yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
	    ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
	    yb = y[i-ifx-20];  z0 += x20*yb;  z1 += x21*yb;
	    ya = y[i-ifx-21];  z0 += x21*ya;  z1 += x22*ya;
	    yb = y[i-ifx-22];  z0 += x22*yb;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==24) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	x11 = x[ifx+11];
	x12 = x[ifx+12];
	x13 = x[ifx+13];
	x14 = x[ifx+14];
	x15 = x[ifx+15];
	x16 = x[ifx+16];
	x17 = x[ifx+17];
	x18 = x[ifx+18];
	x19 = x[ifx+19];
	x20 = x[ifx+20];
	x21 = x[ifx+21];
	x22 = x[ifx+22];
	x23 = x[ifx+23];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
	    ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
	    yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
	    ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
	    yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
	    ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
	    yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
	    ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
	    yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
	    ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
	    yb = y[i-ifx-20];  z0 += x20*yb;  z1 += x21*yb;
	    ya = y[i-ifx-21];  z0 += x21*ya;  z1 += x22*ya;
	    yb = y[i-ifx-22];  z0 += x22*yb;  z1 += x23*yb;
	    ya = y[i-ifx-23];  z0 += x23*ya;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==25) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	x11 = x[ifx+11];
	x12 = x[ifx+12];
	x13 = x[ifx+13];
	x14 = x[ifx+14];
	x15 = x[ifx+15];
	x16 = x[ifx+16];
	x17 = x[ifx+17];
	x18 = x[ifx+18];
	x19 = x[ifx+19];
	x20 = x[ifx+20];
	x21 = x[ifx+21];
	x22 = x[ifx+22];
	x23 = x[ifx+23];
	x24 = x[ifx+24];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
	    ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
	    yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
	    ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
	    yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
	    ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
	    yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
	    ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
	    yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
	    ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
	    yb = y[i-ifx-20];  z0 += x20*yb;  z1 += x21*yb;
	    ya = y[i-ifx-21];  z0 += x21*ya;  z1 += x22*ya;
	    yb = y[i-ifx-22];  z0 += x22*yb;  z1 += x23*yb;
	    ya = y[i-ifx-23];  z0 += x23*ya;  z1 += x24*ya;
	    yb = y[i-ifx-24];  z0 += x24*yb;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==26) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	x11 = x[ifx+11];
	x12 = x[ifx+12];
	x13 = x[ifx+13];
	x14 = x[ifx+14];
	x15 = x[ifx+15];
	x16 = x[ifx+16];
	x17 = x[ifx+17];
	x18 = x[ifx+18];
	x19 = x[ifx+19];
	x20 = x[ifx+20];
	x21 = x[ifx+21];
	x22 = x[ifx+22];
	x23 = x[ifx+23];
	x24 = x[ifx+24];
	x25 = x[ifx+25];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
	    ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
	    yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
	    ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
	    yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
	    ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
	    yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
	    ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
	    yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
	    ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
	    yb = y[i-ifx-20];  z0 += x20*yb;  z1 += x21*yb;
	    ya = y[i-ifx-21];  z0 += x21*ya;  z1 += x22*ya;
	    yb = y[i-ifx-22];  z0 += x22*yb;  z1 += x23*yb;
	    ya = y[i-ifx-23];  z0 += x23*ya;  z1 += x24*ya;
	    yb = y[i-ifx-24];  z0 += x24*yb;  z1 += x25*yb;
	    ya = y[i-ifx-25];  z0 += x25*ya;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==27) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	x11 = x[ifx+11];
	x12 = x[ifx+12];
	x13 = x[ifx+13];
	x14 = x[ifx+14];
	x15 = x[ifx+15];
	x16 = x[ifx+16];
	x17 = x[ifx+17];
	x18 = x[ifx+18];
	x19 = x[ifx+19];
	x20 = x[ifx+20];
	x21 = x[ifx+21];
	x22 = x[ifx+22];
	x23 = x[ifx+23];
	x24 = x[ifx+24];
	x25 = x[ifx+25];
	x26 = x[ifx+26];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
	    ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
	    yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
	    ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
	    yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
	    ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
	    yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
	    ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
	    yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
	    ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
	    yb = y[i-ifx-20];  z0 += x20*yb;  z1 += x21*yb;
	    ya = y[i-ifx-21];  z0 += x21*ya;  z1 += x22*ya;
	    yb = y[i-ifx-22];  z0 += x22*yb;  z1 += x23*yb;
	    ya = y[i-ifx-23];  z0 += x23*ya;  z1 += x24*ya;
	    yb = y[i-ifx-24];  z0 += x24*yb;  z1 += x25*yb;
	    ya = y[i-ifx-25];  z0 += x25*ya;  z1 += x26*ya;
	    yb = y[i-ifx-26];  z0 += x26*yb;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==28) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	x11 = x[ifx+11];
	x12 = x[ifx+12];
	x13 = x[ifx+13];
	x14 = x[ifx+14];
	x15 = x[ifx+15];
	x16 = x[ifx+16];
	x17 = x[ifx+17];
	x18 = x[ifx+18];
	x19 = x[ifx+19];
	x20 = x[ifx+20];
	x21 = x[ifx+21];
	x22 = x[ifx+22];
	x23 = x[ifx+23];
	x24 = x[ifx+24];
	x25 = x[ifx+25];
	x26 = x[ifx+26];
	x27 = x[ifx+27];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
	    ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
	    yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
	    ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
	    yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
	    ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
	    yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
	    ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
	    yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
	    ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
	    yb = y[i-ifx-20];  z0 += x20*yb;  z1 += x21*yb;
	    ya = y[i-ifx-21];  z0 += x21*ya;  z1 += x22*ya;
	    yb = y[i-ifx-22];  z0 += x22*yb;  z1 += x23*yb;
	    ya = y[i-ifx-23];  z0 += x23*ya;  z1 += x24*ya;
	    yb = y[i-ifx-24];  z0 += x24*yb;  z1 += x25*yb;
	    ya = y[i-ifx-25];  z0 += x25*ya;  z1 += x26*ya;
	    yb = y[i-ifx-26];  z0 += x26*yb;  z1 += x27*yb;
	    ya = y[i-ifx-27];  z0 += x27*ya;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==29) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	x11 = x[ifx+11];
	x12 = x[ifx+12];
	x13 = x[ifx+13];
	x14 = x[ifx+14];
	x15 = x[ifx+15];
	x16 = x[ifx+16];
	x17 = x[ifx+17];
	x18 = x[ifx+18];
	x19 = x[ifx+19];
	x20 = x[ifx+20];
	x21 = x[ifx+21];
	x22 = x[ifx+22];
	x23 = x[ifx+23];
	x24 = x[ifx+24];
	x25 = x[ifx+25];
	x26 = x[ifx+26];
	x27 = x[ifx+27];
	x28 = x[ifx+28];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
	    ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
	    yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
	    ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
	    yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
	    ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
	    yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
	    ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
	    yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
	    ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
	    yb = y[i-ifx-20];  z0 += x20*yb;  z1 += x21*yb;
	    ya = y[i-ifx-21];  z0 += x21*ya;  z1 += x22*ya;
	    yb = y[i-ifx-22];  z0 += x22*yb;  z1 += x23*yb;
	    ya = y[i-ifx-23];  z0 += x23*ya;  z1 += x24*ya;
	    yb = y[i-ifx-24];  z0 += x24*yb;  z1 += x25*yb;
	    ya = y[i-ifx-25];  z0 += x25*ya;  z1 += x26*ya;
	    yb = y[i-ifx-26];  z0 += x26*yb;  z1 += x27*yb;
	    ya = y[i-ifx-27];  z0 += x27*ya;  z1 += x28*ya;
	    yb = y[i-ifx-28];  z0 += x28*yb;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    } else if (lx==30) {
	x0 = x[ifx];
	x1 = x[ifx+1];
	x2 = x[ifx+2];
	x3 = x[ifx+3];
	x4 = x[ifx+4];
	x5 = x[ifx+5];
	x6 = x[ifx+6];
	x7 = x[ifx+7];
	x8 = x[ifx+8];
	x9 = x[ifx+9];
	x10 = x[ifx+10];
	x11 = x[ifx+11];
	x12 = x[ifx+12];
	x13 = x[ifx+13];
	x14 = x[ifx+14];
	x15 = x[ifx+15];
	x16 = x[ifx+16];
	x17 = x[ifx+17];
	x18 = x[ifx+18];
	x19 = x[ifx+19];
	x20 = x[ifx+20];
	x21 = x[ifx+21];
	x22 = x[ifx+22];
	x23 = x[ifx+23];
	x24 = x[ifx+24];
	x25 = x[ifx+25];
	x26 = x[ifx+26];
	x27 = x[ifx+27];
	x28 = x[ifx+28];
	x29 = x[ifx+29];
	for (i=ilow; i<=ihigh-1; i+=2) {
	    ya = y[i+1-ifx];  z1 = x0*ya;
	    yb = y[i-ifx];  z0 = x0*yb;  z1 += x1*yb;
	    ya = y[i-ifx-1];  z0 += x1*ya;  z1 += x2*ya;
	    yb = y[i-ifx-2];  z0 += x2*yb;  z1 += x3*yb;
	    ya = y[i-ifx-3];  z0 += x3*ya;  z1 += x4*ya;
	    yb = y[i-ifx-4];  z0 += x4*yb;  z1 += x5*yb;
	    ya = y[i-ifx-5];  z0 += x5*ya;  z1 += x6*ya;
	    yb = y[i-ifx-6];  z0 += x6*yb;  z1 += x7*yb;
	    ya = y[i-ifx-7];  z0 += x7*ya;  z1 += x8*ya;
	    yb = y[i-ifx-8];  z0 += x8*yb;  z1 += x9*yb;
	    ya = y[i-ifx-9];  z0 += x9*ya;  z1 += x10*ya;
	    yb = y[i-ifx-10];  z0 += x10*yb;  z1 += x11*yb;
	    ya = y[i-ifx-11];  z0 += x11*ya;  z1 += x12*ya;
	    yb = y[i-ifx-12];  z0 += x12*yb;  z1 += x13*yb;
	    ya = y[i-ifx-13];  z0 += x13*ya;  z1 += x14*ya;
	    yb = y[i-ifx-14];  z0 += x14*yb;  z1 += x15*yb;
	    ya = y[i-ifx-15];  z0 += x15*ya;  z1 += x16*ya;
	    yb = y[i-ifx-16];  z0 += x16*yb;  z1 += x17*yb;
	    ya = y[i-ifx-17];  z0 += x17*ya;  z1 += x18*ya;
	    yb = y[i-ifx-18];  z0 += x18*yb;  z1 += x19*yb;
	    ya = y[i-ifx-19];  z0 += x19*ya;  z1 += x20*ya;
	    yb = y[i-ifx-20];  z0 += x20*yb;  z1 += x21*yb;
	    ya = y[i-ifx-21];  z0 += x21*ya;  z1 += x22*ya;
	    yb = y[i-ifx-22];  z0 += x22*yb;  z1 += x23*yb;
	    ya = y[i-ifx-23];  z0 += x23*ya;  z1 += x24*ya;
	    yb = y[i-ifx-24];  z0 += x24*yb;  z1 += x25*yb;
	    ya = y[i-ifx-25];  z0 += x25*ya;  z1 += x26*ya;
	    yb = y[i-ifx-26];  z0 += x26*yb;  z1 += x27*yb;
	    ya = y[i-ifx-27];  z0 += x27*ya;  z1 += x28*ya;
	    yb = y[i-ifx-28];  z0 += x28*yb;  z1 += x29*yb;
	    ya = y[i-ifx-29];  z0 += x29*ya;
	    z[i+1] = z1;
	    z[i] = z0;
	}
    }
    if (ihigh>=ilow && (ihigh-ilow)%2==0) {
	ilow = ihigh;
	jlow = ifx;
	jhigh = ilx;
	for (i=ilow; i<=ihigh; ++i) {
	    for (j=jlow,sum=0.0; j<=jhigh; ++j)
		sum += x[j]*y[i-j];
	    z[i] = sum;
	}
    }
	
    /* ROLLING OFF:  ily+ifx < i <= ily+ilx */
    ilow = ily+ifx+1;  if (ilow<ifz) ilow = ifz;
    ihigh = ily+ilx;  if (ihigh>ilz) ihigh = ilz;
    jlow = ilow-ily;
    jhigh = ilx;
    for (i=ilow; i<=ihigh; ++i,++jlow) {
	for (j=jlow,sum=0.0; j<=jhigh; ++j)
	    sum += x[j]*y[i-j];
	z[i] = sum;
    }
	
    /* OFF RIGHT:  ily+ilx < i */
    ilow = ily+ilx+1;  if (ilow<ifz) ilow = ifz;
    ihigh = ilz;
    for (i=ilow; i<=ihigh; ++i)
	z[i] = 0.0;
}
#endif


