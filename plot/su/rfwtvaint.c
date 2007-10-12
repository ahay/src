/* Copyright (c) Colorado School of Mines, 2007.*/
/* All rights reserved.                       */

#include <rsf.h>

#include "intsinc8.h"

void rfwtvaint (int n               /* number of samples in array to rasterize */, 
		float *z            /* [n] array to rasterize */, 
		float zmin          /* z values below zmin will be clipped */, 
		float zmax          /* z values above zmax will be clipped */, 
		float zbase         /* z values between zbase and zmax will be filled */,
		int yzmin           /* horizontal raster coordinate corresponding to zmin */, 
		int yzmax           /* horizontal raster coordinate corresponding to zmax */, 
		int xfirst          /* vertical raster coordinate of z[0] */, 
		int xlast           /* vertical raster coordinate of z[n-1] */,
		int wiggle          /* =0 for no wiggle (VA only); =1 for wiggle (with VA)
				       wiggle 2<=wiggle<=5 for solid/grey coloring of VA option
				       shade of grey: wiggle=2 light grey, wiggle=5 black */, 
		int nbpr            /* number of bytes per row of bits */, 
		unsigned char *bits /* pointer to first (top,left) byte in image */)
/*< Rasterize a float array as wiggle-trace-variable-area.

The raster coordinate of the (top,left) bit in the image is (0,0).
In other words, x increases downward and y increases to the right.
Raster scan lines run from left to right, and from top to bottom.
Therefore, xfirst, xlast, yzmin, and yzmax should not be less than 0.
Likewise, yzmin and yzmax should not be greater than nbpr*8-1, and 
care should be taken to ensure that xfirst and xlast do not cause bits 
to be set outside (off the bottom) of the image. 

Variable area fill is performed on the right-hand (increasing y) side
of the wiggle.  If yzmin is greater than yzmax, then z values between
zmin will be plotted to the right of zmax, and z values between zbase
and zmin are filled.  Swapping yzmin and yzmax is an easy way to 
reverse the polarity of a wiggle.

The variable "endian" must have a value of 1 or 0. If this is
not a case an error is returned. >*/

/******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 07/01/89
Modified:  Craig Artley, Colorado School of Mines, 04/14/92
           Fixed bug in computing yoffset.  Previously, when zmin==zmax
           the rasterized trace was shifted to the left by one trace.
MODIFIED: Paul Michaels, Boise State University, 29 December 2000
          added solid/grey color scheme for peaks/troughs  wiggle=2 option
*****************************************************************************/
{
	int i,y,
		ymin,ymax,ybase;
	int igrey,ideci;
	float yscale,yoffset,zthis;
	register int bit;
	register unsigned char *byte;

	static float *xout  , *yout;
	static int    nx = 0;
	float        *yin   ,  fxin , yinl, yinr, dxin, x0, deltax;
	int           nxin  ,  nxout;

	/* if solid/grey coloring desired      */
	if (wiggle>=2)
	{  igrey=SF_ABS(wiggle); wiggle=1; }
	else
	{  igrey=0; }

	/* Compute the number of raster scan lines. */
	nxout = SF_ABS(xlast - xfirst + 1);

	/* If the # of scan lines has increased from a previous call,
	   allocate more*/
	if ( nxout > nx ) {

	  /* If a previous call allocated space for output values,
	     free them. */
	  if ( nx > 0 ) {
	    free (xout);
	    free (yout);
	  }

	  /* Allocate space for the scan line x values and interpolated
	     z values. */

		xout = (float *)calloc ((size_t)nxout, sizeof(float));
		yout = (float *)calloc ((size_t)nxout, sizeof(float));
	  nx   = nxout;
	}

	nxin   = n       ; /* There are n z-values.  */
	dxin   = 1.0     ; /* We go from index 0 to index n - 1 in steps
			      of 1.0  */
	fxin   = 0.0     ; /* The first index is 0.  */
	yin    = z       ; /* The input array is the z array.  */
	yinl   = z[0]    ; /* Set the values to the left of the array
			      to z[0]. */
	yinr   = z[n - 1]; /* Set the values to the right of the array
			      to z[n-1].*/

	deltax = (float)(n - 1) / (float)(nxout - 1);
	x0     =  0.0;
	if ( xfirst > xlast ) { /* If the z array is to be output backwards, */
		x0 = (float)(n - 1);  /* Then the first output index is n-1, */
		deltax = -deltax;     /*   and we decrement rather than
					   increment. */
	}
	for (i = 0; i < nxout; i++) /* Load the indices of the output values.*/
		xout[i] = x0 + (float)i * deltax;

	ints8r (nxin, dxin, fxin, yin, yinl, yinr, nxout, xout, yout);
		

	/* determine min and max y coordinates */
	ymin = (yzmin<yzmax)?yzmin:yzmax;
	ymax = (yzmax>yzmin)?yzmax:yzmin;

	/* restrict min and max y coordinates */
	ymin = (ymin>0)?ymin:0;
	ymax = (ymax<nbpr*8-1)?ymax:nbpr*8-1;
	
	/* determine y scale factor and offset */
	yscale = (zmax!=zmin)?(yzmax-yzmin)/(zmax-zmin):1.0;
	yoffset = (zmax!=zmin)?yzmin-zmin*yscale:0.5*(yzmin+yzmax);
	
	/* determine base y coordinate */
	ybase = yoffset+zbase*yscale;
	ybase = (ybase>ymin)?ybase:ymin;
	ybase = (ybase<ymax)?ybase:ymax;
	
	/* loop over scan lines */
	for (i = 0; i < nxout; i++) {
		zthis = yout[i];
		y     = yoffset + zthis * yscale;

		/* apply clip */
		if (y < ymin) y = ymin;
		if (y > ymax) y = ymax;
			
		/* determine the bit and byte */
		/* original: bit = 7-y&7; */
		bit = (7-y)&7;

		/* Tony Kocurko: Had been "bits+x*nbpr+(y>>3)".*/
		byte = bits+i*nbpr+(y>>3);

		/* if wiggle or filling, then set the bit */
		if (wiggle || y>ybase) {
/*			if (endian==0)
				*byte |= 1<<(-bit+7);
				else if (endian==1) */
				*byte |= 1<<bit;
/*			else
  fprintf(stderr,"endian must equal either 0 or 1\n"); */
		}

		
		/* while y greater than base, set more bits (SOLID FILL PEAKS) */
		while (y>ybase) {
			y-=1;
			bit+=1;
			if (bit>=8) {
				byte--;
				bit = 0;
			}
/*			if (endian==0)
				*byte |= 1<<(-bit+7);
				else if (endian==1) */
				*byte |= 1<<bit;
/*			else
  fprintf(stderr,"endian must equal either 0 or 1\n"); */
		}  /*  endwhile  */

		/* while y less than base, set more bits (GREY FILL TROUGHS) */
	        if (igrey>0)
	        {
                ideci=6-igrey;
                if (ideci<1) ideci=1;
		
			while (y<ybase) {
				y+=ideci;
				bit-=ideci;
				if (bit<0) {
					byte++;
					bit = 7;
				}
/*				if (endian==0)
					*byte |= 1<<(-bit+7);
					else if (endian==1) */
					*byte |= 1<<bit;
/*				else
  fprintf(stderr,"endian must equal either 0 or 1\n"); */
			}  /* endwhile */
		}  /* endif igrey  */

	} /* next scan line  */
}  /*   end rfwtvaint   */
