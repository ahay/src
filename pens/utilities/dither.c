/*
  Copyright (C) 1987 The Board of Trustees of Stanford University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
 *
 *  source file:   ./filters/utilities/dither.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Steve Cole (SEP), September 1 1987
 *      Added method 4 (digital halftoning).
 * Joe Dellinger March 28 1988
 *	Move invras out of here and into greycorr, where it belonged
 *	in the first place.
 * Steve Cole (SEP), June 10 1988
 *      Replaced if blocks for each dithering method with a single switch.
 * W. Bauske IBM 03-26-91
 *	Apply SysV fixes for RS/6000
 * Joe Dellinger (BP Amoco), October 6 1999
 *	Floyd-Steinberg dithering routine ran off end of errline array,
 *	causing a core dump under Solaris.
 */

/*
 * this subroutine converts a single raster line to single bit
 * using one of the following algorithms:
 *
 * 1 random threshold
 * 2 256 element ordered dither (oriented at 0 degrees)
 * 3 Floyd-Steinberg minimized average error method
 * 4 32 element halftone (oriented at 45 degrees)
 *
 * Steve Cole, April 1987
 *
 */
#include <stdio.h>
#include <stdlib.h>

#include "../include/err.h"
#include "../include/params.h"
#include "../include/extern.h"

#include "error.h"

static int      pix_on = 0, pix_off = 7;
static float   *errline = NULL;
static int      ialloc = 0;

static float    alpha = 0.4375;
static float    beta = 0.1875;
static float    gama = 0.3125;
static float    delta = 0.0625;

static int      dith256[256] = {
    1, 128, 32, 160, 8, 136, 40, 168, 2, 130, 34, 162, 10, 138, 42, 170,
    192, 64, 224, 96, 200, 72, 232, 104, 194, 66, 226, 98, 202, 74, 234, 106,
    48, 176, 16, 144, 56, 184, 24, 152, 50, 178, 18, 146, 58, 186, 26, 154,
    240, 112, 208, 80, 248, 120, 216, 88, 242, 114, 210, 82, 250, 122, 218, 90,
    12, 140, 44, 172, 4, 132, 36, 164, 14, 142, 46, 174, 6, 134, 38, 166,
    204, 76, 236, 108, 196, 68, 228, 100, 206, 78, 238, 110, 198, 70, 230, 102,
    60, 188, 28, 156, 52, 180, 20, 148, 62, 190, 30, 158, 54, 182, 22, 150,
    252, 124, 220, 92, 244, 116, 212, 84, 254, 126, 222, 94, 246, 118, 214, 86,
    3, 131, 35, 163, 11, 139, 43, 171, 1, 129, 33, 161, 9, 137, 41, 169,
    195, 67, 227, 99, 203, 75, 235, 107, 193, 65, 225, 97, 201, 73, 233, 105,
    51, 179, 19, 147, 59, 187, 27, 155, 49, 177, 17, 145, 57, 185, 25, 153,
    243, 115, 211, 83, 251, 123, 219, 91, 241, 113, 209, 81, 249, 121, 217, 89,
    15, 143, 47, 175, 7, 135, 39, 167, 13, 141, 45, 173, 5, 133, 37, 165,
    207, 79, 239, 111, 199, 71, 231, 103, 205, 77, 237, 109, 197, 69, 229, 101,
    63, 191, 31, 159, 55, 183, 23, 151, 61, 189, 29, 157, 53, 181, 21, 149,
    254, 127, 223, 95, 247, 119, 215, 87, 253, 125, 221, 93, 245, 117, 213, 85
};
static int      halftone32[64] = {
    92, 100, 124, 148, 164, 156, 132, 108,
    28, 20, 76, 220, 228, 236, 180, 36,
    4, 12, 84, 212, 252, 244, 172, 44,
    52, 60, 116, 188, 204, 196, 140, 68,
    164, 156, 132, 108, 92, 100, 124, 148,
    228, 236, 180, 36, 28, 20, 76, 220,
    252, 244, 172, 44, 4, 12, 84, 212,
    204, 196, 140, 68, 52, 60, 116, 188
};


void
dithline (unsigned char *inpline, 
	  unsigned char *outline, int npixels, int linenum, int imethod)
/*< dithered line >*/
{
    int             greydata;
    int             i1, ipoint, jpoint;
    float           pixel, pixerr, nexterr;
    int             irand;

    switch (imethod)
    {
/* Random Dither */
	case 1:
	    for (i1 = 0; i1 < npixels; i1++)
	    {
		greydata = inpline[i1];
		irand = (rand () & 255);
		if (greydata > irand)
		{
		    outline[i1] = pix_off;
		}
		else
		{
		    outline[i1] = pix_on;
		}
	    }
	    break;

/* Ordered Dither */
	case 2:
	    for (i1 = 0; i1 < npixels; i1++)
	    {
		greydata = inpline[i1];
		ipoint = i1 % 16;
		jpoint = linenum % 16;
		ipoint = ipoint * 16 + jpoint;
		if (greydata > dith256[ipoint])
		{
		    outline[i1] = pix_off;
		}
		else
		{
		    outline[i1] = pix_on;
		}
	    }
	    break;

/* Floyd-Steinberg */
	case 3:
	    if (ialloc < npixels)
	    {
		if (ialloc > 0)
		{
		    free ((void *) errline);
		    ialloc = 0;
		}
		if ((errline = (float *) malloc ((unsigned) npixels * sizeof (float))) == NULL)
		{
		    ERR (FATAL, name, "Can't allocate space for Floyd-Steinberg\n");
		    return;
		}
		ialloc = npixels;
		for (i1 = 0; i1 < npixels; i1++)
		{
		    errline[i1] = 0.;
		}
	    }
	    nexterr = errline[0];
	    for (i1 = 0; i1 < npixels; i1++)
	    {
		pixel = inpline[i1];
		pixel += nexterr;
		if (pixel < 128)
		{
		    outline[i1] = pix_on;
		    pixerr = pixel;
		}
		else
		{
		    outline[i1] = pix_off;
		    pixerr = pixel - 255;
		}
		if (i1 < npixels - 1)
		{
		    nexterr = errline[i1 + 1] + pixerr * alpha;
		    errline[i1 + 1] = pixerr * delta;
		}
		if (i1 > 0)
		{
		    errline[i1 - 1] += pixerr * beta;
		}
		if (i1 == 0)
		    errline[i1] = pixerr * gama;
		else
		    errline[i1] += pixerr * gama;
	    }
	    break;

/* 32 element halftone at 45 degrees */
	case 4:
	default:
	    for (i1 = 0; i1 < npixels; i1++)
	    {
		greydata = inpline[i1];
		ipoint = i1 % 8;
		jpoint = linenum % 8;
		ipoint = ipoint * 8 + jpoint;
		if (greydata > halftone32[ipoint])
		{
		    outline[i1] = pix_off;
		}
		else
		{
		    outline[i1] = pix_on;
		}
	    }
	    break;
    }
    return;
}
