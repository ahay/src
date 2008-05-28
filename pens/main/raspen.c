/* vplot filter for ppm format output. */
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
 *  source file:   ./filters/raslib/rasattr.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Stewart A. Levin (MRDC/DRL), September 8, 1988
 *      Return DOVPLOT_CONT on exit.
 */
/*
 *
 *  source file:   ./filters/raslib/rasclose.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger (SEP), May 7, 1989
 *	RGB option added.
 * Chuck Karish, August 5, 1989
 * 	For non-SEP users, don't try to open colfile unless it has been
 *	specified!
 */

/*
 *
 *  source file:   ./filters/raslib/raserase.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger (SEP), May 7, 1989
 *	RGB option added. Stew Levin's Color tektronix code modified and
 *	incorporated.
 * Joe Dellinger (SEP), August 5, 1989
 *	Fixed bug that caused only first page of GRF format output to
 *	actually print.
 * Dave Nichols (SEP) May 2 1992
 *      Allow VPLOTSPOOLDIR environment variable to override PEN_SPOOL
 * Dave Nichols (SEP) July 10 1992
 *      Added ppm output.
 */

/*
 *
 *  source file:   ./filters/raslib/rasopen.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger (SEP), Sept 25 1988
 *	"ppi" undocumented, and read as an integer and not a float!
 * Joe Dellinger (SEP), May 7 1989
 *	Cleaned up to make raslib work with Tektronix color plotter.
 *	Incorporated portions of Stew Levin's color plotter code.
 * Joe Dellinger (SEP), August 5, 1989
 *	Use the variable "default_out" to remember the ORIGINAL
 *	value of isatty(fileno(pltout)).
 * Chuck Karish, 5 August 1989
 *	Non-SEP, colfile defaults to empty string instead of "colfile".
 * Joe Dellinger, 17 Jan 1990
 *	ARG!!! I'm changing the default behavior back so it matches
 *	the documentation. I'm also fixing the bug introduced by the
 *	change so that "colfile" was not even settable from the command
 *	line.
 * Dave Nichols, 10 July 1992
 *	Added support for ppm output, (ppmpen, Ppmpen).
 */

/*
 *
 *  source file:   ./filters/raslib/rasreset.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */


/*
 * control graphics attributes
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <ppm.h>

#include "../include/attrcom.h"
#include "../include/closestat.h"
#include "../include/err.h"
#include "../include/extern.h"
#include "../include/params.h"
#include "../include/erasecom.h"

#include "../genlib/genpen.h"
#include "../utilities/util.h"

#include "dovplot.h"
#include "init_vplot.h"

#include "_ras.h"
#include "raspen.h"

char            name[] = "raspen";
#include "rasdoc.h"

struct device   dev =
{

    /* control routines */
    rasopen,		/* open */
    rasreset,		/* reset */
    genmessage,		/* message */
    raserase,		/* erase */
    rasclose,		/* close */
    
    /* high level output */
    rasvector,		/* vector */
    genmarker,		/* marker */
    gentext,		/* text */
    genpatarea,		/* area */
    genraster,		/* raster */
    genpoint,		/* point */
    rasattr,		/* attributes */
    
    /* input */
    gen_dovplot,            /* reader */
    nullgetpoint,		/* getpoint */
    nullinteract,		/* interact */
    
    /* low level output */
    nullplot,		/* plot */
    nullclose,		/* startpoly */
    nullmidpoly,		/* midpoly */
    nullclose		/* endpoly */
};

extern FILE    *pltout;

int             color_table[NCOLOR][3], rascolor;

unsigned char  *image;
extern float    aspect_ratio;
extern float    pixels_per_inch;
int             rasor = 0;
char            colfile[60];

float           o_pixels_per_inch, o_aspect_ratio;
static bool     default_out = true;

extern int      num_col;

void rasattr (int command, int value, int v1, int v2, int v3)
/*< attr >*/
{
    switch (command)
    {
	case SET_COLOR:
	    rascolor = value;
	    break;
	case SET_COLOR_TABLE:
	    color_table[value][0] = v1;
	    color_table[value][1] = v2;
	    color_table[value][2] = v3;
	    break;
	default:
	    break;
    }
}

extern int      color_table[NCOLOR][3];
extern char     colfile[];

void rasclose (int status)
/*< close >*/
{
    switch (status)
    {
	case CLOSE_NORMAL:
	    break;
	case CLOSE_FLUSH:
	    fflush (pltout);
	    break;
	default:
	    break;
    }
}

void raserase (int command)
/*< erase >*/
{
    
    switch (command)
    {
	case ERASE_START:
	default:
	    break;
	case ERASE_MIDDLE:
	case ERASE_END:
/*
 * Output raster file, and then CLEAR it out (erase)
 */
	    ras_write ();
	    zap ();
	    break;
	case ERASE_BREAK:
/*
 * Output raster file, but don't CLEAR it!
 */
	    ras_write ();
	    break;
    }
}

void ras_write (void)
/*< write >*/
{
    static bool called=false;
    pixel * pixrow;
    unsigned char r,g,b;
    unsigned char *ptr;
    register pixel *pP;
    int row, col;

    if(called) return;
    called=true;
    
    ppm_writeppminit( pltout, dev_xmax, dev_ymax, (pixval)255, 0);
    pixrow = ppm_allocrow( dev_xmax );
    for ( row = 0, ptr=image;  row < dev_ymax; ++row )
    {
        for ( col = 0, pP = pixrow; col < dev_xmax; ++col, ++pP ){
	    r = *(ptr++);
	    g = *(ptr++);
	    b = *(ptr++);
	    PPM_ASSIGN( *pP, r, g, b );
	}
        ppm_writeppmrow( pltout, pixrow, dev_xmax, (pixval)255, 0 );
    }
    pm_close( pltout );
}


void rasopen (int argc, char* argv[])
/*< open >*/
{
    char newpath[60];

    txfont = DEFAULT_HARDCOPY_FONT;
    txprec = DEFAULT_HARDCOPY_PREC;

    ppm_init(&argc,argv);


/*
 * physical device parameters for ppm format output
 */
    o_pixels_per_inch = pixels_per_inch = 100.;
    dev_xmax = 10 * pixels_per_inch;
    dev_ymax = 7.5 * pixels_per_inch;
    dev_xmin = 0;
    dev_ymin = 0;
    o_aspect_ratio = aspect_ratio = 1.;

/*
 * device capabilities
 */
    need_end_erase = true;
    smart_clip = false;
    num_col = NCOLOR;

    sf_getfloat ("aspect", &aspect_ratio);
    /* aspect ratio */
    sf_getfloat ("ppi", &pixels_per_inch);
    /* pixels per inch */
    dev_xmax *= pixels_per_inch / o_pixels_per_inch;
    dev_ymax *= (o_aspect_ratio / aspect_ratio) *
     (pixels_per_inch / o_pixels_per_inch);
    sf_getint ("n1", &dev_xmax);
    sf_getint ("n2", &dev_ymax);
    /* image size */

    /*
     * Allocate space for image 
     */
    image = sf_ucharalloc (dev_xmax * dev_ymax * 3);

    default_out = (bool) isatty(fileno(pltout));

    if (default_out)
    {
	sprintf (newpath, "%s", "raster_file");
	pltout = fopen (newpath, "w");
	if (pltout == NULL)
	    ERR (FATAL, name, "can't open file %s\n", newpath);
    }
}

void rasreset (void)
/*< reset >*/
{
    int             value;
    zap ();
    for (value = 0; value < NCOLOR; value++)
    {
	color_table[value][0] = -1;
    }
    for (value = 0; value < 8; value++)
    {
	color_table[value][0] = MAX_GUN * ((value & 2) / 2);
	color_table[value][1] = MAX_GUN * ((value & 4) / 4);
	color_table[value][2] = MAX_GUN * ((value & 1) / 1);
    }
    color_table[0][0] =256;
    color_table[0][1] =0;
    color_table[0][2] =0;
}

void zap (void)
/*< Zero image >*/
{
    unsigned char  *p;
    for (p = image; p < &image[dev_xmax * dev_ymax * 3]; p++)
	*p = 0;
}

#define  RGB

#define WRITEIT(A,B) \
Image3 (A, B, 0) = color_table[rascolor][0]; \
Image3 (A, B, 1) = color_table[rascolor][1]; \
Image3 (A, B, 2) = color_table[rascolor][2]

void rasvector (int x1, int y1, int x2, int y2, int nfat, int dashon)
/*< vector >*/
{
    int             test, tmp, x, y;
    double          slope, fx, fx3, fy, fy3;

/*
 * Vector rasterizes the line defined by the endpoints (x1,y1) and (x2,y2).
 * If 'nfat' is nonzero then draw parallel lines to fatten the line, by
 * recursive calls to vector.
 */

    if (nfat < 0)
	return;

    if (dashon)
    {
	dashvec (x1, y1, x2, y2, nfat, dashon);
	return;
    }

    if (nfat)
    {
	if (clip (&x1, &y1, &x2, &y2))
	    return;

	fatvec (x1, y1, x2, y2, nfat, dashon);
	return;
    }

    if (clip (&x1, &y1, &x2, &y2))
	return;

/* Beware checks out of bounds, since the coordinate system may have rotated */

    test = (abs (x2 - x1) >= abs (y2 - y1));

    if (test)
    {
	if (x1 == x2)
	{
	    /* Just a point */
	    WRITEIT (x1, y1);
	    return;
	}
	else
	if (x1 > x2)
	{
	    tmp = x1;
	    x1 = x2;
	    x2 = tmp;
	    tmp = y1;
	    y1 = y2;
	    y2 = tmp;
	}
	slope = (double) (y2 - y1) / (double) (x2 - x1);
	fy3 = y1;

	for (x = x1, fy = fy3; x < x2; x++, fy += slope)
	{
	    y = fy + .5;	/* OK rounding, since always positive */
	    WRITEIT (x, y);
	}
	WRITEIT (x2, y2);
	return;
    }
    else
    {
	/* y1 can't equal y2 here */
	if (y1 > y2)
	{
	    tmp = x1;
	    x1 = x2;
	    x2 = tmp;
	    tmp = y1;
	    y1 = y2;
	    y2 = tmp;
	}
	slope = (double) (x2 - x1) / (double) (y2 - y1);
	fx3 = x1;

	for (y = y1, fx = fx3; y < y2; y++, fx += slope)
	{
	    x = fx + .5;
	    WRITEIT (x, y);
	}
	WRITEIT (x2, y2);
	return;
    }
}
