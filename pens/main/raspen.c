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
#include "../include/enum.h"

#include "../genlib/genpen.h"
#include "../utilities/util.h"

#include "dovplot.h"
#include "init_vplot.h"

#define	Image(IX,IY)		image[(IX)+xmax*(ymax-1-(IY))]
#define	Image3(IX,IY,RGB)	image[(RGB)+(IX)*3+xmax*3*(ymax-1-(IY))]
#define	Min(IX,IY)		((IX) < (IY) ? (IX) : (IY))
#define	Max(IX,IY)		((IX) > (IY) ? (IX) : (IY))
#define NCOLOR 256		/* number of colors */

#include "raspen.h"

char            name[] = "raspen";
#include "rasdoc.h"

#include "_device.h"

extern FILE    *pltout;

int             color_table[NCOLOR][3], rascolor;

extern char     colfile[];

unsigned char  *image;
extern float    aspect_ratio;
extern float    pixels_per_inch;
int             rasor = 0;
char            colfile[60];

static bool     default_out = true;
static int xmax, ymax;

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
    
    ppm_writeppminit( pltout, xmax, ymax, (pixval)255, 0);
    pixrow = ppm_allocrow( xmax );
    for ( row = 0, ptr=image;  row < ymax; ++row )
    {
        for ( col = 0, pP = pixrow; col < xmax; ++col, ++pP ){
	    r = *(ptr++);
	    g = *(ptr++);
	    b = *(ptr++);
	    PPM_ASSIGN( *pP, r, g, b );
	}
        ppm_writeppmrow( pltout, pixrow, xmax, (pixval)255, 0 );
    }
    pm_close( pltout );
}

void opendev (int argc, char* argv[])
/*< open >*/
{
    char newpath[60];
    float pixels_per_inch, aspect_ratio;

    ppm_init(&argc,argv);

    dev.txfont = DEFAULT_HARDCOPY_FONT;
    dev.txprec = DEFAULT_HARDCOPY_PREC;
    dev.brake = BREAK_IGNORE;

    dev.reset = rasreset;
    dev.erase = raserase;
    dev.close = rasclose;

    dev.area = genpatarea;
    dev.plot = rasplot;
    dev.attributes = rasattr;

    /* physical device parameters for ppm format output */
    pixels_per_inch = dev.pixels_per_inch = 100.;
    dev.xmax = VP_STANDARD_HEIGHT * dev.pixels_per_inch;
    dev.ymax = VP_SCREEN_RATIO * dev.xmax;
    aspect_ratio = dev.aspect_ratio = 1.;

    /* device capabilities */
    dev.need_end_erase = true;
    dev.num_col = NCOLOR;

    sf_getfloat ("aspect", &dev.aspect_ratio);
    /* aspect ratio */
    sf_getfloat ("ppi", &dev.pixels_per_inch);
    /* pixels per inch */
    dev.xmax *= dev.pixels_per_inch / pixels_per_inch;
    dev.ymax *= (aspect_ratio / dev.aspect_ratio) *
     (dev.pixels_per_inch / pixels_per_inch);
    sf_getint ("n1", &dev.xmax);
    sf_getint ("n2", &dev.ymax);
    /* image size */

    xmax = dev.xmax+1;
    ymax = dev.ymax+1;

    /*
     * Allocate space for image 
     */
    image = sf_ucharalloc (xmax * ymax * 3);

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
    for (p = image; p < &image[xmax * ymax * 3]; p++)
	*p = 0;
}

#define  RGB

#define WRITEIT(A,B) \
Image3 (A, B, 0) = color_table[rascolor][0]; \
Image3 (A, B, 1) = color_table[rascolor][1]; \
Image3 (A, B, 2) = color_table[rascolor][2]


static void rasline (int x1, int y1, int x2, int y2)
/*< Rasterizes the line defined by the endpoints (x1,y1) and (x2,y2). >*/
{
    int             test, tmp, x, y;
    double          slope, fx, fx3, fy, fy3;

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
	else if (x1 > x2)
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

void rasplot(int x, int y, int draw)
/*< plot >*/
{
    static int oldx = 0, oldy = 0;

    if (draw)
    {
	rasline(oldx, oldy, x, y);
    } else {
	dev.lost = 0;
    }

    oldx = x;
    oldy = y;
}
