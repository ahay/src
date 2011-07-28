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

#ifdef _PPM
#include <ppm.h>
#endif

#ifdef _TIFF
#include <tiffio.h>
#endif

#ifdef _JPEG
#include <jpeglib.h>
#endif

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

static bool     default_out = true, light = false;
static int xmax, ymax;

#ifdef _TIFF
static TIFF *tiffout;
static char *tiffname;
#endif

#ifdef _JPEG
static struct jpeg_compress_struct *jpeg;
static struct jpeg_error_mgr jpeg_err;
#endif

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

#ifdef _PPM

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

#endif

#ifdef _TIFF

void ras_write (void)
/*< write >*/
{
    static bool called=false;
    int y, nbuf;
    unsigned char *p;
    FILE *tiffin;
    char buf[BUFSIZ];

    if(called) return;
    called=true;

    p = image;
    for (y = 0; y < ymax; y++) {
	if (TIFFWriteScanline(tiffout, p, y, 0) < 0) 
	    ERR (FATAL, name, "Trouble writing TIFF file");
	p += xmax * 3;
    }

    TIFFClose(tiffout);
    tiffin = fopen(tiffname,"rb");

    while (1) {
	nbuf = fread(buf,1,BUFSIZ,tiffin);
	if (nbuf <= 0) break;
	fwrite(buf,1,nbuf,pltout);
    }

    fclose(tiffin);
    fclose(pltout);
    
    unlink(tiffname);
}

#endif

#ifdef _JPEG

void ras_write (void)
/*< write >*/
{
    static bool called=false;
    JSAMPROW scan;

    if(called) return;
    called=true;

    jpeg_start_compress(jpeg, TRUE);
    while ((int) jpeg->next_scanline < ymax) {
	scan = image + jpeg->next_scanline * xmax * 3;
	jpeg_write_scanlines(jpeg, &scan, 1);
    }
    jpeg_finish_compress(jpeg);

    jpeg_destroy_compress(jpeg);
    fclose(pltout);
}

#endif

void opendev (int argc, char* argv[])
/*< open >*/
{
    char newpath[60];
    const char *color;
    float pixels_per_inch, aspect_ratio;

    dev.txfont = DEFAULT_HARDCOPY_FONT;
    dev.txprec = DEFAULT_HARDCOPY_PREC;
    dev.brake = BREAK_IGNORE;

    dev.reset = rasreset;
    dev.erase = raserase;
    dev.close = rasclose;

    dev.area = genpatarea;
    dev.plot = rasplot;
    dev.attributes = rasattr;

    pixels_per_inch = dev.pixels_per_inch = 100.;
    xmax = VP_STANDARD_HEIGHT * dev.pixels_per_inch;
    ymax = VP_SCREEN_RATIO * xmax;
    aspect_ratio = dev.aspect_ratio = 1.;

    /* device capabilities */
    dev.need_end_erase = true;
    dev.num_col = NCOLOR;

    sf_getfloat ("aspect", &dev.aspect_ratio);
    /* aspect ratio */
    sf_getfloat ("ppi", &dev.pixels_per_inch);
    /* pixels per inch */
    xmax *= dev.pixels_per_inch / pixels_per_inch;
    ymax *= (aspect_ratio / dev.aspect_ratio) *
	(dev.pixels_per_inch / pixels_per_inch);
    sf_getint ("n1", &xmax);
    sf_getint ("n2", &ymax);
    /* image size */

    dev.xmax = xmax-1;
    dev.ymax = ymax-1;

    if (NULL == (color = sf_getstring("bgcolor"))) color="black";
    /* background color */
    light = (bool) (color[0] == 'w' || color[0] == 'l');

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

#ifdef _PPM
    ppm_init(&argc,argv);
#endif

#ifdef _TIFF
    fclose(sf_tempfile(&tiffname,"w"));
    tiffout = TIFFOpen(tiffname,"wb");

    if (tiffout == NULL)
	ERR (FATAL, name, "can't open file %s\n", tiffname);

    TIFFSetField(tiffout,TIFFTAG_IMAGEWIDTH,xmax);
    TIFFSetField(tiffout,TIFFTAG_IMAGELENGTH,ymax);
    TIFFSetField(tiffout,TIFFTAG_SAMPLESPERPIXEL,3);
    TIFFSetField(tiffout,TIFFTAG_BITSPERSAMPLE,8);
    TIFFSetField(tiffout,TIFFTAG_ORIENTATION,ORIENTATION_TOPLEFT);
    TIFFSetField(tiffout,TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField(tiffout,TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tiffout, TIFFTAG_ROWSPERSTRIP, 1);
#endif

#ifdef _JPEG
    jpeg = (struct jpeg_compress_struct *) sf_alloc(1,sizeof(*jpeg));
    jpeg->err = jpeg_std_error(&jpeg_err);
    jpeg_create_compress(jpeg);
    jpeg_stdio_dest(jpeg, pltout);
    jpeg->image_width = xmax;
    jpeg->image_height = ymax;
    jpeg->input_components = 3;
    jpeg->in_color_space = JCS_RGB;
    jpeg_set_defaults(jpeg);
#endif
}

void rasreset (void)
/*< reset >*/
{
    int             value, r, g, b;
    zap ();
    if (light) {
	color_table[0][0] =255;
	color_table[0][1] =255;
	color_table[0][2] =255;
    } else {
	color_table[0][0] =0;
	color_table[0][1] =0;
	color_table[0][2] =0;
    }
    for (value = 1; value < 8; value++)
    {
	r = MAX_GUN * ((value & 2) / 2);
	g = MAX_GUN * ((value & 4) / 4);
	b = MAX_GUN * ((value & 1) / 1);

	if (light) {
	    color_table[value][0] = 255-r;
	    color_table[value][1] = 255-g;
	    color_table[value][2] = 255-b;
	} else {
	    color_table[value][0] = r;
	    color_table[value][1] = g;
	    color_table[value][2] = b;
	}
    }
    for (value = 8; value < NCOLOR; value++)
    {
	color_table[value][0] = -1;
    }   
}

void zap (void)
/*< Zero image >*/
{
    unsigned char  *p;
    for (p = image; p < &image[xmax * ymax * 3]; p++)
	*p = light? 255:0;
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
