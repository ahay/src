/* vplot filter for LibGD. */
/*
  Copyright (C) 2009 The University of Texas at Austin
  
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

#include <gd.h>

#include "../include/attrcom.h"
#include "../include/extern.h"
#include "../include/params.h"
#include "../include/err.h"
#include "../include/erasecom.h"
#include "../include/closestat.h"
#include "../include/enum.h"

#include "../genlib/genpen.h"
#include "../utilities/util.h"

#include "dovplot.h"
#include "init_vplot.h"

#include "gdpen.h"

char            name[] = "gdpen";
#include "gddoc.h"

#include "_device.h"

#define PPI 100.0  /* pixels per inch */
#define NCOLOR 256 /* number of colors */
#define MAXVERT 1000

static gdImagePtr image, oldimage;
static bool default_out = true;
static int color_table[NCOLOR], gdcolor, black, delay;
static char *image_type;

void opendev (int argc, char* argv[])
/*< open >*/
{
    float pixels_per_inch, aspect_ratio;
    char newpath[60];
    int value;

    dev.txfont = DEFAULT_HARDCOPY_FONT;
    dev.txprec = DEFAULT_HARDCOPY_PREC;
    dev.brake = BREAK_IGNORE;

    /* control routines */
    dev.reset = gdreset;
    dev.erase = gderase;
    dev.close = gdclose;

    dev.area = gdarea;
    dev.attributes = gdattr;
    dev.plot = gdplot;

    pixels_per_inch = dev.pixels_per_inch = PPI;
    dev.xmax = VP_STANDARD_HEIGHT * dev.pixels_per_inch;
    dev.ymax = VP_SCREEN_RATIO * dev.xmax;
    aspect_ratio = dev.aspect_ratio = 1.;

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

    dev.need_end_erase = true;
    dev.smart_clip= true; 
    dev.num_col = NCOLOR;

    image = gdImageCreate(dev.xmax,dev.ymax);
    black = gdImageColorAllocate(image, 0, 0, 0);

    default_out = (bool) isatty(fileno(pltout));

    if (default_out)
    {
	sprintf (newpath, "%s", "image_file");
	pltout = fopen (newpath, "wb");
	if (pltout == NULL)
	    ERR (FATAL, name, "can't open file %s\n", newpath);
    }

    for (value = 0; value < NCOLOR; value++)
    {
	color_table[value] = -1;
    }

    if (NULL == (image_type = sf_getstring("type"))) 
	image_type = "png";
    /* image type (png, jpeg, gif) */

    if (!sf_getint("delay",&delay)) delay=10;
    /* GIF animation delay (if type=="gif") */
}

void gdreset (void)
/*< reset >*/
{
    int value, color;

    /* reset color table */
    color = color_table[0];
    if (-1 != color) gdImageColorDeallocate(image, color);
    color_table[0] = black;

    for (value = 1; value < 8; value++)
    {
	color = color_table[value];
	if (-1 != color) gdImageColorDeallocate(image, color);
	color_table[value] = gdImageColorAllocate(image, 
						  MAX_GUN * ((value & 2) / 2),
						  MAX_GUN * ((value & 4) / 4),
						  MAX_GUN * ((value & 1) / 1));
    }
 
    for (value = 8; value < NCOLOR; value++)
    {
	color = color_table[value];
	if (-1 != color) {
	    gdImageColorDeallocate(image, color);
	    color_table[value] = -1;
	}
    }
}

static void gd_write (void)
{
    static bool called=false;

    switch (image_type[0]) {
	case 'j':
	case 'J':
	    if(called) return;
	    gdImageJpeg(image, pltout, -1);
	    break;
	case 'g':
	case 'G':
	    if (called) {
		gdImageGifAnimAdd(image, pltout, 0, 0, 0, delay, 1, oldimage);
	    } else {
		gdImageGifAnimBegin(image, pltout, 1, 0);
		gdImageGifAnimAdd(image, pltout, 0, 0, 0, delay, 1, NULL);
	    }
	    oldimage = image;
	    image = gdImageCreate(dev.xmax,dev.ymax);
	    gdImagePaletteCopy(image, oldimage);
	    break;
	case 'p':
	case 'P':
	default:
	    if(called) return;
	    gdImagePng(image, pltout);
	    break;
    }

    called=true;
}

void gderase (int command)
/*< erase >*/
{
    switch (command)
    {
	case ERASE_START:
	    gdImageFilledRectangle(image, 0, 0, dev.xmax, dev.ymax, black);
	    break;
	case ERASE_MIDDLE:
	    gd_write();
	    gdImageFilledRectangle(image, 0, 0, dev.xmax, dev.ymax, black);
	    break;
	case ERASE_END:
	    gd_write();
	    if (image_type[0]=='g' || image_type[0]=='G') 
		gdImageGifAnimEnd(pltout);
	    fclose(pltout);
	    break;
	case ERASE_BREAK:
	    gd_write();
	    break;
	default:
	    break;
    }
}

void gdclose (int status)
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

void gdattr (int command, int value, int v1, int v2, int v3)
/*< attr >*/
{
    int xmin, ymin, xmax, ymax;

    switch (command)
    {
	case SET_COLOR:
	    gdcolor = color_table[value];
	    break;
	case SET_COLOR_TABLE:
	    color_table[value] = gdImageColorAllocate(image, v1, v2, v3);
	    break;
	case SET_WINDOW:
	    xmin = value;
	    ymin = dev.ymax - v3;
	    xmax = v2;
	    ymax = dev.ymax - v1;
	    gdImageSetClip(image, xmin, ymin, xmax, ymax);
	    break;
	default:
	    break;
    }
}

void gdplot(int x, int y, int draw)
/*< plot >*/
{
    static int oldx = 0, oldy = 0;

    if (draw)
    {
	gdImageLine(image, oldx, dev.ymax - oldy, x, dev.ymax - y, gdcolor);
    } else {
	dev.lost = 0;
    }

    oldx = x;
    oldy = y;
}

void gdarea (int npts, struct vertex *head)
/*< area >*/
{
    int i;
    gdPoint vlist[MAXVERT];

    /* translate data structures */
    for (i = 0; i < npts && i < MAXVERT; i++)
    {
	vlist[i].x = head->x;
	vlist[i].y = dev.ymax - head->y;
	head = head->next;
    }

    gdImageFilledPolygon(image, vlist, npts, gdcolor);
}
