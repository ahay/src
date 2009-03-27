/* vplot filter for Cairo Graphics. */
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

#include <cairo/cairo.h>
#include <cairo/cairo-svg.h> 
#include <cairo/cairo-pdf.h> 

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

#include "crpen.h"

char            name[] = "crpen";
#include "crdoc.h"

#include "_device.h"

#define PPI 100.0  /* pixels per inch */
#define NCOLOR 256 /* number of colors */
#define MAXVERT 1000

static cairo_surface_t *surface;
static cairo_t *cr;

static float color_table[NCOLOR][3], nx, ny;
static int type;

static cairo_status_t cr_fwrite(void *closure, const unsigned char *data, unsigned int length)
{
    if (length != fwrite(data,1,length,pltout)) ERR(FATAL, name, "trouble writing");
    return 0;
}

void opendev (int argc, char* argv[])
/*< open >*/
{
    float pixels_per_inch, aspect_ratio;
    char newpath[60], *image_type;

    dev.txfont = DEFAULT_HARDCOPY_FONT;
    dev.txprec = DEFAULT_HARDCOPY_PREC;
    dev.brake = BREAK_IGNORE;

    /* control routines */
    dev.reset = crreset;
    dev.erase = crerase;
    dev.close = crclose;

    dev.attributes = crattr;
    dev.plot = crplot;

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

    nx = dev.xmax+1;
    ny = dev.ymax+1;

    if (NULL == (image_type = sf_getstring("type"))) 
	image_type = "png";
    /* image type (png, pdf, svg) */

    switch (image_type[0]) {
	case 'p':
	case 'P':
	default:
	    if ('d'==image_type[1] || 'D'==image_type[1]) {
		type=1;
		surface = cairo_pdf_surface_create_for_stream(cr_fwrite, NULL, nx, ny);
	    } else {
		type=0;
		surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, nx, ny);
	    } 
	    break;
	case 's':
	case 'S':
	    type=2;
	    surface = cairo_svg_surface_create_for_stream(cr_fwrite, NULL, nx, ny);
	    break;
    }

    cr = cairo_create (surface);

    if (isatty(fileno(pltout)))
    {
	sprintf (newpath, "%s", "image_file");
	pltout = fopen (newpath, "wb");
	if (pltout == NULL)
	    ERR (FATAL, name, "can't open file %s\n", newpath);
    }
}

static void cr_clear(void)
{
    cairo_set_source_rgb (cr, 
			  color_table[0][0], 
			  color_table[0][1], 
			  color_table[0][2]);
    cairo_rectangle (cr, 0, 0, nx, ny);
    cairo_fill(cr);
}

void crreset (void)
/*< reset >*/
{
    int value;

    cr_clear();
    color_table[0][0] =0.0;
    color_table[0][1] =0.0;
    color_table[0][2] =0.0;
    for (value = 1; value < 8; value++)
    {
	color_table[value][0] = ((value & 2) / 2);
	color_table[value][1] = ((value & 4) / 4);
	color_table[value][2] = ((value & 1) / 1);
    }
    for (value = 8; value < NCOLOR; value++)
    {
	color_table[value][0] = -1;
    }   
}


static void cr_write (void)
{
    static bool called=false;

    if(called) return;

    if (type) {
	cairo_show_page(cr);
    } else {
	cairo_surface_write_to_png_stream (surface,cr_fwrite,NULL);
    }
    
    cairo_surface_destroy(surface);
    cairo_destroy(cr);

    called=true;
}


void crerase (int command)
/*< erase >*/
{
    switch (command)
    {
	case ERASE_START:
	    cr_clear();
	    break;
	case ERASE_MIDDLE:
	    cairo_stroke(cr);
	    cr_write();
	    cr_clear();
	    break;
	case ERASE_END:
	    cairo_stroke(cr);
	    cr_write();
	    fclose(pltout);
	    break;
	case ERASE_BREAK:
	    cairo_stroke(cr);
	    cr_write();
	    break;
	default:
	    break;
    }
}

void crclose (int status)
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

void crattr (int command, int value, int v1, int v2, int v3)
/*< attr >*/
{
    switch (command)
    {
	case SET_COLOR:
	    cairo_stroke(cr);
	    cairo_set_source_rgb (cr, 
				  color_table[value][0], 
				  color_table[value][1], 
				  color_table[value][2]);
	    break;
	case SET_COLOR_TABLE:
	    color_table[value][0] = v1/255.0;
	    color_table[value][1] = v2/255.0;
	    color_table[value][2] = v3/255.0;
	    break;
	default:
	    break;
    }
}

void crplot(int x, int y, int draw)
/*< plot >*/
{
    if (draw) {
	cairo_line_to (cr,x,dev.ymax-y);
    } else {
	cairo_move_to (cr,x,dev.ymax-y);
	dev.lost = 0;
    }
}

