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

#ifdef _SVG
#include <cairo/cairo-svg.h>
#endif
 
#ifdef _PDF
#include <cairo/cairo-pdf.h> 
#endif

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

#define PPI 72.0  /* pixels per inch */
#define NCOLOR 256 /* number of colors */
#define MAXVERT 1000

static cairo_surface_t *surface;
static cairo_t *cr;

static bool light, fill;
static int color_table[NCOLOR][3], nx, ny;
static float magnify;

static cairo_status_t cr_fwrite(void *closure, 
				const unsigned char *data, unsigned int length)
{
    if (length == fwrite(data,1,length,pltout))
	return CAIRO_STATUS_SUCCESS;
    return CAIRO_STATUS_WRITE_ERROR;
}

void opendev (int argc, char* argv[])
/*< open >*/
{
    int value;
    char newpath[60], *color;

    dev.txfont = DEFAULT_HARDCOPY_FONT;
    dev.txprec = DEFAULT_HARDCOPY_PREC;
    dev.brake = BREAK_IGNORE;

    /* control routines */
    dev.reset = crreset;
    dev.erase = crerase;
    dev.close = crclose;

    dev.attributes = crattr;
    dev.vector = crvector;
    dev.plot = crplot;
    dev.area = crarea;

    dev.pixels_per_inch = PPI;
    nx = VP_STANDARD_HEIGHT * dev.pixels_per_inch; 
    ny = VP_SCREEN_RATIO * VP_STANDARD_HEIGHT * dev.pixels_per_inch;
    /* image size */

    dev.aspect_ratio = 1.;

    sf_getfloat ("aspect", &dev.aspect_ratio);
    /* aspect ratio */
    sf_getfloat ("ppi", &dev.pixels_per_inch);
    /* pixels per inch */
    
    magnify = dev.pixels_per_inch / PPI;
    nx *= magnify;
    ny *= (1.0 / dev.aspect_ratio) * magnify;

    sf_getint ("n1", &nx);    
    sf_getint ("n2", &ny); 
    /* image size */

    dev.need_end_erase = true;
    dev.smart_clip= true; 
    dev.num_col = NCOLOR;

    dev.xmax = nx-1;
    dev.ymax = ny-1;

    if (NULL == (color = sf_getstring("bgcolor"))) color="light";
    /* background color (black,white,light,dark) 
       'light' and 'dark' cause the background to be transparent (in PNG and GIF) */

    fill = (bool)  ('l' != color[0] && 'd' != color[0]);
    light = (bool) ('b' != color[0] && 'd' != color[0]);
 
#ifdef _PDF
    surface = cairo_pdf_surface_create_for_stream(cr_fwrite, NULL, nx, ny);
#endif

#ifdef _PNG
    surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, nx, ny);
#endif
	    
#ifdef _SVG
    surface = cairo_svg_surface_create_for_stream(cr_fwrite, NULL, nx, ny);
#endif

    cr = cairo_create (surface);
    cairo_set_fill_rule (cr, CAIRO_FILL_RULE_EVEN_ODD);

    cairo_set_line_width (cr, magnify);
    cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND); 

    if (isatty(fileno(pltout)))
    {
	sprintf (newpath, "%s", "image_file");
	pltout = fopen (newpath, "wb");
	if (pltout == NULL)
	    ERR (FATAL, name, "can't open file %s\n", newpath);
    }

    for (value = 0; value < NCOLOR; value++)
    {
	color_table[value][0] = -1;
    } 
}

void crreset (void)
/*< reset >*/
{
    int value;
    const int red[]   = {255, 255,   0,   0, 255, 255,   0, 0};
    const int green[] = {255, 255, 255, 255,   0,   0,   0, 0};
    const int blue[]  = {255,   0, 255,   0, 255,   0, 255, 0};

    for (value = 0; value < 8; value++)
    {
	if (light) {
	    color_table[value][0] = red[value];
	    color_table[value][1] = green[value];
	    color_table[value][2] = blue[value];
	} else {
	    color_table[value][0] = 255-red[value];
	    color_table[value][1] = 255-green[value];
	    color_table[value][2] = 255-blue[value];
	}
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

#ifdef _PNG
    cairo_surface_write_to_png_stream (surface,cr_fwrite,NULL);
#else
    cairo_show_page(cr);
#endif

    called=true;
}

static void fillall(void)
{
    cairo_new_path(cr);
    cairo_rectangle (cr,0,0,dev.xmax+1,dev.ymax+1);
    cairo_set_source_rgb (cr, 
			  color_table[0][0]/255.0, 
			  color_table[0][1]/255.0, 
			  color_table[0][2]/255.0);
    cairo_fill(cr);
}

void crerase (int command)
/*< erase >*/
{
    switch (command)
    {
	case ERASE_START:
	    if (fill) fillall();
	    break;
	case ERASE_MIDDLE:
	    cr_write();
	    if (fill) fillall();
	    break;
	case ERASE_END:
	    cr_write();
	    cairo_surface_destroy(surface);
	    cairo_destroy(cr);
	    fclose(pltout); 
	    break;
	case ERASE_BREAK:
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
    int x0, y0, width, height;

    switch (command)
    {
	case SET_COLOR:
	    cairo_set_source_rgb (cr, 
				  color_table[value][0]/255.0, 
				  color_table[value][1]/255.0, 
				  color_table[value][2]/255.0);
	    break;
	case SET_COLOR_TABLE:
	    color_table[value][0] = v1;
	    color_table[value][1] = v2;
	    color_table[value][2] = v3;
	    break;
	case SET_WINDOW:
	    x0 = value;
	    width = SF_MAX(v2-value+1,0);
	    y0 = dev.ymax-v3;
	    height = SF_MAX(v3-v1+1,0);

	    cairo_rectangle (cr,x0,y0,width,height);
	    cairo_clip(cr);
	default:
	    break;
    }
}

void crplot(int x, int y, int draw)
/*< plot >*/
{
    static int oldx, oldy;
 
    if (draw) {
	cairo_move_to (cr,oldx,dev.ymax-oldy);
	cairo_line_to (cr,x,dev.ymax-y);
	cairo_stroke(cr);
    } else {
	cairo_move_to (cr,x,dev.ymax-y);
	dev.lost = 0;
    }
    
    oldx=x;
    oldy=y;
}

void crarea (int npts, struct vertex *head)
/*< area >*/
{
    int i;

    cairo_new_path(cr);
    cairo_move_to (cr,head->x,dev.ymax - head->y);
    head = head->next;
    for (i = 1; i < npts; i++)
    {
	cairo_line_to (cr,head->x,dev.ymax - head->y);
	head = head->next;
    }
    cairo_fill(cr);
}



void crvector (int x1, int y1, int x2, int y2, int nfat, int dashon)
/*< device-dependent vector >*/
{
    if (nfat < 0)
	return;
    
    if (nfat > 1) cairo_set_line_width (cr, nfat*magnify);
    
    genvector(x1,y1,x2,y2,0,dashon);
    
    if (nfat > 1) cairo_set_line_width (cr, magnify);
}

