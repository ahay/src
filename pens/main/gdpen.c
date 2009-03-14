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

#include "../genlib/genpen.h"

#include "dovplot.h"
#include "init_vplot.h"

#include "gdpen.h"

char            name[] = "gdpen";
#include "gddoc.h"

static gdImagePtr image;
static bool     default_out = true;

struct device   dev =
{

    /* control routines */
    gdopen,		/* open */
    gdreset,		/* reset */
    genmessage,		/* message */
    gderase,		/* erase */
    gdclose,		/* close */
    
    /* high level output */
    genvector,		/* vector */
    genmarker,		/* marker */
    gentext,		/* text */
    genpatarea,		/* area */
    genraster,		/* raster */
    genpoint,		/* point */
    gdattr,		/* attributes */
    
    /* input */
    gen_dovplot,            /* reader */
    nullgetpoint,		/* getpoint */
    nullinteract,		/* interact */
    
    /* low level output */
    gdplot,		/* plot */
    nullclose,		/* startpoly */
    nullmidpoly,		/* midpoly */
    nullclose		/* endpoly */
};

#define NCOLOR 256 /* number of colors */

void gdopen (int argc, char* argv[])
/*< open >*/
{
    char newpath[60];

    txfont = DEFAULT_HARDCOPY_FONT;
    txprec = DEFAULT_HARDCOPY_PREC;

    if (!sf_getint ("n1", &dev_xmax)) dev_xmax = VP_STANDARD_HEIGHT * pixels_per_inch;
    if (!sf_getint ("n2", &dev_ymax)) dev_ymax = VP_SCREEN_RATIO * dev_xmax;
    /* image size */

    /* Why not VP_STANDARD_HEIGHT and VP_SCREEN_RATIO ??? */

    dev_xmin = 0;
    dev_ymin = 0;

    need_end_erase = true;
    smart_clip = false;
    num_col = NCOLOR;

    image = gdImageCreate(dev_xmax,dev_ymax);

    default_out = (bool) isatty(fileno(pltout));

    if (default_out)
    {
	sprintf (newpath, "%s", "image_file");
	pltout = fopen (newpath, "w");
	if (pltout == NULL)
	    ERR (FATAL, name, "can't open file %s\n", newpath);
    }
}
