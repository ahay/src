/*
 * Copyright 1987 the Board of Trustees of the Leland Stanford Junior
 * University. Official permission to use this software is included in
 * the documentation. It authorizes you to use this file for any
 * non-commercial purpose, provided that this copyright notice is not
 * removed and that any modifications made to this file are commented
 * and dated in the style of my example below.
 */

/*
 *
 *  source file:   ./filters/genlib/genraster.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger Feb 16 1988
 *	Make number of arguments consistent between smart and dumb forms
 *	of dev.raster.
 */

#include <stdio.h>
#include "../include/err.h"
#include "../include/extern.h"
#include "../include/enum.h"
#include "../include/attrcom.h"

extern bool      overlay;
extern int  cur_color, need_devcolor;

void genraster (int count, int out_of, int xpos, int ypos, int length, int orient, 
		unsigned char **raster, int dummy1, int dummy2)
/*< device-independent raster >*/
{
int             ii, sign=0, xy=0, xrpos=0, yrpos=0;
int             color, start;
static int      cur_color_save;

    switch (orient)
    {
    case 0:
	xrpos = xpos;
	yrpos = ypos;
	sign = 1;
	xy = 0;
	break;
    case 1:
	xrpos = ypos;
	yrpos = xpos;
	sign = -1;
	xy = 1;
	break;
    case 2:
	xrpos = xpos;
	yrpos = ypos;
	sign = -1;
	xy = 0;
	break;
    case 3:
	xrpos = ypos;
	yrpos = xpos;
	sign = 1;
	xy = 1;
	break;
    }

    start = xrpos;
    color = raster[0][0];

    if (count == 0)
    {
	/*
	 * First time remember the color so we can restore it the last time. 
	 */
	cur_color_save = cur_color;
    }

    for (ii = 0; ii <= length; ii++)
    {
	if (ii == length || raster[0][ii] != color)
	{
	    if (!((overlay == 1) && (color == 0)))
	    {
		if (cur_color != color || need_devcolor)
		{
		    cur_color = color;
		    dev.attributes (SET_COLOR, cur_color, 0, 0, 0);
		    need_devcolor = NO;
		}
		if (xy)
		    dev.vector (yrpos, start, yrpos, xrpos + sign * (ii - 1), 0, 0);
		else
		    dev.vector (start, yrpos, xrpos + sign * (ii - 1), yrpos, 0, 0);
	    }
	    if (ii != length)
	    {
		color = raster[0][ii];
		start = xrpos + sign * ii;
	    }
	}
    }

    if (count == out_of - 1 && cur_color != cur_color_save)
    {
	cur_color = cur_color_save;
	need_devcolor = YES;
    }
}
