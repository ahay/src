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
