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
 *  source file:   ./filters/genlib/genraster1.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger Feb 16 1988
 *	Make number of arguments consistent between smart and dumb forms
 *	of dev.raster.
 */

#include <stdio.h>
#include <stdlib.h>

#include "../include/extern.h"
#include "../include/err.h"
#include "../include/enum.h"
#include "../include/attrcom.h"
#include "../include/params.h"

#define	MAXMEM 200000
#define OFF 0
#define ON 1

extern bool      overlay;
extern int cur_color, need_devcolor, num_col_8;
/*
 * Less than or equal to this, better to use point mode 
 */
int             break_point = 4;

/*
 * A more efficient genraster version for devices that can draw a
 * string of points quickly. Break the plot up by color, and by
 * points and vectors.
 */

void genraster1 (int count, int out_of, int xpos, int ypos, int length, int orient, 
		 unsigned short *raster, int dummy1, int dummy2, int byte2)
{
int             ii, jj, yy, kk, ll;
int             xstart=0, state;
static int      cur_color_save;
static int      num_lines;
static unsigned char *array;
static int      ylength, ystart, jstart;
static int      color_used[MAX_COL + 1];
int             xsign=0, ysign=0, xy=0, xrpos=0, yrpos=0;

    switch (orient)
    {
    case 0:
	xrpos = xpos;
	yrpos = ypos;
	xsign = 1;
	ysign = 1;
	xy = 0;
	break;
    case 1:
	xrpos = ypos;
	yrpos = xpos;
	xsign = -1;
	ysign = 1;
	xy = 1;
	break;
    case 2:
	xrpos = xpos;
	yrpos = ypos;
	xsign = -1;
	ysign = -1;
	xy = 0;
	break;
    case 3:
	xrpos = ypos;
	yrpos = xpos;
	xsign = 1;
	ysign = -1;
	xy = 1;
	break;
    }

    if (count == 0)
    {
	/*
	 * First time remember the color so we can restore it the last time. 
	 */
	cur_color_save = cur_color;

	/*
	 * Also find out how many lines we can do at once. 
	 */
	num_lines = MAXMEM / length;
	if (num_lines < 1)
	    num_lines = 1;
	array = (unsigned char *) malloc ((unsigned) num_lines * length * sizeof (unsigned char));
	/*
	 * See whether we need to do color 0 or not 
	 */
	jstart = 0;
	if (overlay == 1)
	    jstart = 1;
    }

    if (ylength == 0)
    {
	/*
	 * Just starting a block. Remember where it is. 
	 */
	ystart = yrpos;
	for (ii = 0; ii <= num_col_8; ii++)
	{
	    color_used[ii] = 0;
	}
    }

/*
 * Save it.
 */
    for (ii = 0; ii < length; ii++)
    {
	array[length * ylength + ii] = raster[ii];
	color_used[raster[ii]] = 1;
    }
    ylength++;

    if (ylength >= num_lines || count == out_of - 1)
    {
/*
 * Plot it. Loop by color
 */

	need_devcolor = NO;

	for (kk = 0; kk < 2; kk++)
	{
/*
 * For maximum efficiency of output, better to have the kk loop
 * inside the jj loop. However, while watching on the screen better
 * to have "most important" stuff come out first.
 */
	    for (jj = jstart; jj < num_col_8; jj++)
	    {
		if (color_used[jj] == 0)
		    continue;

		cur_color = jj;
		dev.attributes (SET_COLOR, cur_color, 0, 0, 0);

		for (yy = 0; yy < ylength; yy++)
		{
		    state = OFF;
		    for (ii = 0; ii <= length; ii++)
		    {
			if (ii != length && array[length * yy + ii] == jj && state == OFF)
			{
			    xstart = xrpos + xsign * ii;
			    state = ON;
			    continue;
			}
			if ((ii == length || array[length * yy + ii] != jj) && state == ON)
			{
			    switch (kk)
			    {
			    case 0:
				if (xsign * (xrpos + xsign * ii - xstart) > break_point)
				{
				    if (xy)
					dev.vector (ystart - ysign * yy, xstart, ystart - ysign * yy, xrpos + xsign * (ii - 1), 0, 0);
				    else
					dev.vector (xstart, ystart - ysign * yy, xrpos + xsign * (ii - 1), ystart - ysign * yy, 0, 0);
				}
				break;
			    case 1:
				if (xsign * (xrpos + xsign * ii - xstart) <= break_point)
				{
				    for (ll = 0; ll < xsign * (xrpos + xsign * ii - xstart); ll++)
				    {
					if (xy)
					    dev.point (ystart - ysign * yy, xstart + xsign * ll);
					else
					    dev.point (xstart + xsign * ll, ystart - ysign * yy);
				    }
				}
				break;
			    }
			    state = OFF;
			}
		    }
		}

	    }
	}
	ylength = 0;
	if (count == out_of - 1)
	{
	    free ((char *) array);
	    if (cur_color != cur_color_save)
	    {
		cur_color = cur_color_save;
		need_devcolor = YES;
	    }
	}
    }
}
