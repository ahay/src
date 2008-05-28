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
 *  source file:   ./filters/genlib/genpatarea.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger, Feb 16 1988
 *	Make number of arguments to dev.attributes and dev.raster consistent.
 */

#include <stdio.h>
#include <stdlib.h>

#include "../include/pat.h"
#include "../include/vertex.h"
#include "../include/params.h"
#include "../include/extern.h"
#include "../include/enum.h"
#include "../include/attrcom.h"
#include "../include/err.h"

#include "../utilities/util.h"


/* Make the damn modulo function work for negative numbers */
#define MODULO(A,B)	(((A)>=0)?((A)%(B)):(((A)%(B)+(B))%(B)))

/*
 * This routine turns polygons filled with a pattern into calls to
 * dev.raster. (one line at a time).
 */
extern int      need_devcolor, cur_color;

void genpatarea (int npts, struct vertex  *head)
/*< patarea >*/
{
    register int    y, i, ii;
    register int    xstr, xend;
    int             ncross;
    int             vminx, vmaxx, vminy, vmaxy;
    struct vertex  *yhead, *v;
    int            *crosses;
    unsigned char  **rasline;
    static int      cur_color_save;
    int             color;
 
/*
 * Save color so we can restore it the last time.
 */
    cur_color_save = cur_color;

    /*
     * allocate storage for scan line cross points 
     */
    crosses = (int *) malloc ((unsigned) npts * sizeof (int));

    /*
     * allocate storage for raster line 
     */

    rasline = (unsigned char **) malloc (sizeof (unsigned char*));
    if (rasline == NULL)
	ERR (FATAL, name, "cannot alloc memory to load raster image");
    rasline[0] = (unsigned char *) malloc ((xwmax - xwmin + 1) *  sizeof (unsigned char));
    if (rasline[0] == NULL)
	ERR (FATAL, name, "cannot alloc memory to load raster image");

    /*
     * double link the vertices. (head) is set to the node with the maximum
     * x-value so that intersect() will not eliminate 'head' while casting
     * off vertices. 
     */
    vminx = head->x;
    vmaxx = head->x;
    vminy = head->y;
    vmaxy = head->y;
    yhead = head;

    v = head;
    for (i = 0; i < npts; i++)
    {
	if (v->x > vmaxx)
	{
	    vmaxx = v->x;
	}
	if (v->y > vmaxy)
	{
	    vmaxy = v->y;
	    yhead = v;
	}
	if (v->x < vminx)
	    vminx = v->x;
	if (v->y < vminy)
	    vminy = v->y;
	v++;
    }

    if (vmaxx > xwmax)
	vmaxx = xwmax;
    if (vminx < xwmin)
	vminx = xwmin;
    if (vmaxy > ywmax)
	vmaxy = ywmax;
    if (vminy < ywmin)
	vminy = ywmin;

    if ((pat[ipat] .ydim > 0) && (pat[ipat] .xdim > 0))
    {
	/* stretch polygon in y-direction */
	v = yhead;
	do
	{
	    v->y = 2 * (v->y) + 1;
	    v = v->next;
	} while (v != yhead);

	for (y = vminy; y <= vmaxy; y++)
	{
	    ncross = intersect (2 * y, crosses, yhead, 1);
	    sort (crosses, ncross);
	    for (i = 0; i < ncross; i += 2)
	    {
		xstr = crosses[i];
		xend = crosses[i + 1];

		if (xstr < xwmin && xend < xwmin)
		    continue;
		if (xstr > xwmax && xend > xwmax)
		    continue;

		if (xstr < xwmin)
		    xstr = xwmin;
		if (xend > xwmax)
		    xend = xwmax;

		if (pat[ipat] .xdim == 1)
		{
/* Faster to fill it with one vector */
		    color =
		     pat[ipat] .patbits[
					((pat[ipat] .ydim - 1) - (MODULO (y, pat[ipat] .ydim))) * pat[ipat] .xdim
		     ];
		    if (cur_color != color || need_devcolor)
		    {
			cur_color = color;
			dev.attributes (SET_COLOR, cur_color, 0, 0, 0);
			need_devcolor = NO;
		    }
		    dev.vector (xstr, y, xend, y, 0, 0);
		}
		else
		{
		    for (ii = xstr; ii <= xend; ii++)
		    {
			rasline[0][ii - xstr] =
			 pat[ipat] .patbits[
					    ((pat[ipat] .ydim - 1) - (MODULO (y, pat[ipat] .ydim))) * pat[ipat] .xdim
					    + MODULO (ii, pat[ipat] .xdim)
			 ];
		    }
		    if (xstr <= xend)
			dev.raster (0, 1, xstr, y, xend - xstr + 1, 0, rasline, 0, 0);
		}
	    }
	}
	/* shrink in y */
	v = yhead;
	do
	{
	    v->y = ((v->y - 1) / 2);
	    v = v->next;
	} while (v != yhead);
    }

    free (*rasline);
    free (rasline);
    free (crosses);

    if (cur_color != cur_color_save)
    {
	cur_color = cur_color_save;
	need_devcolor = YES;
    }
}
