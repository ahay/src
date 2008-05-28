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
 *  source file:   ./filters/genlib/vecarea.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include <stdio.h>
#include <stdlib.h>

#include "../include/pat.h"
#include "../include/vertex.h"
#include "../include/params.h"
#include "../include/extern.h"

#include "../utilities/util.h"

void vecarea (int npts, struct vertex  *head)
/*< device-independent area >*/
{
register int    x, y, i;
register int    xstr, xend, ystr, yend;
int             skip;
int             ncross;
int             vminx, vmaxx, vminy, vmaxy;
struct vertex  *xhead, *yhead, *v;
int            *crosses;

    /*
     * allocate storage for scan line cross points 
     */
    crosses = (int *) malloc ((unsigned) npts * sizeof (int));

    /*
     * double link the vertices. (head) is set to the node with the maximum
     * x-value so that intersect() will not eliminate 'head' while casting
     * off vertices. 
     */
    vminx = head->x;
    vmaxx = head->x;
    vminy = head->y;
    vmaxy = head->y;
    xhead = head;
    yhead = head;

    v = head;
    for (i = 0; i < npts; i++)
    {
	if (v->x > vmaxx)
	{
	    vmaxx = v->x;
	    xhead = v;
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

    if ((pat[ipat] .ydim > 0) || (pat[ipat] .xdim == 1))
    {
	/* stretch polygon in y-direction */
	v = yhead;
	do
	{
	    v->y = 2 * (v->y) + 1;
	    v = v->next;
	} while (v != yhead);

	skip = (pat[ipat] .xdim == 1) ? 1 : pat[ipat] .ydim;
	for (y = vminy; y <= vmaxy; y += skip)
	{
	    ncross = intersect (2 * y, crosses, yhead, 1);
	    sort (crosses, ncross);
	    for (i = 0; i < ncross; i += 2)
	    {
		xstr = crosses[i];
		xend = crosses[i + 1];
		dev.vector (xstr, y, xend, y, 0, 0);
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

    if ((pat[ipat] .xdim > 1) && (pat[ipat] .ydim > 1))
    {
	/*
	 * expand in x 
	 */
	v = xhead;
	do
	{
	    v->x = 2 * v->x + 1;
	    v = v->next;
	} while (v != xhead);

	skip = pat[ipat] .xdim;
	for (x = vminx; x <= vmaxx; x += skip)
	{
	    ncross = intersect (2 * x, crosses, xhead, 0);
	    sort (crosses, ncross);
	    for (i = 0; i < ncross; i += 2)
	    {
		ystr = crosses[i];
		yend = crosses[i + 1];
		dev.vector (x, ystr, x, yend, 0, 0);
	    }
	}

	/*
	 * shrink in x 
	 */
	v = xhead;
	do
	{
	    v->x = ((v->x - 1) / 2);
	    v = v->next;
	} while (v != xhead);
    }
    free ((char *) crosses);
}
