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
 *  source file:   ./filters/utilities/clip.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 *
 */

#include "solve.h"

extern int      xwmin, xwmax, ywmin, ywmax;

#define code(x,y) (x<xwmin?1:(x>xwmax?2:0))|(y<ywmin?4:(y>ywmax?8:0))

int clip (int *x1, int *y1, int *x2, int *y2)
/*< Cohen-Sutherland Clipping routine >*/
{
register int    c1, c2, temp;
int             swap;
    c1 = code (*x1, *y1);
    c2 = code (*x2, *y2);
    swap = 0;
    if (!(c1 || c2))
	return (0);		/* line completely in bounds */
    while (c1 | c2)
    {
	if (c1 & c2)
	    return (1);		/* line completely out of bounds */
	if (!c1)		/* interchange endpoints */
	{
	    temp = *x1;
	    *x1 = *x2;
	    *x2 = temp;
	    temp = *y1;
	    *y1 = *y2;
	    *y2 = temp;
	    temp = c1;
	    c1 = c2;
	    c2 = temp;
	    swap = ~swap;
	}
	if (c1 < 4)		/* move endpoint in x */
	{
	    temp = (c1 & 2 ? xwmax : xwmin);
	    *y1 = solve (temp, *x1, *y1, *x2, *y2);
	    *x1 = temp;
	}
	else			/* move endpoint in y */
	{
	    temp = (c1 & 8 ? ywmax : ywmin);
	    *x1 = solve (temp, *y1, *x1, *y2, *x2);
	    *y1 = temp;
	}
	c1 = code (*x1, *y1);
    }
    if (swap)			/* put endpoints in order */
    {
	temp = *x1;
	*x1 = *x2;
	*x2 = temp;
	temp = *y1;
	*y1 = *y2;
	*y2 = temp;
    }
    return (0);
}
