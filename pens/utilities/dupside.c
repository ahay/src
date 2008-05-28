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
 *  source file:   ./filters/utilities/dupside.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include "../include/vertex.h"
/*^*/

int dupside (register struct vertex *base)
/*< Determine if other sides in the polygon are
 * identical to the side specified by the
 * vertices v and v->b. >*/
{
register struct vertex *v;
register int    x1, x2, y1, y2;

    x1 = base->x;
    x2 = base->last->x;
    y1 = base->y;
    y2 = base->last->y;
    v = base->next;
    do
    {
	if (x1 == v->x && y1 == v->y && x2 == v->last->x && y2 == v->last->y)
	    return (1);
	if (x2 == v->x && y2 == v->y && x1 == v->last->x && y1 == v->last->y)
	    return (1);
	v = v->next;
    }
    while (v != base);
    return (0);
}
