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
 *  source file:   ./filters/utilities/vecoutline.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include <stdio.h>
#include "../include/params.h"
#include "../include/pat.h"
#include "../include/vertex.h"
#include "../include/extern.h"

#include "dupside.h"

void vecoutline (struct vertex  *head)
/*< Draw the outline of the polygon pointed to by 'head'.  If any side of
 * the polygon is identical to any other (that is they have identical
 * endpoints), then neither line is drawn.  This allows among things a
 * doughnut to be defined by a single polygon without the line that
 * connects the inner and outer being plotted. >*/
{
    register int    xlast, ylast;
    register struct vertex *v;

    xlast = head->last->x;
    ylast = head->last->y;
    v = head;
    do
    {
	if (!dupside (v))
	    dev.vector (v->x, v->y, xlast, ylast, afat, 0);
	xlast = v->x;
	ylast = v->y;
    } while ((v = v->next) != head);
}
