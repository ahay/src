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
 *  source file:   ./filters/utilities/intersect.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include	"../include/vertex.h"

#include "solve.h"

int intersect (int x, int *crosses, struct vertex  *head, int scany)
/*< intersect >*/
{
register struct vertex *v;
int             ncross, y;

    ncross = 0;
    v = head;
    if (scany)
    {
	do
	{
	    if (v->y > x && v->last->y > x)
		continue;
	    if (v->y < x && v->last->y < x)
		continue;

	    y = solve (x, v->y, v->x, v->last->y, v->last->x);
	    crosses[ncross++] = y;
	} while ((v = v->next) != head);
    }
    else
    {
	do
	{
	    if (v->x > x && v->last->x > x)
		continue;
	    if (v->x < x && v->last->x < x)
		continue;

	    y = solve (x, v->x, v->y, v->last->x, v->last->y);
	    crosses[ncross++] = y;
	} while ((v = v->next) != head);
    }
    return (ncross);
}
