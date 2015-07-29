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
 *  source file:   ./filters/genlib/genarea.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include <stdio.h>
#include "../include/vertex.h"
#include "../include/params.h"
#include "../include/extern.h"

#include "polysubs.h"
#include "polyfix.h"

extern int      Allgone;

void genarea (int npts, struct vertex  *head)
/*< Device Independent Polygon treatment.
 * Do a first-pass sort of clipping using
 * polysubs, and then finish the job by calling polyfix and
 * polystart. >*/
{
struct vertex  *v;
int             firstpoint, i;

    Allgone = 1;		/* Assume none left unless polyfix tells us */
    if (!dev.smart_clip)
    {
	firstpoint = 2;
	xminclip (0, 0, &firstpoint);	/* Tell them all to get ready */
    }

    firstpoint = 1;

    v = head;
    for (i = 0; i < npts; i++)
    {
	if (!dev.smart_clip)
	{
	    xminclip (v->x, v->y, &firstpoint);
	}
	else
	{
	    polyfix (v->x, v->y, &firstpoint);
	}
	v++;
    }
    if (!dev.smart_clip)
    {
	firstpoint = -1;	/* Means this was the last point! */
	xminclip (0, 0, &firstpoint);
    }
    if (Allgone == 0)		/* If still 1, means there's nothing left! */
    {
	polystart ();
    }
}
