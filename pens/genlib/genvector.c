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
 *  source file:   ./filters/genlib/genvector.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger (SOEST), June 23 1992
 *	Check the status of "lost" after every call to dev.plot,
 *	it might happen to change right after the first call of a
 *	vector pair.
 */

/*
 * Generic vector routine for devices that don't support fatness.
 * This version tries to be smart and minimize the "motion of the pen",
 * and tries to prolong strings of "draws" where possible.
 */
#include <stdio.h>
#include "../include/extern.h"

#include "../utilities/util.h"

#define MOVE 0
#define DRAW 1

#define ABS(a)   ((a) >= 0  ? (a) : (-(a)))

void genvector (int x1, int y1, int x2, int y2, int nfat, int dashon)
/*< device-independent vector >*/
{
static int      xlst, ylst;
int             d1, d2;

    if (nfat < 0)
	return;

    if (dashon)
    {
	dashvec (x1, y1, x2, y2, nfat, dashon);
	return;
    }

    if (nfat)			/* Recursively calls itself to make fat lines */
    {
	fatvec (x1, y1, x2, y2, nfat, dashon);
	return;
    }

    /*
     * Do clipping
     */
    if (!dev.smart_clip && clip (&x1, &y1, &x2, &y2))
	    return;
    /*
     * Important special case: Zero-length vector at the end of what you've
     * already plotted. Don't need to do anything.
     */
    if (x1 == x2 && y1 == y2 && !dev.lost && x1 == xlst && y1 == ylst)
    {
	return;
    }

    /*
     * Minimize movement of "pen"
     */
    if (!dev.lost)
    {
	d1 = ABS (x1 - xlst) + ABS (y1 - ylst);
	d2 = ABS (x2 - xlst) + ABS (y2 - ylst);
	if (d2 < d1)
	{
	    d1 = x1;
	    d2 = y1;
	    x1 = x2;
	    y1 = y2;
	    x2 = d1;
	    y2 = d2;
	}
    }

    if ((x1 != xlst) || (y1 != ylst) || dev.lost)
    {
	/* Make sure it is a move, not a draw */
	if (!dev.lost && ABS (x1 - xlst) <= 1 && ABS (y1 - ylst) <= 1)
	{
	    /*
	     * We're within one pixel, so go ahead and draw a vector to the
	     * new point. This avoids having to leave and re-enter vector
	     * mode.
	     */
	    dev.plot (x1, y1, DRAW);

	    /*
	     * It is remotely possible that JUST NOW, when we did that DRAW,
	     * we caused the device to become "lost"... So check for that
	     * unfortunate unlikely case, and fix it. (I've since learned why
	     * a plethora of externals isn't good programming style...)
	     */
	    if (dev.lost)
		dev.plot (x1, y1, MOVE);
	    /*
	     * If the device is STILL lost it's a BUG!
	     */
	}
	else
	{
	    dev.plot (x1, y1, MOVE);
	}
    }

    dev.plot (x2, y2, DRAW);
    xlst = x2;
    ylst = y2;
}
