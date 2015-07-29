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
 *  source file:   ./filters/utilities/dashvec.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger Mar 2 1987
 *	Simplify a few expressions that made the IBM RT's barf.
 */

#include <stdio.h>
#include <math.h>

#include <rsf.h>

#include "../include/extern.h"
#include "../include/round.h"

static double dashmod (double a, double b);

/*
 * dash1, gap1, dash2, gap2
 * in inches.
 */

void dashvec (int x1, int y1, int x2, int y2, int nfat, int dashon)
/*< Utility routine to make dashed lines. 
 * Should ONLY be called if dashon > 0 >*/
{
    double          dash_ahead, dash_behind, dashdist, dist, sine, cosine;
    double          deltax, deltay, lambda1, lambda2;
    int             i, nextdash;
    int             xv1, xv2, yv1, yv2;

/*
 * If not a dashed line, should never have even been called.
 * Can't dash a single point!
 */
    if (!dashon || (x1 == x2 && y1 == y2))
	return;

/*
 * find current position in dashes
 */
    dashpos = (float) dashmod ((double) dashpos, (double) dashsum);
    dash_behind = dashpos;

    i = 0;
    dash_ahead = 0.;
    while (dash_ahead < dash_behind)
    {
	dash_ahead += dashes[i];
	i++;
    }
    nextdash = i - 1;
/*
 * compute distances, properly scaled
 */
    deltax = x2 - x1;
    deltay = y2 - y1;
    lambda1 = 1. / dev.pixels_per_inch;
    lambda2 = dev.aspect_ratio / dev.pixels_per_inch;
    dist = sqrt (lambda1 * lambda1 * deltax * deltax + lambda2 * lambda2 * deltay * deltay);
    sine = deltay / dist;
    cosine = deltax / dist;

/*
 * draw the dashed line
 */
    for (dashdist = 0.;
	 dash_ahead - dashpos < dist;
	 dashdist += dash_ahead - dash_behind,
	     dash_behind = dash_ahead,
	     nextdash++,
	     dash_ahead += dashes[nextdash % (dashon * 2)])
    {
	if (nextdash % 2 == 0)
	{
	    xv1 = ROUND (x1 + cosine * (dashdist));
	    yv1 = ROUND (y1 + sine * (dashdist));
	    xv2 = ROUND (x1 + cosine * (dashdist + (dash_ahead - dash_behind)));
	    yv2 = ROUND (y1 + sine * (dashdist + (dash_ahead - dash_behind)));

	    dev.vector (xv1, yv1, xv2, yv2, nfat, 0);
	}
    }

    if (nextdash % 2 == 0)
    {
	xv1 = ROUND (x1 + cosine * (dashdist));
	yv1 = ROUND (y1 + sine * (dashdist));

	dev.vector (xv1, yv1, x2, y2, nfat, 0);
    }

/*
 * Increment to new position in dash pattern.
 */
    dashpos += dashmod (dist, (double) dashsum);
}

/*
 * mod subroutine
 */
static double dashmod (double a, double b)
{
    return (a - floor (a / b) * b);
}
