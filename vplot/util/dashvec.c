/*
 * Copyright 1987 the Board of Trustees of the Leland Stanford Junior
 * University. Official permission to use this software is included in
 * the documentation. It authorizes you to use this file for any
 * non-commercial purpose, provided that this copyright notice is not
 * removed and that any modifications made to this file are commented
 * and dated in the style of my example below.
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

/*
 * Utility routine to make dashed lines.
 * Should ONLY be called if dashon > 0
 */

#include <stdio.h>
#include <math.h>

#include "_round.h"
#include "device.h"

static double mod (double a, double b);

void dashvec (device dev,
	      int x1, int y1, int x2, int y2, int nfat, int dashon)
{
    double  dash_ahead, dash_behind, dashdist, dist, sine, cosine;
    double  deltax, deltay, lambda1, lambda2;
    int i, nextdash;
    int xv1, xv2, yv1, yv2;

/*
 * If not a dashed line, should never have even been called.
 * Can't dash a single point!
 */
    if (!dashon || (x1 == x2 && y1 == y2))
	return;

/*
 * find current position in dashes
 */
    dev->dashpos = (float) mod ((double) dev->dashpos, (double) dev->dashsum);
    dash_behind = dev->dashpos;

    i = 0;
    dash_ahead = 0.;
    while (dash_ahead < dash_behind)
    {
	dash_ahead += dev->dashes[i];
	i++;
    }
    nextdash = i - 1;
/*
 * compute distances, properly scaled
 */
    deltax = x2 - x1;
    deltay = y2 - y1;
    lambda1 = 1. / dev->pixels_per_inch;
    lambda2 = dev->aspect_ratio / dev->pixels_per_inch;
    dist = hypot (lambda1 * deltax, lambda2 * deltay);
    sine = deltay / dist;
    cosine = deltax / dist;

/*
 * draw the dashed line
 */
    for (dashdist = 0.;
	 dash_ahead - dev->dashpos < dist;
	 dashdist += dash_ahead - dash_behind,
	     dash_behind = dash_ahead,
	     nextdash++,
	     dash_ahead += dev->dashes[nextdash % (dashon * 2)])
    {
	if (nextdash % 2 == 0)
	{
	    xv1 = ROUND (x1 + cosine * (dashdist));
	    yv1 = ROUND (y1 + sine * (dashdist));
	    xv2 = ROUND (x1 + cosine * (dashdist + (dash_ahead - dash_behind)));
	    yv2 = ROUND (y1 + sine * (dashdist + (dash_ahead - dash_behind)));
	    
	    dev->vector (xv1, yv1, xv2, yv2, nfat, 0);
	}
    }
    
    if (nextdash % 2 == 0)
    {
	xv1 = ROUND (x1 + cosine * (dashdist));
	yv1 = ROUND (y1 + sine * (dashdist));
	
	dev->vector (xv1, yv1, x2, y2, nfat, 0);
    }

    dev->dashpos += mod (dist, (double) dev->dashsum);
}

static double mod (double a, double b)
{
    return (a - floor (a / b) * b);
}
