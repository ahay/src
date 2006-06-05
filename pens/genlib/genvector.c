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

extern int      smart_clip;
extern int      lost;

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
    if (!smart_clip)
	if (clip (&x1, &y1, &x2, &y2))
	    return;
    /*
     * Important special case: Zero-length vector at the end of what you've
     * already plotted. Don't need to do anything.
     */
    if (x1 == x2 && y1 == y2 && !lost && x1 == xlst && y1 == ylst)
    {
	return;
    }

    /*
     * Minimize movement of "pen"
     */
    if (!lost)
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

    if ((x1 != xlst) || (y1 != ylst) || lost)
    {
	/* Make sure it is a move, not a draw */
	if (!lost && ABS (x1 - xlst) <= 1 && ABS (y1 - ylst) <= 1)
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
	    if (lost)
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
