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
 *  source file:   ./filters/utilities/vptodev.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger (SEP), Oct 28 1988
 *	Put in "txsquare=no" option.
 * Joe Dellinger (SEP), May 11 1989
 *	Better to keep text within bounding box instead of
 *	preserving area.
 */

/*
 * convert vplot coordinates to device coordinates, or vice versa
 */
#include <stdio.h>
#include <math.h>
#include "../include/extern.h"
#include "../include/round.h"

extern int      no_stretch_text;

vptodevxy (x, y, outx, outy)
    int             x, y;
    int            *outx, *outy;
{
float           tempx, tempy, temp;


    tempx = (float) (x - xorigin) * xscale;
    tempy = (float) (y - yorigin) * yscale;

    temp = mxx * tempx + mxy * tempy;
    tempy = myx * tempx + myy * tempy;
    tempx = temp;

    tempx = tempx * hdevscale + dev_xmin + hshift;
    tempy = tempy * vdevscale + dev_ymin + vshift;

    *outx = ROUND (tempx);
    *outy = ROUND (tempy);
}

vptodevw (x1, y1, x2, y2, x1out, y1out, x2out, y2out)
    int             x1, y1, x2, y2;
    int            *x1out, *y1out, *x2out, *y2out;
{
int             x11, y11, x12, y12, x21, y21, x22, y22, a, b;

    vptodevxy (x1, y1, &x11, &y11);
    vptodevxy (x1, y2, &x12, &y12);
    vptodevxy (x2, y1, &x21, &y21);
    vptodevxy (x2, y2, &x22, &y22);

    a = (x11 > x12 ? x11 : x12);
    b = (x22 > x21 ? x22 : x21);
    *x2out = (a > b ? a : b);

    a = (y11 > y12 ? y11 : y12);
    b = (y22 > y21 ? y22 : y21);
    *y2out = (a > b ? a : b);

    a = (x11 < x12 ? x11 : x12);
    b = (x22 < x21 ? x22 : x21);
    *x1out = (a < b ? a : b);

    a = (y11 < y12 ? y11 : y12);
    b = (y22 < y21 ? y22 : y21);
    *y1out = (a < b ? a : b);
}

devtovpxy (x, y, outx, outy)
    int             x, y;
    int            *outx, *outy;
{
float           tempx, tempy, temp;

    tempx = (float) (x - dev_xmin - hshift) / hdevscale;
    tempy = (float) (y - dev_ymin - vshift) / vdevscale;

    temp = mxx * tempx - mxy * tempy;
    tempy = -myx * tempx + myy * tempy;
    tempx = temp;

    tempx = tempx / xscale + xorigin;
    tempy = tempy / yscale + yorigin;

    *outx = ROUND (tempx);
    *outy = ROUND (tempy);
}

devtovpw (x1, y1, x2, y2, x1out, y1out, x2out, y2out)
    int             x1, y1, x2, y2;
    int            *x1out, *y1out, *x2out, *y2out;
{
int             x11, y11, x12, y12, x21, y21, x22, y22, a, b;

    devtovpxy (x1, y1, &x11, &y11);
    devtovpxy (x1, y2, &x12, &y12);
    devtovpxy (x2, y1, &x21, &y21);
    devtovpxy (x2, y2, &x22, &y22);

    a = (x11 > x12 ? x11 : x12);
    b = (x22 > x21 ? x22 : x21);
    *x2out = (a > b ? a : b);

    a = (y11 > y12 ? y11 : y12);
    b = (y22 > y21 ? y22 : y21);
    *y2out = (a > b ? a : b);

    a = (x11 < x12 ? x11 : x12);
    b = (x22 < x21 ? x22 : x21);
    *x1out = (a < b ? a : b);

    a = (y11 < y12 ? y11 : y12);
    b = (y22 < y21 ? y22 : y21);
    *y1out = (a < b ? a : b);
}

vptodevxy_text (x, y, outx, outy)
    int             x, y;
    int            *outx, *outy;
{
float           xscale_save, yscale_save;

    if (no_stretch_text)
    {
	xscale_save = xscale;
	yscale_save = yscale;

	if (fabs (xscale) < fabs (yscale))
	{
	    yscale = xscale;
	}
	else
	{
	    xscale = yscale;
	}
/*
 * This makes more sense geometrically, but it seems better to
 * keep text within the stretched bounding box.
 *
 *	xscale = sqrt (xscale_save * yscale_save);
 *	yscale = xscale;
 */
    }

    vptodevxy (x, y, outx, outy);

    if (no_stretch_text)
    {
	xscale = xscale_save;
	yscale = yscale_save;
    }
}
