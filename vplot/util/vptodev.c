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
#include <math.h>

#include "_round.h"
#include "device.h"

void vptodevxy (device dev, int x, int y, int *outx, int *outy)
/*< convert vplot coordinates to device coordinates >*/
{
    float tempx, tempy, temp;
    
    tempx = (float) (x - dev->xorigin) * dev->xscale;
    tempy = (float) (y - dev->yorigin) * dev->yscale;
    
    temp  = dev->mxx * tempx + dev->mxy * tempy;
    tempy = dev->myx * tempx + dev->myy * tempy;
    tempx = temp;

    tempx = tempx * dev->hscale + dev->xmin + dev->hshift;
    tempy = tempy * dev->vscale + dev->ymin + dev->vshift;

    *outx = ROUND (tempx);
    *outy = ROUND (tempy);
}

void vptodevw (device dev, int x1, int y1, int x2, int y2, 
	       int *x1out, int *y1out, int *x2out, int *y2out)
/*< convert vplot coordinates to device coordinates >*/
{
    int x11, y11, x12, y12, x21, y21, x22, y22, a, b;

    vptodevxy (dev, x1, y1, &x11, &y11);
    vptodevxy (dev, x1, y2, &x12, &y12);
    vptodevxy (dev, x2, y1, &x21, &y21);
    vptodevxy (dev, x2, y2, &x22, &y22);

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

void devtovpxy (device dev, int x, int y, int *outx, int *outy)
/*< convert device coordinates to vplot coordinates >*/
{
    float tempx, tempy, temp;

    tempx = (float) (x - dev->xmin - dev->hshift) / dev->hscale;
    tempy = (float) (y - dev->ymin - dev->vshift) / dev->vscale;

    temp  =  dev->mxx * tempx - dev->mxy * tempy;
    tempy = -dev->myx * tempx + dev->myy * tempy;
    tempx = temp;

    tempx = tempx / dev->xscale + dev->xorigin;
    tempy = tempy / dev->yscale + dev->yorigin;

    *outx = ROUND (tempx);
    *outy = ROUND (tempy);
}

void devtovpw (device dev, int x1, int y1, int x2, int y2, 
	       int *x1out, int *y1out, int *x2out, int *y2out)
/*< convert device coordinates to vplot coordinates >*/
{
    int x11, y11, x12, y12, x21, y21, x22, y22, a, b;

    devtovpxy (dev, x1, y1, &x11, &y11);
    devtovpxy (dev, x1, y2, &x12, &y12);
    devtovpxy (dev, x2, y1, &x21, &y21);
    devtovpxy (dev, x2, y2, &x22, &y22);

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

void vptodevxy_text (device dev, int x, int y, int *outx, int *outy)
/*< convert vplot coordinates to device coordinates >*/
{
    float xscale, yscale;

    if (dev->no_stretch_text) {
	xscale = dev->xscale;
	yscale = dev->yscale;

	if (fabsf (xscale) < fabsf (yscale)) {
	    yscale = xscale;
	} else {
	    xscale = yscale;
	}

/*
 * This makes more sense geometrically, but it seems better to
 * keep text within the stretched bounding box.
 *
 *	dev->xscale = sqrt (xscale * yscale);
 *	dev->yscale = xscale;
 */
    }

    vptodevxy (dev, x, y, outx, outy);

    if (dev->no_stretch_text) {
	dev->xscale = xscale;
	dev->yscale = yscale;
    }
}
