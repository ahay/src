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

#include <stdio.h>
#include <math.h>
#include "../include/extern.h"
#include "../include/round.h"

extern int      no_stretch_text;

void vptodevxy (int x, int y, int *outx, int *outy)
/*< convert vplot coordinates to device coordinates >*/
{
    float tempx, tempy, temp;

    tempx = (float) (x - dev.xorigin) * xscale;
    tempy = (float) (y - dev.yorigin) * yscale;

    temp = mxx * tempx + mxy * tempy;
    tempy = myx * tempx + myy * tempy;
    tempx = temp;

    tempx = tempx * hdevscale + dev.xmin + hshift;
    tempy = tempy * vdevscale + dev.ymin + vshift;

    *outx = ROUND (tempx);
    *outy = ROUND (tempy);
}

void vptodevw (int x1, int y1, int x2, int y2, 
	       int *x1out, int *y1out, int *x2out, int *y2out)
/*< convert vplot coordinates to device coordinates >*/
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

void devtovpxy (int x, int y, int *outx, int *outy)
/*< convert device coordinates to vplot coordinates >*/
{
float           tempx, tempy, temp;

    tempx = (float) (x - dev.xmin - hshift) / hdevscale;
    tempy = (float) (y - dev.ymin - vshift) / vdevscale;

    temp = mxx * tempx - mxy * tempy;
    tempy = -myx * tempx + myy * tempy;
    tempx = temp;

    tempx = tempx / xscale + dev.xorigin;
    tempy = tempy / yscale + dev.yorigin;

    *outx = ROUND (tempx);
    *outy = ROUND (tempy);
}

void devtovpw (int x1, int y1, int x2, int y2, 
	       int *x1out, int *y1out, int *x2out, int *y2out)
/*< convert device coordinates to vplot coordinates >*/
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

void vptodevxy_text (int x, int y, int *outx, int *outy)
/*< convert vplot coordinates to device coordinates >*/
{
    float           xscale_save=0., yscale_save=0.;

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
