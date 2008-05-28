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
 *  source file:   ./filters/genlib/polysubs.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include <stdio.h>
#include "../include/extern.h"
#define OUT 0
#define IN  1
#define UNSET -1

#include "polyfix.h"

static void yminclip (int xin, int yin, int *first);
static void xmaxclip (int xin, int yin, int *first);
static void ymaxclip (int xin, int yin, int *first);

int inter (int x1, int x2, int y1, int y2, int x)
{
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1);
}

void xminclip (int xin, int yin, int *first)
/*< Do a simple-minded polygon clipping. If it goes out, draw it to where it
 * crossed the edge. When it comes back in, draw it from where it hit the
 * edge. This is complicated if you have to deal with several edges, but
 * very easy for one edge. So, since I'm lazy, I did it four times, each
 * routine reading in points, clipping, and sending the output on for
 * further clipping. Finally, it sends it on to polyfix which does
 * complicated processing. JAD 5-11-84 >*/
{
static int      xstart, ystart;
static int      ostatus;
int             status;
static int      firstout;
static int      xold, yold;

    if (*first == 2)
    {
	ostatus = UNSET;
	firstout = 2;
	yminclip (0, 0, &firstout);
	return;
    }

    if (*first == -1)
    {
	if (ostatus == UNSET)
	{
	    /* We never got anything! */
	    return;
	}
	/* finish up */
	xin = xstart;
	yin = ystart;
    }

    status = IN;
    if (xin < xwmin)
	status = OUT;

    if (*first == 1)
    {
	/* This is the first time we have been called */
	xstart = xin;
	ystart = yin;
	firstout = 1;
	*first = 0;
	ostatus = status;
	xold = xin;
	yold = yin;
	return;
    }
/* Not our first time */

    switch (status)
    {
    case IN:
	switch (ostatus)
	{
	case IN:
	    /* in this time, in last time */
	    yminclip (xin, yin, &firstout);
	    break;
	case OUT:
	    /* out last time, in now */
	    /* find where we came in! */
	    yminclip (xwmin, inter (xold, xin, yold, yin, xwmin), &firstout);
	    yminclip (xin, yin, &firstout);
	    break;
	}
	break;
    case OUT:
	switch (ostatus)
	{
	case IN:
	    /* in last time, out now */
	    /* find where we went out */
	    yminclip (xwmin, inter (xold, xin, yold, yin, xwmin), &firstout);
	    break;
	case OUT:
	    /* out last time, still out */
	    /* don't output anything */
	    break;
	}
	break;
    }
    if (*first == -1)
    {
	firstout = -1;
	yminclip (0, 0, &firstout);
    }
    else
    {
	xold = xin;
	yold = yin;
	ostatus = status;
    }
}


static void yminclip (int xin, int yin, int *first)
{
static int      xstart, ystart;
static int      ostatus;
int             status;
static int      firstout;
static int      xold, yold;

    if (*first == 2)
    {
	ostatus = UNSET;
	firstout = 2;
	xmaxclip (0, 0, &firstout);
	return;
    }

    if (*first == -1)
    {
	if (ostatus == UNSET)
	{
	    /* We never got anything! */
	    return;
	}
	/* finish up */
	xin = xstart;
	yin = ystart;
    }

    status = IN;
    if (yin < ywmin)
	status = OUT;

    if (*first == 1)
    {
	/* This is the first time we have been called */
	xstart = xin;
	ystart = yin;
	firstout = 1;
	*first = 0;
	ostatus = status;
	xold = xin;
	yold = yin;
	return;
    }
/* Not our first time */

    switch (status)
    {
    case IN:
	switch (ostatus)
	{
	case IN:
	    /* in this time, in last time */
	    xmaxclip (xin, yin, &firstout);
	    break;
	case OUT:
	    /* out last time, in now */
	    /* find where we came in! */
	    xmaxclip (inter (yold, yin, xold, xin, ywmin), ywmin, &firstout);
	    xmaxclip (xin, yin, &firstout);
	    break;
	}
	break;
    case OUT:
	switch (ostatus)
	{
	case IN:
	    /* in last time, out now */
	    /* find where we went out */
	    xmaxclip (inter (yold, yin, xold, xin, ywmin), ywmin, &firstout);
	    break;
	case OUT:
	    /* out last time, still out */
	    /* don't output anything */
	    break;
	}
	break;
    }
    if (*first == -1)
    {
	firstout = -1;
	xmaxclip (0, 0, &firstout);
    }
    else
    {
	xold = xin;
	yold = yin;
	ostatus = status;
    }
}

static void xmaxclip (int xin, int yin, int *first)
{
static int      xstart, ystart;
static int      ostatus;
int             status;
static int      firstout;
static int      xold, yold;

    if (*first == 2)
    {
	ostatus = UNSET;
	firstout = 2;
	ymaxclip (0, 0, &firstout);
	return;
    }

    if (*first == -1)
    {
	if (ostatus == UNSET)
	{
	    /* We never got anything! */
	    return;
	}
	/* finish up */
	xin = xstart;
	yin = ystart;
    }

    status = IN;
    if (xin > xwmax)
	status = OUT;

    if (*first == 1)
    {
	/* This is the first time we have been called */
	xstart = xin;
	ystart = yin;
	firstout = 1;
	*first = 0;
	ostatus = status;
	xold = xin;
	yold = yin;
	return;
    }
/* Not our first time */

    switch (status)
    {
    case IN:
	switch (ostatus)
	{
	case IN:
	    /* in this time, in last time */
	    ymaxclip (xin, yin, &firstout);
	    break;
	case OUT:
	    /* out last time, in now */
	    /* find where we came in! */
	    ymaxclip (xwmax, inter (xold, xin, yold, yin, xwmax), &firstout);
	    ymaxclip (xin, yin, &firstout);
	    break;
	}
	break;
    case OUT:
	switch (ostatus)
	{
	case IN:
	    /* in last time, out now */
	    /* find where we went out */
	    ymaxclip (xwmax, inter (xold, xin, yold, yin, xwmax), &firstout);
	    break;
	case OUT:
	    /* out last time, still out */
	    /* don't output anything */
	    break;
	}
	break;
    }
    if (*first == -1)
    {
	firstout = -1;
	ymaxclip (0, 0, &firstout);
    }
    else
    {
	xold = xin;
	yold = yin;
	ostatus = status;
    }
}


static void ymaxclip (int xin, int yin, int *first)
{
static int      xstart, ystart;
static int      ostatus;
int             status;
static int      firstout;
static int      xold, yold;

    if (*first == 2)
    {
	ostatus = UNSET;
	return;
    }

    if (*first == -1)
    {
	if (ostatus == UNSET)
	{
	    /* We never got anything! */
	    return;
	}
	/* finish up */
	xin = xstart;
	yin = ystart;
    }

    status = IN;
    if (yin > ywmax)
	status = OUT;

    if (*first == 1)
    {
	/* This is the first time we have been called */
	xstart = xin;
	ystart = yin;
	firstout = 1;
	*first = 0;
	ostatus = status;
	xold = xin;
	yold = yin;
	return;
    }
/* Not our first time */

    switch (status)
    {
    case IN:
	switch (ostatus)
	{
	case IN:
	    /* in this time, in last time */
	    polyfix (xin, yin, &firstout);
	    break;
	case OUT:
	    /* out last time, in now */
	    /* find where we came in! */
	    polyfix (inter (yold, yin, xold, xin, ywmax), ywmax, &firstout);
	    polyfix (xin, yin, &firstout);
	    break;
	}
	break;
    case OUT:
	switch (ostatus)
	{
	case IN:
	    /* in last time, out now */
	    /* find where we went out */
	    polyfix (inter (yold, yin, xold, xin, ywmax), ywmax, &firstout);
	    break;
	case OUT:
	    /* out last time, still out */
	    /* don't output anything */
	    break;
	}
	break;
    }
    if (*first == -1)
    {
	/* We're done! */
	return;
    }
    else
    {
	xold = xin;
	yold = yin;
	ostatus = status;
    }
}
