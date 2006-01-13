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
 *  source file:   ./filters/utilities/greycorr.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger March 28 1988
 *	Do invras here.
 */

/*
 * Utility to modify color tables for plotting grey rasters.
 */

#include <stdio.h>
#include "../include/extern.h"
#include "../include/round.h"

int
greycorr (colornum)
    int             colornum;
{
float           newval;
extern float    greyc, pixc;
extern bool      invras;
int             outval;

    newval = colornum;

    if (invras)
	newval = 255 - newval;

/* 
 * correction to simulate nonlinearity of graphics displays
 */
    if (greyc != 1.)
    {
	newval /= 255.;
	newval = (-2. + 2. * greyc) * newval * newval * newval + 3. * (1. - greyc) * newval * newval + greyc * newval;
	newval *= 255.;
	if (newval < 0)
	    newval = 0.;
	if (newval > 255.)
	    newval = 255.;
    }

/*
 * correction for pixel overlap on hardcopy devices
 */
    if (pixc != 1.)
    {
	if (newval < pixc * 128.)
	{
	    newval /= pixc;
	}
	else
	{
	    newval = 128. + (newval - pixc * 128.) / (2. - pixc);
	}
	if (newval < 0)
	    newval = 0.;
	if (newval > 255.)
	    newval = 255.;
    }

    outval = ROUND (newval);

    return outval;
}
