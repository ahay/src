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
 *  source file:   ./filters/Tests/texttest.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include <stdio.h>
#include <math.h>

#include <rsfplot.h>

#define NP 7

int main (void)
{
    float           xp, yp;
    int             i;

    vp_init();
/*
 * Decide what "style" plot this is to be.
 */
    vp_style (VP_STANDARD);

/*
 * text
 */
    for (i = 0; i < NP; i++)
    {
	xp = 1. + (7. / NP) * i;
	yp = xp;
	vp_color (i + 1);
	vp_tfont (i * 2, VP_STRING, OVLY_NORMAL);
	vp_text (xp, yp, 15, 0, "This is a test...");
/*	vp_utext(xp,yp,15,0,"This is a test..."); */
    }

/*
 * Finish up
 */
    return 0;
}
