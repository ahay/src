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
 *  source file:   ./filters/Tests/polytest.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include <stdio.h>
#include <math.h>

#include <rsfplot.h>

#define NP 10

int main (void)
{
    int             mtype, msize;
    float           xarray[NP], yarray[NP];
    float           dash[2], gap[2];
    int             i, j;
    int             pattern[100];

    vp_init();

/*
 * Decide what "style" plot this is to be.
 */
    vp_style (VP_STANDARD);

/*
 * Draw polymarkers.
 */
    for (i = 0; i < NP; i++)
    {
	xarray[i] = 4.25 + i % 3;
	yarray[i] = (1. - (float) i / NP) * 8.;
    }
    msize = 20;
    mtype = 23;
    vp_color (VP_RED);
    vp_pmark (NP, mtype, msize, xarray, yarray);

/*
 * Draw dashed polyline.
 */
    dash[0] = (.5);
    dash[1] = (.1);
    gap[0] = (.1);
    gap[1] = (.1);
    vp_setdash (dash, gap, 2);
    vp_color (VP_WHITE);

    for (i = 0; i < NP; i++)
    {
	xarray[i] = 4.25 + (i + 1) % 3;
	yarray[i] = (1. - (float) i / NP) * 8.;
    }
    vp_pline (xarray, yarray, NP);

    vp_color (VP_BLUE);

    for (i = 0; i < NP; i++)
    {
	xarray[i] = 4.25 + (i + 2) % 3;
	yarray[i] = (1. - (float) i / NP) * 8.;
    }
    vp_pline (xarray, yarray, NP);

    for (i = 0; i < 10; i++)
	for (j = 0; j < 10; j++)
	{
	    pattern[i + j * 10] = (int)
	     (.25 * ((i - 6.) * (i - 3.) + (j - 6.) * (j - 3.)));
	    if (pattern[i + j * 10] > 7 || pattern[i + j * 10] < 0)
		pattern[i + j * 10] = 7;
	}
    vp_patload (10, 10, 10, VP_GREEN, pattern);

    i = 0;
    pattern[i++] = 1;
    pattern[i++] = VP_RED;
    pattern[i++] = 0. * HATCHPERIN;
    pattern[i++] = (.5) * HATCHPERIN;

    pattern[i++] = 1;
    pattern[i++] = VP_WHITE;
    pattern[i++] = (.25) * HATCHPERIN;
    pattern[i++] = (.5) * HATCHPERIN;

    pattern[i++] = 1;
    pattern[i++] = VP_BLUE;
    pattern[i++] = 0. * HATCHPERIN;
    pattern[i++] = (.5) * HATCHPERIN;

    pattern[i++] = 1;
    pattern[i++] = VP_WHITE;
    pattern[i++] = (.25) * HATCHPERIN;
    pattern[i++] = (.5) * HATCHPERIN;

    vp_hatchload (30, 2, VP_RED, pattern);

    vp_color (VP_GREEN);
    for (i = 0; i < NP; i++)
    {
	xarray[i] = 9. + 2. * cos (i * 2. * 3.14159 / NP);
	yarray[i] = 4. + 2. * sin (i * 2. * 3.14159 / NP);
    }
    vp_fill (xarray, yarray, NP);

    vp_color (VP_RED);
    for (i = 0; i < NP; i++)
    {
	xarray[i] = 2. + 2. * cos (i * 2. * 3.14159 / NP);
	yarray[i] = 4. + 2. * sin (i * 2. * 3.14159 / NP);
    }
    vp_fill (xarray, yarray, NP);

    vp_fat (5);
    vp_setdash (NULL, NULL, 0);
    vp_pline (xarray, yarray, NP);

    return 0;
}
