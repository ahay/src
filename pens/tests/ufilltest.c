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
 *  source file:   ./filters/Tests/ufilltest.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include <stdio.h>
#include <math.h>

#include <rsfplot.h>

#define NP 4

int main (void)
{
    float           xarray[NP], yarray[NP];
    int             angle, numhatch;
    int             hatcharray[8];

    vp_init();

/*
 * Decide what "style" plot this is to be.
 */
    vp_style (VP_STANDARD);

/* 
 * solid fill
 */
    xarray[0] = 2.;
    yarray[0] = 2.;
    xarray[1] = 2.;
    yarray[1] = 6.;
    xarray[2] = 6.;
    yarray[2] = 6.;
    xarray[3] = 6.;
    yarray[3] = 2.;

    vp_color (VP_RED);
    vp_ufill (xarray, yarray, NP);
/*
 * define a hatching pattern and fill.
 */
    angle = 30;
    numhatch = 1;
    hatcharray[0] = 1;
    hatcharray[1] = VP_CYAN;
    hatcharray[2] = 0;
    hatcharray[3] = 20;
    hatcharray[4] = 1;
    hatcharray[5] = VP_WHITE;
    hatcharray[6] = 0;
    hatcharray[7] = 10;

    xarray[0] = 2.;
    yarray[0] = 2.;
    xarray[1] = 2.;
    yarray[1] = 6.;
    xarray[2] = 6.;
    yarray[2] = 6.;
    xarray[3] = 6.;
    yarray[3] = 2.;

    vp_hatchload (angle, numhatch, VP_CYAN, hatcharray);
    vp_color (VP_CYAN);
    vp_ufill (xarray, yarray, NP);

    return 0;
}
