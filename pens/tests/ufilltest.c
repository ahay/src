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
