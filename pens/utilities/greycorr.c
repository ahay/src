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
 *  source file:   ./filters/utilities/greycorr.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger March 28 1988
 *	Do invras here.
 */

#include <stdio.h>
#include "../include/extern.h"
#include "../include/round.h"

int greycorr (int colornum)
/*< Utility to modify color tables for plotting grey rasters. >*/
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
