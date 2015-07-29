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
