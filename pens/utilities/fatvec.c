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
 *  source file:   ./filters/utilities/fatvec.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include <stdio.h>
#include <math.h>
#include "../include/extern.h"

void fatvec (int x1, int y1, int x2, int y2, int nfat, int dashon)
/*< Utility routine to make fat vectors from several thin ones.
 * Should ONLY be called if nfat > 0 and dashon = 0
 *
 * Algorithm by Glenn Kroeger
 * Changes added by Joe Dellinger to make it more efficient when plotting >*/
{
register int    i;
register int    fplus, fminus;
static int      lastdir = 0;

    lastdir = 1 - lastdir;

    if (dev.aspect_ratio != 1. && (y2 != y1 || x2 != x1) && nfat)
    {
	nfat = (float) nfat *sqrt (
	    ((y2 - y1) * (y2 - y1) + ((x2 - x1) * (x2 - x1))
	     /               (dev.aspect_ratio * dev.aspect_ratio)) /
	    ((y2 - y1) * (y2 - y1) + ((x2 - x1) * (x2 - x1)))
	    );
    }

    fminus = (nfat / 2);
    fplus = (nfat + 1) / 2;

    if (x1 <= x2)
    {
	if (y2 > y1)
	{
	    if (x1 == x2)
	    {
		for (i = -fminus; i <= fplus; i++)
		    dev.vector (x1 + i, y1 - fminus, x2 + i, y2 + fplus, 0, 0);
	    }
	    else
	    {

		if (lastdir)
		{
		    for (i = (fminus + fplus); i > 0; i--)
		    {
			dev.vector (x1 - fminus + i, y1 - fminus, x2 + fplus,
				    y2 + fplus - i, 0, 0);
		    }
		    dev.vector (x1 - fminus, y1 - fminus, x2 + fplus, y2 + fplus, 0, 0);
		    for (i = 1; i < (fminus + fplus + 1); i++)
		    {
			dev.vector (x1 - fminus, y1 - fminus + i, x2 + fplus - i,
				    y2 + fplus, 0, 0);
		    }
		}
		else
		{
		    for (i = (fminus + fplus); i > 0; i--)
		    {
			dev.vector (x1 - fminus, y1 - fminus + i, x2 + fplus - i,
				    y2 + fplus, 0, 0);
		    }
		    dev.vector (x1 - fminus, y1 - fminus, x2 + fplus, y2 + fplus, 0, 0);
		    for (i = 1; i < (fminus + fplus + 1); i++)
		    {
			dev.vector (x1 - fminus + i, y1 - fminus, x2 + fplus,
				    y2 + fplus - i, 0, 0);
		    }
		}
	    }
	}
	else
	if (y2 == y1)
	{
	    for (i = -fminus; i <= fplus; i++)
		dev.vector (x1 - fminus, y1 + i, x2 + fplus, y2 + i, 0, 0);
	}
	else
	{
	    if (x1 == x2)
	    {
		for (i = -fminus; i <= fplus; i++)
		    dev.vector (x1 + i, y1 + fplus, x2 + i, y2 - fminus, 0, 0);
	    }
	    else
	    {

		if (lastdir)
		{
		    for (i = (fminus + fplus); i > 0; i--)
		    {
			dev.vector (x1 - fminus + i, y1 + fplus, x2 + fplus,
				    y2 - fminus + i, 0, 0);
		    }
		    dev.vector (x1 - fminus, y1 + fplus, x2 + fplus, y2 - fminus, 0, 0);
		    for (i = 1; i < (fminus + fplus + 1); i++)
		    {
			dev.vector (x1 - fminus, y1 + fplus - i, x2 + fplus - i,
				    y2 - fminus, 0, 0);
		    }
		}
		else
		{
		    for (i = (fminus + fplus); i > 0; i--)
		    {
			dev.vector (x1 - fminus, y1 + fplus - i, x2 + fplus - i,
				    y2 - fminus, 0, 0);
		    }
		    dev.vector (x1 - fminus, y1 + fplus, x2 + fplus, y2 - fminus, 0, 0);
		    for (i = 1; i < (fminus + fplus + 1); i++)
		    {
			dev.vector (x1 - fminus + i, y1 + fplus, x2 + fplus,
				    y2 - fminus + i, 0, 0);
		    }
		}
	    }
	}
    }
    else
    {
	if (y2 > y1)
	{
	    if (lastdir)
	    {
		for (i = (fminus + fplus); i > 0; i--)
		{
		    dev.vector (x1 + fplus, y1 - fminus + i, x2 - fminus + i, y2 + fplus, 0, 0);
		}
		dev.vector (x1 + fplus, y1 - fminus, x2 - fminus, y2 + fplus, 0, 0);
		for (i = 1; i < (fminus + fplus + 1); i++)
		{
		    dev.vector (x1 + fplus - i, y1 - fminus, x2 - fminus, y2 + fplus - i, 0, 0);
		}
	    }
	    else
	    {
		for (i = (fminus + fplus); i > 0; i--)
		{
		    dev.vector (x1 + fplus - i, y1 - fminus, x2 - fminus, y2 + fplus - i, 0, 0);
		}
		dev.vector (x1 + fplus, y1 - fminus, x2 - fminus, y2 + fplus, 0, 0);
		for (i = 1; i < (fminus + fplus + 1); i++)
		{
		    dev.vector (x1 + fplus, y1 - fminus + i, x2 - fminus + i, y2 + fplus, 0, 0);
		}
	    }
	}
	else
	if (y2 == y1)
	{
	    for (i = -fminus; i <= fplus; i++)
		dev.vector (x1 + fplus, y1 + i, x2 - fminus, y2 + i, 0, 0);
	}
	else
	{
	    if (lastdir)
	    {
		for (i = (fminus + fplus); i > 0; i--)
		{
		    dev.vector (x1 + fplus, y1 + fplus - i, x2 - fminus + i, y2 - fminus, 0, 0);
		}
		dev.vector (x1 + fplus, y1 + fplus, x2 - fminus, y2 - fminus, 0, 0);
		for (i = 1; i < (fminus + fplus + 1); i++)
		{
		    dev.vector (x1 + fplus - i, y1 + fplus, x2 - fminus, y2 - fminus + i, 0, 0);
		}
	    }
	    else
	    {
		for (i = (fminus + fplus); i > 0; i--)
		{
		    dev.vector (x1 + fplus - i, y1 + fplus, x2 - fminus, y2 - fminus + i, 0, 0);
		}
		dev.vector (x1 + fplus, y1 + fplus, x2 - fminus, y2 - fminus, 0, 0);
		for (i = 1; i < (fminus + fplus + 1); i++)
		{
		    dev.vector (x1 + fplus, y1 + fplus - i, x2 - fminus + i, y2 - fminus, 0, 0);
		}
	    }
	}
    }
}
