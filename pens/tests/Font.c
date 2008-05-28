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
 *  source file:   ./filters/Tests/Font.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include <stdio.h>
#include <stdlib.h>

#include <rsfplot.h>

#define NUMFONTS 17
char           *fontnames[] =
{
 "pen",
 "romans",
 "romand",
 "romanc",
 "romant",
 "italicc",
 "italict",
 "scripts",
 "scriptc",
 "greeks",
 "greekc",
 "cyrilc",
 "gothgbt",
 "gothgrt",
 "gothitt",
 "math",
 "misc",
 ""
};

int main (int argc, char *argv[])
{
    int             ii, xx, yy, font;
    
    if (argc == 1)
    {
	fprintf (stderr, "usage: Font font_number | plas | ?pen\n");
	exit (1);
    }
    font = atoi (argv[1]);

    printf ("S s\n");

    printf ("m 2400 30\n");
    printf ("F 0 1 0\n");
    if (font < NUMFONTS)
	printf ("T 5 0\n\\f2 Font %d: %s\n", font, fontnames[font]);
    else
	printf ("T 5 0\n\\f2 Font %d: hardware\n", font);

    for (ii = 30; ii < 200; ii++)
    {
	xx = (ii - 30) % 10;
	yy = (ii - 30) / 10;

	xx = xx * 550 + 300;
	yy = yy * 255 + 300;

	printf ("m %d %d\n", xx, yy);
	printf ("F 0 1 0\nJ %d %d\n", TH_RIGHT, TV_NORMAL);
	printf ("T 3 0\n%d:\\v%d :\n", ii, ii);
	printf ("m %d %d\n", xx + 10, yy);
	printf ("F %d 1 0\nJ %d %d\n", font, TH_NORMAL, TV_NORMAL);
	printf ("T 8 0\n\\v%d \n", ii);
    }

    exit(0);
}
