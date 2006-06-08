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
