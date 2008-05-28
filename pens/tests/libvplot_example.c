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
 *  source file:   ./filters/Tests/libvplot_example.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include    <math.h>
#include    <stdio.h>

#include    <rsfplot.h>

/*
 * Sample program to demonstrate simple libvplot usage.
 *      Excerpted from 'legend.c', by C. R. Karish:
 *      format entries for a legend for a geologic map.
 *  To compile:
 *      cc -o example vplot_example.c -lvplot -lm
 *  To run:
 *      example < infile > outfile
 *      tube outfile
 *  Sample input:       (remove asterisks)
 *      unit "Qa"  "Quaternary alluvium"
 *      unit "Tb"  "Tertiary basalt"
 *      unit "Kg"  "biotite monzogranite"
 *      unit "Krp"  "rhyolite porphyry"
 *      unit "Krd"  "biotite rhyodacite"
 */

/*
 * external parameters, used by many subroutines
 */
static float           boxsz[2];	/* size of the symbol box, in inches */
static float           tab[2];	/* locations of tab stops, in inches */
static int             charsz;	/* text size used for explanations (in vplot size
			 * units) */
static float           vspace;	/* vertical spacing for rock unit entries */
static float           vpos;	/* present location on the page (inches) */

static void mapunit (char *symbol, char *text);
static void feedpage (void);
static void drawbox (float xpos, float ypos);

#define MAXLINE 4096

int main (void)
{
char            command[30], line[125], symbol[30], text[90];
int             count;

    /*
     * set global geometric parameters 
     */

    vpos = 10.3;		/* start at 10.3 inches from bottom of page */
    vspace = 0.90;		/* base line spacing = 0.90 inches */
    charsz = 10;		/* default text size = .30 inches */
    boxsz[0] = 1.0;		/* default box: 1.0 x 0.70 inches */
    boxsz[1] = 0.70;
    tab[0] = 1.0;		/* default tabs: 1.0 and 3.5 inches */
    tab[1] = 3.5;

    /*
     * set up a vplot output stream and coordinate system 
     */
    vp_init();
    vp_style (VP_STANDARD);	/* origin at lower left; y axis vertical */
    vp_clip (0., 0., 8.5, 11.0);/* set clipping window */
    /* (wider than screen!) */
    vp_tfont (1, VP_STROKE, 0);	/* use Roman simplex, full precision, */
    /* no special overlay processing */
    vp_fat ((int) (charsz / 2));/* set line width, scaled to text size */
    vp_orig (0., 0.);		/* set the user origin to coincide with */
    /* the device origin (same as default) */

    /*
     * main loop: read and interpret input file 
     */
    while (fgets (line,MAXLINE,stdin) != NULL)
    {
	if (sscanf (line, "%s", command) == 0)	/* skip blank lines */
	    continue;
	if (command[0] == '#')	/* ... and comments */
	    continue;
	if ((count = sscanf (line,	/* parse strings from input */
			   "%*s \"%[^\"]\" \"%[^\"]\"", symbol, text)) != 2)
	{
	    fprintf (stderr, "unit count = %d\n", count);
	    fprintf (stderr, "%s\n", line);	/* complain if argument */
	    continue;		/* count is wrong */
	}
	mapunit (symbol, text);	/* produce some output */
    }				/* end of main loop */
    
    return 0;
}

static void mapunit (char *symbol, char *text)	
/* produce an entry for a map legend */
{
    if (vpos < 1.0)		/* leave 1" bottom margin */
	feedpage ();
    vpos -= vspace / 2;		/* space down a half line */
    if (strcmp (symbol, "nobox"))	/* strcmp==TRUE means `no match' */
    {
	drawbox (tab[0], vpos);	/* draw the box */
	vp_tjust (TH_CENTER, TV_HALF);	/* center the symbol in the box */
	vp_fat ((int) ((charsz - 1) / 2));	/* pick a good text fatness */
	vp_text (tab[0] + (boxsz[0] / 2),
		 vpos, charsz - 1, 0, symbol);	/* write the map symbol */
    }
    vp_fat ((int) (charsz / 2));
    vp_tjust (TH_LEFT, TV_HALF);/* reset left justification */
    vp_text (tab[1], vpos, charsz, 0, text);	/* write the explanation */
    vpos -= vspace / 2;		/* space down a half line */
}


/*
 * draw a box, boxsz[0] by boxsz[1] inches, left edge at xpos,
 *      centered vertically at ypos:
 */

static void drawbox (float xpos, float ypos)
{
    vp_fat ((int) (charsz / 2));/* set linewidth, same as for text */
    vp_move (xpos, ypos - (boxsz[1] / 2));	/* move to the upper left
						 * corner */
    vp_draw (xpos + boxsz[0], ypos - (boxsz[1] / 2));	/* draw it. */
    vp_draw (xpos + boxsz[0], ypos + (boxsz[1] / 2));
    vp_draw (xpos, ypos + (boxsz[1] / 2));
    vp_draw (xpos, ypos - (boxsz[1] / 2));
}

static void feedpage (void)
{
    vp_erase ();		/* feed a page, and reset defaults */
    vp_clip (0., 0., 8.0, 11.0);/* reset our clipping window */
    vp_tfont (1, VP_STROKE, 0);	/* reset our font */
    vp_fat ((int) (charsz / 2));/* reset our fatness */
    vp_orig (0., 0.);		/* reset our origin */
    vpos = 10.3;		/* start at the top of the page */
}
