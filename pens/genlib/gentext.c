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
 *  source file:   ./filters/genlib/gentext.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */
/*
 * Joe Dellinger Oct 18 1987
 * 	Keep track of fatness as a float, not an int. Round when needed.
 * Joe Dellinger Jan 16 1988
 *	Allow user-defined fonts. As a check, require that they have a
 *	magic sequence on the front. If they ask for a font >= NUMGENFONT,
 *	modulo back into range.
 * Joe Dellinger Feb 16 1988
 *	Make number of arguments to dev.attributes consistent.
 * Joe Dellinger May 4 1988
 *	Explicitly check to make sure the ERROR font is compile-time
 *	loaded.
 * Joe Dellinger May 5 1990
 *	Shorts aren't worth the savings in file size for runtime-loaded
 *	fonts. Just make everything ints.
 * Joe Dellinger Oct 26 1991
 *	Fall back on looking for vplot binary font files in
 *	SYSTEM_FONT_DIRECTORY only if the environmental variable
 *	VPLOTFONTDIR is not set. Note both of these will normally
 *	be set to something ending in a "/"!
 *	Also, don't reject a file with a bad magic number out of
 *	hand, check to see if the problem is byte-swapping. If it
 *	is, just go ahead and silently fix it.
 */

/*
 * VPLOT soft text plotting
 *
 * Keywords: vplot text vector hershey font
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <limits.h>

#include <rsfplot.h>

#include	"../include/extern.h"
#include	"../include/err.h"
#include	"../include/enum.h"
#include	"../include/params.h"
#include	"../include/font_definitions.h"
#include	"../include/attrcom.h"
#include	"../include/round.h"

#include "../utilities/util.h"

#include "gentext.h"

#define		NMARK   8	/* Maximum number of marks to use */
#define		MAXPOLY 100	/* Maximum number of points in polygon */

/*
 * The font to use for undefined fonts.
 * It should NOT be a runtime-loaded font!
 */
#define ERRFONT		0
/*
 * The glyph to use for undefined glyphs.
 * It must be a glyph in the font ERRFONT.
 * Needless to say, this glyph at least had better exist or you're
 * in real trouble.
 */
/* We use Glyph 30 in font 0, which is a '?' with a square around it. */
#define ERRGLYPH 	(30-font[0].dim[START])

/* Fraction of height of capital letter to use as vertical padding for box */
#define VSPACE_FRAC  .20
/* Fraction of inter-letter space to use for horizontal padding of box */
#define HSPACE_FRAC  0.5

#define EOC 0x8000	/* END OF CHARACTER BIT */
#define DRAWBIT 0x4000	/* DRAW BIT */
#define POLYBITS (EOC | DRAWBIT)	/* Polygon */
#define XBIT 0x0040
#define YBIT 0x2000
#define UNDEFINED	-1

#define SIGN_X (glyph_stroke & XBIT)
#define SIGN_Y (glyph_stroke & YBIT)
#define DRAW (glyph_stroke & DRAWBIT)
#define INPOLY ((glyph_stroke & POLYBITS) == POLYBITS)

#define COLOR_MAP(A) color_set[A][MAP];

/*
 * Correct for difference in size between what you want and what you've got
 */
#define SIZE_FACTOR(A) (((double)tsize/100.)*(A)/(double)(font[ttxfont].dim[CAP]-font[ttxfont].dim[BASE]))

/*
 * When you change sizes or fonts in midstream, what level do you line up on?
 */
#define ALIGN_HEIGHT	(font[ttxfont].dim[BASE])

#define CONTROL 001
#define BS	010
#define CR	015
#define NL	012

static double   xorigin_f, yorigin_f, xold_f, yold_f;
static int      ttxfont, cur_color_save;
extern int      cur_color, ipat, need_devcolor;
extern bool overlay;
static bool overlay_save;
extern int      color_set[MAX_COL + 1][_NUM_PRIM];

extern void drawpolygon (int npts, int *x, int *y);

static void load_font (int ifont);
static void mov (double hadd, double vadd);

double   path_orient_dx, path_orient_dy;
double   up_orient_dx, up_orient_dy;

#include <stdlib.h>

static void swab_font (char *stuff, int bytecount)
{
    int             icount;
    char            temp;

    for (icount = 0;
	 icount < bytecount - (int) sizeof (int) + 1;
	 icount += sizeof (int))
    {
	temp = stuff[icount + 0];
	stuff[icount + 0] = stuff[icount + 3];
	stuff[icount + 3] = temp;
	temp = stuff[icount + 1];
	stuff[icount + 1] = stuff[icount + 2];
	stuff[icount + 2] = temp;
    }
}

void gentext (char *string, float pathx, float pathy, float upx, float upy)
/*< interpret characters into vectors >*/
{
    double          fpathx, fpathy, fupx, fupy;
    double          up, path;
    double          xp, yp;
    float           tfat;
    int             add;
    int             ixp, iyp;
    int            *istring;
    unsigned int   *glyphptr, glyph_stroke;
    double          xtxshift, ytxshift;
    int             a, b, ii, jj, kk;
    int             string_length;
    double          last_widthl, last_widthr;
    double          widthl_1st_char, char_width, total_width = 0.;
    int             tsize, first, ttxfont_save=0, txfont_checked;
    int             ghost = 0;
    double          th_symb, th_symb_s=0, tv_symb, tv_symb_s=0;
    double          char_width_s[NMARK];
    double          total_width_s[NMARK];
    int             mark_flag[NMARK];
    double          last_widthl_s[NMARK], last_widthr_s[NMARK];
    double          xold_f_s[NMARK];
    double          yold_f_s[NMARK];
    double          maxtop, minbot;
    double          vspace, hspace, vline;
    int             flag;
    char           *charp;
    int             polycount, xxx[MAXPOLY], yyy[MAXPOLY];
    int            *ligp;
    static int      one_error = YES;
    int             linecount;

    if (*string == '\0')
	return;

/*
 * Set the initial parameters
 */

/*
 * Convert the input float path vectors to doubles.
 */
    fpathx = (double) pathx;
    fpathy = (double) pathy;
    fupx = (double) upx;
    fupy = (double) upy;

    path = hypot ((double) fpathx, (double) fpathy);
    up = hypot ((double) fupx, (double) fupy);

    if (path == 0. || up == 0.)
    {
/* Text got squashed away to nothing */
	return;
    }

    path_orient_dx = fpathx / path;
    path_orient_dy = fpathy / path;
    up_orient_dx = fupx / up;
    up_orient_dy = fupy / up;

/*
 * We didn't bomb out right away, so save things we may change so
 * they can be restored at exit
 */
    cur_color_save = cur_color;
    overlay_save = overlay;

    overlay = false;

    for (ii = 0; ii < NMARK; ii++)
	mark_flag[ii] = 0;

/* Text starts out in the default font at 100% of the requested size */
    tsize = 100;

    if (font[ERRFONT].load == NO)
    {
	ERR (FATAL,name,
	     "(gentext) Font %d, the ERRFONT, must be loaded at compile time!",
	     ERRFONT);
    }

    txfont_checked = dev.txfont;

    if (txfont_checked >= 0)
    {
	ttxfont = txfont_checked % NUMGENFONT;
	if (font[ttxfont].load == NO)
	{
	    load_font (txfont_checked);
	}
	txfont_checked = ttxfont;
    }
    else
    {
	txfont_checked = ERRFONT;
    }

    ttxfont = txfont_checked;
    total_width = 0.;
    char_width = font[ttxfont].dim[SPACE] * SIZE_FACTOR (path);
    th_symb = char_width / 2.;
    last_widthl = font[ttxfont].dim[SPACE] * .5 * SIZE_FACTOR (path);
    last_widthr = font[ttxfont].dim[SPACE] * .5 * SIZE_FACTOR (path);
    widthl_1st_char = font[ttxfont].dim[SPACE] * .5 * SIZE_FACTOR (path);
    tv_symb = (font[ttxfont].dim[HALF] - ALIGN_HEIGHT) * SIZE_FACTOR (up);
    maxtop = (font[ttxfont].dim[TOP] - ALIGN_HEIGHT) * SIZE_FACTOR (up);
    minbot = (font[ttxfont].dim[BOTTOM] - ALIGN_HEIGHT) * SIZE_FACTOR (up);

/* These used in making the bounding box */
    vline = (font[ttxfont].dim[TOP] - font[ttxfont].dim[BOTTOM] + font[ttxfont].dim[LINE]) * SIZE_FACTOR (up);
    vspace = VSPACE_FRAC * (font[ttxfont].dim[CAP] - font[ttxfont].dim[BASE]) * SIZE_FACTOR (up);
    hspace = (HSPACE_FRAC * font[ttxfont].dim[SPACE] + font[ttxfont].dim[LETTER]) * SIZE_FACTOR (path);

/*
 * Parse ligatures, control sequences, etc.
 * each object gets 2 ints;
 * The first tells what sort of object it is.
 * Positive == printing character
 * Negative == command
 * The second is there for any parameters that are associated.
 * For normal characters, it gives the glyph number in this font.
 * For commands it gives any parameter the command may have.
 */
    string_length = strlen (string);
    istring = (int *) calloc ((unsigned) 2 * (string_length + 1), sizeof (int));

    for (ii = 0, charp = string; (charp - string) < string_length; ii++, charp++)
    {
	switch ((int) (*charp))
	{
/* Check for special ASCII characters first */
	    case ' ':
	    case NL:
	    case CR:
	    case BS:
		istring[2 * ii] = -(int) (*charp);
		continue;
		break;
/* Check for \ commands */
	    case '\\':
		charp++;
		switch (*charp)
		{
/* \\ just falls through and makes a \ with no special properties */
		    case '\\':
			break;
/* \ commands with no arguments */
		    case '-':
		    case '>':
		    case '<':
		    case '^':
		    case '_':
		    case 'g':
		    case 'G':
		    case 'n':
		    case 'h':
			istring[2 * ii] = -(int) (*charp);
			continue;
			break;
/* \ commands with arguments */
		    case 's':
		    case 'f':
		    case 'F':
		    case 'k':
		    case 'r':
		    case 'm':
		    case 'M':
		    case 'v':
		    case 'c':
			istring[2 * ii] = -(int) (*charp);
			charp++;
/* default value of the argument is 0 if they just leave a space */
			istring[2 * ii + 1] = 0;
/* read the argument */
			sscanf (charp, "%d ", &istring[2 * ii + 1]);
/* skip past it and check for syntax */
			do
			{
			    if ((*charp >= '0' && *charp <= '9') ||
				*charp == '-' || *charp == '+')
			    {
				charp++;
			    }
			    else
				ERR (FATAL, name, "In text \\%c must be followed by an integer and then a space.", (char) (-istring[2 * ii]));
			} while (*charp != ' ');

			if (istring[2 * ii] == -(int) ('v'))
			{
/*
 * The \v command.
 * Make an ordinary character with the proper value.
 */
			    istring[2 * ii] = istring[2 * ii + 1];
			    istring[2 * ii + 1] = 0;
			}
			else
			    if (istring[2 * ii] == -(int) ('F'))
			    {
/* Font change command */
				if (istring[2 * ii + 1] >= 0)
				{
				    ttxfont = istring[2 * ii + 1] % NUMGENFONT;

/* Perform serif font substitution once and for all here, if needed */
				    if (! serifs_OK)
				    {
					if (ttxfont == DEFAULT_HARDCOPY_FONT) ttxfont = DEFAULT_SANSSERIF_FONT;
					if (ttxfont == GREEK_SERIF_FONT) ttxfont = GREEK_SANSSERIF_FONT;
				    }

/* On this first pass through, load all the fonts we're going to need */
				    if (font[ttxfont].load == NO)
				    {
					load_font (istring[2 * ii + 1]);
				    }
				    istring[2 * ii + 1] = ttxfont;
				}
				else
				{
/* \F-1  means the default font again. */
				    if (istring[2 * ii + 1] == -1)
					istring[2 * ii + 1] = txfont_checked;
				    else
					istring[2 * ii + 1] = ERRFONT;
				}
			    }
			    else
				if (istring[2 * ii] == -(int) ('c'))
				{
/* Color change command */
				    if (istring[2 * ii + 1] == -1)
				    {
/*
 * They want to return to the original text color.
 * This has already been checked to be within range and
 * properly mapped, so just use it!
 */
					istring[2 * ii + 1] = cur_color_save;
				    }
				    else
				    {
/*
 * Map from the color asked for to the colors that are available.
 * Normally only dovplot is allowed to do this, but this is an
 * unusual case where dovplot can't do it for us.
 */
					if (istring[2 * ii + 1] > MAX_COL || istring[2 * ii + 1] < 0)
					    ERR (FATAL, name, "(gentext) bad color number %d (max %d, min 0)",
						 istring[2 * ii + 1], MAX_COL);
					istring[2 * ii + 1] = COLOR_MAP (istring[2 * ii + 1]);
				    }
				}
			continue;
			break;
		    default:
			ERR (WARN, name, "(gentext) Unknown command \\%c.", *charp);
			charp--;
			break;
		}
	    default:
		break;
	}
/* Normal character */
	istring[2 * ii] = (int) (*charp);
    }
    string_length = ii;

/* Ligatures */
    if (dev.txprec > 1)
    {
	ttxfont = txfont_checked;
/*
 * Turning things into ligatures can only make the string shorter.
 * ii keeps track of where we are without ligatures,
 * kk keeps track of where we are with ligatures included.
 * The string is copied back into itself. Since ii >= kk, there is
 * no recursion problem.
 */
	for (ii = 0, kk = 0; ii < string_length; ii++, kk++)
	{
	    if (istring[2 * ii] < 0)
	    {
/*
 * The only special command we care about for constructing ligatures
 * is the font change command, since ligatures are font dependent.
 * The commands WILL break up a ligature, but other than that aren't
 * interpreted at all here.
 */
		if (-istring[2 * ii] == 'F')
		    ttxfont = istring[2 * ii + 1];
		istring[2 * kk] = istring[2 * ii];
		istring[2 * kk + 1] = istring[2 * ii + 1];
		continue;
	    }

/*
 * Take the first ligature that matches. This means that longer ligatures
 * MUST be listed first in the font data!
 */
/*
 * Loop over ligatures.
 * Each ligature has 1 number at the beginning giving the number of characters
 * in this ligature. The next number gives the glyph that is drawn for this
 * ligature. The next several numbers give the glyphs that are combined.
 * ligp points to the ligature we are currently searching for.
 * ligp += 2 + ligp[0] moves to the beginning of the next ligature:
 * 2 places for the 2 numbers at the beginning plus the ligp[0] characters
 * making up the ligature.
 */
	    for (ligp = font[ttxfont].lig; ligp[0] > 0; ligp += 2 + ligp[0])
	    {
/* Is there enough room before the end to possibly make this? */
		if (ii + ligp[0] - 1 < string_length)
		{
/* Loop over the characters in the ligature */
		    for (jj = 0; jj < ligp[0]; jj++)
		    {
/* Didn't match. Stop looking on this one. */
			if (ligp[jj + 2] != istring[2 * (ii + jj)])
			    goto failed;
		    }
/* Got to the end and so it worked. Put in the glyph for the ligature */
		    istring[2 * kk] = ligp[1];
/* skip past the ligp[0] characters in the original string that went into it */
		    ii += ligp[0] - 1;
		    goto success;
		}
	    failed:
		continue;
	    }
/* No ligatures for this one. Copy it across unchanged */
	    istring[2 * kk] = istring[2 * ii];
/*
 * Don't need to look at any more ligatures for this character
 * (Ligatures don't nest)
 */
	success:
	    continue;
	}
/* Update the length of the string */
	string_length = kk;
    }


/************************************************************************/
/************************************************************************/

/*
 * This section conducts a "dry run" through the text string to determine
 * its length.
 *
 * Each character has a left half-width (widthl) and a right half-width (widthr).
 * These give the left and right half-widths of the character's bounding box
 * away from the character's origin.
 * The vertical dimensions of each character's bounding box are a function of
 * the font only; font[font_number].dim[TOP] and font[font_number].dim[BOTTOM]
 * give the coordinates in font units of the top and bottom of the character's box.
 *
 * Each character is also separated from its neighbors by an inter-letter space
 * (font[font_number].dim[LETTER]). This is effectively tacked onto the beginning
 * of each new character, with the exception of the first.
 *
 * When we actually output vectors we will start with the center of the
 * 1st character as our origin. (So we have to remember the left half-width
 * in order to be able to compensate for the fact that we measure the length
 * of the string starting from the left hand edge of the first character.)
 * Variables like "char_width" keep track of the latest character's width
 * in case we have to back up over it again.
 *
 * Multiplying by SIZE_FACTOR converts from the Font's units to Vplot's.
 * (Horizontal scales with "path", vertical scales with "up". These simply
 * give the length of the path and up vectors.)
 * All variables except for those in the font structure itself are in Vplot units.
 *
 * Left and right the total width of everything (characters and inter-spaces
 * between them) is summed into total_width. This is used to do the horizontal
 * text justification.
 *
 * Up and down the highest top and lowest bottom to date are saved in "maxtop" and
 * "minbot". These are used for vertical text justification. "ALIGN_HEIGHT"
 * gives the effective origin for vertical glyph positioning.
 *
 * The "symb" variables keep track of the symbol position of the latest character.
 */
    first = 1;
    flag = 1;
    linecount = 1;
    ttxfont = txfont_checked;
    for (ii = 0; ii < string_length; ii++)
    {
/* 
 * Figure the lenth of the message string.
 * Justification is based on the first line of text only
 * (That's what "flag" is for)
 */

/*
 * Check for special characters
 */
	if (istring[2 * ii] < 0)
	{
	    switch (-istring[2 * ii])
	    {
		case 'n':
		case NL:
		    linecount++;
		case CR:
		    flag = 0;
		    break;
		case 'h':
		case BS:
		    total_width -= font[ttxfont].dim[LETTER] *
			SIZE_FACTOR (path) + char_width;
		    break;
		case 'F':
/* Change the font */
		    ttxfont = istring[2 * ii + 1];
		    break;
		case 's':
/* Change the size. This affects the SIZE_FACTOR. */
		    tsize = istring[2 * ii + 1];
		    break;
		case 'k':
		    if (flag)
		    {
/*
 * Add in the blank space created by horizontal 'k'earning.
 * This is measured in percent of the width of a space in this font.
 *
 * Similar vertical movements are ignored for the purposes of justification.
 */
			total_width += font[ttxfont].dim[SPACE]
			    * SIZE_FACTOR (path) * (istring[2 * ii + 1] / 100.);
		    }
		    break;
		case 'm':
		    if (istring[2 * ii + 1] < 0 || istring[2 * ii + 1] >= NMARK)
			ERR (FATAL, name,
			     "(gentext) Too high a mark number %d", istring[2 * ii + 1]);
/* Save all relevant parameters as they are at this instant */
		    if (flag)
		    {
/* Vertical symbol alignment position */
			tv_symb_s = tv_symb;
/* Horizontal symbol alignment position */
			th_symb_s = th_symb;
/* Width of this character (in case the next thing is a backspace) */
			char_width_s[istring[2 * ii + 1]] = char_width;
/* The width so far up to this point */
			total_width_s[istring[2 * ii + 1]] = total_width;
		    }
		    mark_flag[istring[2 * ii + 1]] = 1;
		    break;
		case 'M':
		    if (istring[2 * ii + 1] < 0 || istring[2 * ii + 1] >= NMARK)
			ERR (FATAL, name,
			     "(gentext) Too high a mark number %d", istring[2 * ii + 1]);
/* Make sure it isn't junk */
		    if (!mark_flag[istring[2 * ii + 1]])
			ERR (FATAL, name,
			     "(gentext) Attempt to use undefined mark number %d",
			     istring[2 * ii + 1]);
/*
 * Restore the parameters previously saved. All events after that point
 * are now ignored for the purposes of justification.
 */
		    if (flag)
		    {
			tv_symb = tv_symb_s;
			th_symb = th_symb_s;
			char_width = char_width_s[istring[2 * ii + 1]];
			total_width = total_width_s[istring[2 * ii + 1]];
		    }
		    break;
		case '-':		/* Nothing */
		    break;
		case '>':		/* Forward one inter-letter space */
		    if (flag)
			total_width += font[ttxfont].dim[LETTER]
			    * SIZE_FACTOR (path);
		    break;
		case '<':		/* Remove one inter-letter space */
		    if (flag)
			total_width -= font[ttxfont].dim[LETTER]
			    * SIZE_FACTOR (path);
		    break;
		case '^':		/* Up a half letter */
		case '_':		/* Down a half letter */
		case 'g':		/* Make text invisible */
		case 'G':		/* Make text visible again */
		    break;
		case ' ':		/* Space */
		    if (flag)
		    {
			char_width = font[ttxfont].dim[SPACE]
			    * SIZE_FACTOR (path);
			th_symb = char_width / 2.;
			tv_symb = (font[ttxfont].dim[HALF] - ALIGN_HEIGHT)
			    * SIZE_FACTOR (up);
			total_width += char_width;
			if (first)
			{
/* If it is the first character, remember the left half width */
			    widthl_1st_char = font[ttxfont].dim[SPACE] * .5
				* SIZE_FACTOR (path);
/* No longer at the first printable character */
			    first = 0;
			}
			else
			{
/* else add inter-letter space between it and the previous character */
			    total_width += font[ttxfont].dim[LETTER]
				* SIZE_FACTOR (path);
			}
		    }
		    break;
	    }
	    continue;
	}

	if (flag)
	{
/*
 * There are 2 ways a glyph can be undefined: it can be outside the range of
 * the font, OR it can have no data associated with it
 */
	    if (istring[2 * ii] >= font[ttxfont].dim[START] && istring[2 * ii] <= font[ttxfont].dim[END])
	    {
/* Find the glyph number, and save it so we don't have to recalculate it */
		istring[2 * ii + 1] = istring[2 * ii] - font[ttxfont].dim[START];
		if (font[ttxfont].saddr[istring[2 * ii + 1]] != UNDEFINED)
		{
/* OK glyph */
/* In case it's the last one, save its vertical symbol position */
		    tv_symb = (font[ttxfont].symbol[istring[2 * ii + 1]] - ALIGN_HEIGHT)
			* SIZE_FACTOR (up);
/* And in case we back up later its width */
		    char_width = (font[ttxfont].swidthl[istring[2 * ii + 1]] +
				  font[ttxfont].swidthr[istring[2 * ii + 1]])
			* SIZE_FACTOR (path);
/* and horizontal symbol position */
		    th_symb = font[ttxfont].swidthr[istring[2 * ii + 1]]
			* SIZE_FACTOR (path);
/* See if it sets a new record high */
		    if ((font[ttxfont].dim[TOP] - ALIGN_HEIGHT) * SIZE_FACTOR (up) > maxtop)
			maxtop = (font[ttxfont].dim[TOP] - ALIGN_HEIGHT) * SIZE_FACTOR (up);
/* Or a record low */
		    if ((font[ttxfont].dim[BOTTOM] - ALIGN_HEIGHT)
			* SIZE_FACTOR (up) < minbot)
			minbot = (font[ttxfont].dim[BOTTOM] - ALIGN_HEIGHT)
			    * SIZE_FACTOR (up);
/* Add it into the total width */
		    total_width += char_width;
		    if (first)
		    {
/* If it's the first remember its left half width */
			widthl_1st_char = font[ttxfont].swidthl[istring[2 * ii + 1]]
			    * SIZE_FACTOR (path);
		    }
		    else
		    {
/* or if not first add in the space between it and the previous glyph */
			total_width += font[ttxfont].dim[LETTER]
			    * SIZE_FACTOR (path);
		    }
		}
		else
		{
/* Second way to be undefined. Turn it into a "special" character */
		    istring[2 * ii] = UNDEFINED;
		}
	    }
	    else
	    {
/* First way to be undefined. Turn it into a "special" character */
		istring[2 * ii] = UNDEFINED;
	    }

	    if (istring[2 * ii] == UNDEFINED)
	    {
/*
 * If it is undefined, use the special "ERROR" glyph and then
 * treat that just like we would treat a regular character.
 */
		ttxfont_save = ttxfont;
		ttxfont = ERRFONT;
		istring[2 * ii + 1] = ERRGLYPH;
		char_width = (font[ttxfont].swidthl[ERRGLYPH] +
			      font[ttxfont].swidthr[ERRGLYPH])
		    * SIZE_FACTOR (path);
		th_symb = font[ttxfont].swidthr[ERRGLYPH] * SIZE_FACTOR (path);
		tv_symb = (font[ttxfont].symbol[ERRGLYPH] - ALIGN_HEIGHT)
		    * SIZE_FACTOR (up);
		if ((font[ttxfont].dim[TOP] - ALIGN_HEIGHT) * SIZE_FACTOR (up) > maxtop)
		    maxtop = (font[ttxfont].dim[TOP] - ALIGN_HEIGHT) * SIZE_FACTOR (up);
		if ((font[ttxfont].dim[BOTTOM] - ALIGN_HEIGHT) * SIZE_FACTOR (up) < minbot)
		    minbot = (font[ttxfont].dim[BOTTOM] - ALIGN_HEIGHT) * SIZE_FACTOR (up);
		total_width += char_width;
		if (first)
		{
		    widthl_1st_char = font[ttxfont].swidthl[ERRGLYPH] * SIZE_FACTOR (path);
		}
		else
		    total_width += font[ttxfont].dim[LETTER] * SIZE_FACTOR (path);
		ttxfont = ttxfont_save;
	    }

/* We printed something, so we aren't at the first character anymore */
	    first = 0;
	}
	else
	{
/*
 * If we're past the first line of text, do the few things that aren't related
 * to justification (looking for undefined glyphs, finding the glyph numbers)
 */
	    if (istring[2 * ii] >= font[ttxfont].dim[START] && istring[2 * ii] <= font[ttxfont].dim[END])
	    {
		istring[2 * ii + 1] = istring[2 * ii] - font[ttxfont].dim[START];
		if (font[ttxfont].saddr[istring[2 * ii + 1]] == UNDEFINED)
		{
		    istring[2 * ii] = UNDEFINED;
		}
	    }
	    else
	    {
		istring[2 * ii] = UNDEFINED;
	    }
	    if (istring[2 * ii] == UNDEFINED)
	    {
		istring[2 * ii + 1] = ERRGLYPH;
	    }
	}
    }

/*
 *  Set the proper alignment from the calculated length of the 
 *  text string. Remember that when we plot zero will be in the center
 *  of the first character and not the left hand edge of the first character,
 *  so we have to use widthl_1st_char to compensate for that.
 */
    switch (txalign.hor)
    {
	case TH_SYMBOL:
	    xtxshift = total_width - widthl_1st_char - th_symb;
	    break;
	case TH_CENTER:
	    xtxshift = total_width / 2. - widthl_1st_char;
	    break;
	case TH_RIGHT:
	    xtxshift = total_width - widthl_1st_char;
	    break;
	case TH_NORMAL:
	case TH_LEFT:
	default:
	    xtxshift = -widthl_1st_char;
	    break;
    }

    tsize = 100;
    ttxfont = txfont_checked;
    tfat = fat;

/*
 * CAP, HALF, and BASE are calculated based on font and size of default
 * TOP and BOTTOM are based on highest TOP and lowest BOTTOM of all
 * glyphs in string.
 */
    switch (txalign.ver)
    {
	case TV_SYMBOL:
	    ytxshift = tv_symb;
	    break;
	case TV_TOP:
	    ytxshift = maxtop;
	    break;
	case TV_CAP:
	    ytxshift = (font[ttxfont].dim[CAP] - ALIGN_HEIGHT) * SIZE_FACTOR (up);
	    break;
	case TV_HALF:
	    ytxshift = (font[ttxfont].dim[HALF] - ALIGN_HEIGHT) * SIZE_FACTOR (up);
	    break;
	case TV_BOTTOM:
	    ytxshift = minbot;
	    break;
	case TV_NORMAL:
	case TV_BASE:
	default:
	    ytxshift = (font[ttxfont].dim[BASE] - ALIGN_HEIGHT) * SIZE_FACTOR (up);
	    break;
    }


/************************************************************************/
/************************************************************************/


/*
 * This part of the code draws the characters.
 *
 * The complexity arises because when we do each character we have to
 * be at its CENTER in the left-right direction and at its ALIGN_HEIGHT
 * in the up-down direction.
 * So to move from one character to the next we have to move over by
 * the right half-width of the previous character, an inter-letter space,
 * and then the left half-width of the character we're on!
 * This helps to make the code confusing.
 *
 * We also have to always be ready to unexpectedly back up over the last
 * character, so we also have to keep around the left half-width of the
 * previous character as well.
 */
    xold_f = xold;
    yold_f = yold;

/*
 * The "mov" routine moves us around in the cock-eyed coordinate system
 * determined by the up and path vectors. There's no problem if these
 * vectors aren't orthogonal!
 */
    mov (-xtxshift, -ytxshift);
    xorigin_f = xold_f;
    yorigin_f = yold_f;

    if (dev.txovly)
    {
	mov (-widthl_1st_char - hspace, minbot - vline * (linecount - 1) - vspace);
	xxx[0] = ROUND (xold_f);
	yyy[0] = ROUND (yold_f);
	mov (0., maxtop - minbot + vline * (linecount - 1) + 2. * vspace);
	xxx[1] = ROUND (xold_f);
	yyy[1] = ROUND (yold_f);
	mov (total_width + 2. * hspace, 0.);
	xxx[2] = ROUND (xold_f);
	yyy[2] = ROUND (yold_f);
	mov (0., -(maxtop - minbot + vline * (linecount - 1) + 2. * vspace));
	xxx[3] = ROUND (xold_f);
	yyy[3] = ROUND (yold_f);

	if (dev.txovly == 2 || dev.txovly == 3)
	{
	    if (cur_color != 0 || need_devcolor)
	    {
		cur_color = 0;
		dev.attributes (SET_COLOR, cur_color, 0, 0, 0);
		need_devcolor = NO;
	    }
	    drawpolygon (4, xxx, yyy);
	    if (cur_color != cur_color_save || need_devcolor)
	    {
		cur_color = cur_color_save;
		dev.attributes (SET_COLOR, cur_color, 0, 0, 0);
		need_devcolor = NO;
	    }
	}

	if (dev.txovly == 1 || dev.txovly == 3)
	{
	    dev.vector (xxx[0], yyy[0], xxx[1], yyy[1], ROUND (tfat), 0);
	    dev.vector (xxx[1], yyy[1], xxx[2], yyy[2], ROUND (tfat), 0);
	    dev.vector (xxx[2], yyy[2], xxx[3], yyy[3], ROUND (tfat), 0);
	    dev.vector (xxx[3], yyy[3], xxx[0], yyy[0], ROUND (tfat), 0);
	}

	xold_f = xorigin_f;
	yold_f = yorigin_f;
    }

    first = 1;

/*
 * This is where the actual drawing of the characters takes place.
 * Loop over all the characters in the string
 */
    for (ii = 0; ii < string_length; ii++)
    {
/*
 * Check for special characters first
 */
	if (istring[2 * ii] < 0)
	{
	    switch (-istring[2 * ii])
	    {			/* standard carriage controls */
		case 'h':
		case BS:
		    mov (-font[ttxfont].dim[LETTER] * SIZE_FACTOR (path) -
			 (last_widthl + last_widthr), 0.);
		    break;
		case NL:
		case 'n':
		    xold_f = xorigin_f;
		    yold_f = yorigin_f;
		    mov (0., -(
			     (font[ttxfont].dim[TOP] - font[ttxfont].dim[BOTTOM] + font[ttxfont].dim[LINE])
			     * SIZE_FACTOR (up)));
		    xorigin_f = xold_f;
		    yorigin_f = yold_f;
		    first = 1;
		    break;
		case CR:
		    xold_f = xorigin_f;
		    yold_f = yorigin_f;
		    first = 1;
		    break;
		case 'F':
		    ttxfont = istring[2 * ii + 1];
		    break;
		case 'c':
		    if (cur_color != istring[2 * ii + 1] || need_devcolor)
		    {
			cur_color = istring[2 * ii + 1];
			dev.attributes (SET_COLOR, cur_color, 0, 0, 0);
			need_devcolor = NO;
		    }
		    break;
		case 's':
		    tsize = istring[2 * ii + 1];
		    break;
		case 'f':
		    tfat += istring[2 * ii + 1] * fatmult;
		    break;
		case 'k':
/* Horizontal motion */
		    mov (font[ttxfont].dim[SPACE] * SIZE_FACTOR (path) *
			 (istring[2 * ii + 1] / 100.), 0.);
		    break;
		case 'r':
/* Vertical motion */
		    mov (0., (font[ttxfont].dim[CAP] - font[ttxfont].dim[BASE])
			 * SIZE_FACTOR (up) * (istring[2 * ii + 1] / 100.));
		    break;
		case 'm':
/*
 * Save the current position. No need to check mark number valid; this
 * was checked when we did the justification
 */
		    last_widthl_s[istring[2 * ii + 1]] = last_widthl;
		    last_widthr_s[istring[2 * ii + 1]] = last_widthr;
		    xold_f_s[istring[2 * ii + 1]] = xold_f;
		    yold_f_s[istring[2 * ii + 1]] = yold_f;
		    break;
		case 'M':
/*
 * Restore the current position
 */
		    last_widthl = last_widthl_s[istring[2 * ii + 1]];
		    last_widthr = last_widthr_s[istring[2 * ii + 1]];
		    xold_f = xold_f_s[istring[2 * ii + 1]];
		    yold_f = yold_f_s[istring[2 * ii + 1]];
		    break;
		case 'G':
		    ghost = 0;
		    break;
		case 'g':
		    ghost = 1;
		    break;
		case '^':
/* Up half a character */
		    mov (0.,
			 (font[ttxfont].dim[CAP] - font[ttxfont].dim[BASE])
			 * (.5) * SIZE_FACTOR (up));
		    break;
		case '_':
/* Down half a character */
		    mov (0., -(
			     (font[ttxfont].dim[CAP] - font[ttxfont].dim[BASE])
			     * (.5) * SIZE_FACTOR (up)));
		    break;
		case '-':
		    break;
		case '>':
/* Right an inter-letter space */
		    mov (font[ttxfont].dim[LETTER] * SIZE_FACTOR (path), 0.);
		    break;
		case '<':
/* Left an inter-letter space */
		    mov (-(font[ttxfont].dim[LETTER] * SIZE_FACTOR (path)), 0.);
		    break;
		case -(UNDEFINED):
		    /* Don't overload them with error messages */
		    if (one_error)
		    {
			ERR (WARN, name,
			     "(gentext) Attempt(s) to use undefined glyph(s) in font %d",
			     ttxfont);
			one_error = NO;
		    }
/* Switch to use the ERROR glyph, and the treat it as a regular glyph */
		    ttxfont_save = ttxfont;
		    ttxfont = ERRFONT;
		    goto not_special;
		    break;
		case ' ':
		default:
		    if (!first)
		    {
			mov (
			    (font[ttxfont].dim[SPACE] * .5 +
			     font[ttxfont].dim[LETTER])
			    * SIZE_FACTOR (path)
			    + last_widthr
			    ,0.);
		    }
		    else
		    {
			first = 0;
		    }
		    last_widthl = font[ttxfont].dim[SPACE] * .5 * SIZE_FACTOR (path);
		    last_widthr = font[ttxfont].dim[SPACE] * .5 * SIZE_FACTOR (path);
		    break;
	    }
	}
	else
	{
	not_special:
/*
 *  Printable character.
 *  Pull out the actual strokes that make up each glyph.
 *  First get the address of the character from the address array
 *  Get the address by adding the offset (add) to the base array address
 *  (font[ttxfont].svec).
 */
	    add = font[ttxfont].saddr[istring[2 * ii + 1]];
	    glyphptr = font[ttxfont].svec + add;
/*
 * Now that we have the address of the fonts,
 * we position the pen at the beginning of the next character
 * (Unless it's the first printable character in which case we are already
 *  there)
 */
	    if (!first)
	    {
		mov (
		    (font[ttxfont].swidthl[istring[2 * ii + 1]] +
		     font[ttxfont].dim[LETTER])
		    * SIZE_FACTOR (path)
		    + last_widthr
		    ,0.);
	    }
	    else
		first = 0;

/* Save the left and right half-widths of this glyph */
	    last_widthl = font[ttxfont].swidthl[istring[2 * ii + 1]]
		* SIZE_FACTOR (path);
	    last_widthr = font[ttxfont].swidthr[istring[2 * ii + 1]]
		* SIZE_FACTOR (path);

/*
 * Calculate where to position each character in high precision.
 */
	    xnew = ROUND (xold_f);
	    ynew = ROUND (yold_f);
/*
 * This loop contains the structure for the actual drawing of the characters
 * We go through this block until an "END OF CHARACTER" is read
 *
 *  Strokes are kept in a packed format to save
 *  space. Each stroke is packed into an unsigned
 *  short int with the following format: 
 *
 *         edsyyyyyysxxxxxx
 *         ||||____|||____|--> The x-coordinate value
 *         |||   |  `--------> Set if X < 0			
 *         |||   `-----------> The y-coordinate value
 *         ||`---------------> Set if Y < 0
 *         |`----------------> Draw bit, set if
 *         |                    command is draw.
 *         |                    Clear, if move.
 *         `-----------------> End of Character Bit
 *
 *  This is enough bits per coordinate to accomodate all the "Hershey" fonts.
 *
 *  Polygons are also encoded into this scheme. If the EOC and DRAW bits
 *  are simultaneously on, then we are inside a polygon. The last point of
 *  the polygon will only have the draw flag on, and at that point the entire
 *  polygon will be outputted.
 *
 *  Unfortunately it was later discovered that using shorts and ints mixed
 *  together was a bad idea, so everything is now ints even though the
 *  top part of some ints are unused!!! So we get the worst of both worlds.
 */

	    polycount = 0;

	    while ((glyph_stroke = *glyphptr++) != EOC)
	    {
		a = glyph_stroke & 077;
		if (SIGN_X)
		    a = -a;
		b = (glyph_stroke >> 7) & 077;
		if (SIGN_Y)
		    b = -b;
		b -= ALIGN_HEIGHT;

/*
 * Here is the correct place to insert code to rotate a glyph.
 * You want to do that before it gets distorted by the global coordinate
 * transformation. The "ALIGN_HEIGHT" defines where the vertical origin
 * is for a glyph in a font. Note that we are in the font's coordinate
 * system units at this point.
 */
		xp = xold_f;
		yp = yold_f;
/*
 * Cock-eyed coordinate system.
 * "up" is in the direction of the up vector,
 * "right" is in the direction of the path vector.
 * These can be screwy, and thus distort the glyph.
 *
 * "path" and "up" contain the magnitude, and the
 * "orient_dx"'s and "orient_dy"'s contain the direction cosines
 */
		xp += SIZE_FACTOR (path) * a * path_orient_dx +
		    SIZE_FACTOR (up) * b * up_orient_dx;
		yp += SIZE_FACTOR (path) * a * path_orient_dy +
		    SIZE_FACTOR (up) * b * up_orient_dy;
		ixp = ROUND (xp);
		iyp = ROUND (yp);

		if (polycount > 0 && !INPOLY)
		{
/*
 * If we just WERE in a polygon, but are not now, then we must have just
 * finished one. Plot it out.
 */
		    xxx[polycount] = ixp;
		    yyy[polycount] = iyp;
		    polycount++;
		    if (!ghost)
		    {
			drawpolygon (polycount, xxx, yyy);
		    }
/*
 * Done with this one, reset the vertex counter
 */
		    polycount = 0;
		}
		else
		    if (INPOLY)
		    {
/*
 * We're still saving up the polygon. Save this vertex.
 */
			xxx[polycount] = ixp;
			yyy[polycount] = iyp;
			polycount++;
			if (polycount > MAXPOLY - 1)
			{
			    ERR (FATAL, name,
				 "(gentext) Too many points in polygon.");
			}
		    }
		    else
		    {
/* No polygons, just output the vector */
			if (DRAW && !ghost)
			    dev.vector (xnew, ynew, ixp, iyp, ROUND (tfat), 0);
		    }

		xnew = ixp;
		ynew = iyp;
	    }

/*
 * If we switched to the error font just to plot the error glyph,
 * then switch back to the correct font
 */
	    if (istring[2 * ii] == UNDEFINED)
	    {
		ttxfont = ttxfont_save;
	    }
	}
    }

/* Restore the correct color, if necessary */
    if (cur_color != cur_color_save)
    {
	cur_color = cur_color_save;
	need_devcolor = YES;
    }

/* Restore overlay mode */
    overlay = overlay_save;

/*
 * If they jump back into text they can continue right where they left off
 * (As long as they do nothing to move the pen between now and then.)
 */
    mov (last_widthr + font[ttxfont].dim[LETTER] * SIZE_FACTOR (path),
	 ytxshift);
    xold = ROUND (xold_f);
    yold = ROUND (yold_f);
    free ((char *) istring);
}

static void mov (double hadd, double vadd)
{
    xold_f += hadd * path_orient_dx + vadd * up_orient_dx;
    yold_f += hadd * path_orient_dy + vadd * up_orient_dy;
}

#define FONTCHECK_STRING "Vplot Binary fonT  \n"

static void load_font (int ifont)
{
    int             fd, length;
    char            filename[PATH_MAX], *fname;
    char            string[80];
    char           *newfont;
    int             offs[7];
    static int      done[NUMGENFONT];
    char           *stringptr;
    int             fontcheck;
    int             need_swab, file_ok;
    size_t          nread;
    int             vplotfontdirset;
    const char     *rsfroot;

    if (done[ttxfont])
    {
/*
 * Fonts are modulo NUMGENFONT. The first try in each slot determines
 * what happens. If the first try bombed, well, too bad, you get the
 * error font for all fonts that land there modulo NUMGENFONT forever more.
 */
	ttxfont = ERRFONT;
	return;
    }

/*
 * One way or another we're going to be done looking for this font
 * when we exit this routine.
 */
    done[ttxfont] = YES;


/*
 * Here is where we look for the font.
 * If the font is one of the "hardwired" ones, we look in the
 * SYSTEM_FONT_DIRECTORY for it by name [unless the environmental variable
 * "VPLOTFONTDIR" has been set in which case we look there instead].
 * If the font is not one of the "hardwired" ones, we look for
 * font number XX in the file "fontXX" in the current directory.
 *
 * In either case we also getpar for "fontXX" to see if the user has specified
 * an overriding file name on the command line.
 */
    if (ttxfont < NUM_FONTS) {
	vplotfontdirset = YES;

	if ((stringptr = getenv ("VPLOTFONTDIR")) != NULL) {
	    sprintf (filename, "%s%s.bin", stringptr, font[ttxfont].name);
	} else {
	    if ((rsfroot = getenv("RSFROOT")) == NULL) rsfroot="/usr";
	    sprintf (filename, "%s/include/%s.bin", rsfroot, font[ttxfont].name);
	}
    } else {
	vplotfontdirset = -1;
/*
 * The default place to look for a user-specified vplot binary font.
 */
	sprintf (filename, "./font%d.bin", ifont);
    }

/*
 * getpar parameter fontXX=name overrides all else. (Unless the font was
 * loaded at compile time in which case you don't get the opportunity.)
 */
    sprintf (string, "font%d", ifont);
    fname = sf_getstring(string);
    if (NULL != fname) {
	vplotfontdirset = -1;
    } else {
	fname = filename;
    }

/*
 * Try opening it.
 */
    if ((fd = open (fname, O_RDONLY)) == -1)
    {
/* Huh, it wasn't there. Complain and give up. */
	ERR (WARN, name,
	     "(gentext) Couldn't find font %d, file \"%s\".",
	     ttxfont, fname);
	if (vplotfontdirset == NO)
	{
	    ERR (COMMENT, name,
		 "(gentext) Perhaps you need to setenv VPLOTFONTDIR on this machine?");
	}
	else
	    if (vplotfontdirset == YES)
	    {
		ERR (COMMENT, name,
		     "(gentext) Is your environmental variable VPLOTFONTDIR correct?");
	    }
	    else
	    {
		ERR (COMMENT, name,
		     "(gentext) Is the command-line parameter font%d=\"%s\" correct?",
		     ifont, fname);
	    }
	ttxfont = ERRFONT;
	return;
    }


/*
 * At least there is a file there. But is it the right kind of file?
 * First check to make sure it has the magic sequence "Vplot Binary fonT  \n"
 * at the beginning. If it doesn't, it's junk, so don't read it in.
 */
    nread = read (fd, string, strlen (FONTCHECK_STRING));

    if (nread != strlen (FONTCHECK_STRING) ||
	strncmp (FONTCHECK_STRING, string, strlen (FONTCHECK_STRING)) != 0)
    {
	close (fd);
	ERR (WARN, name,
	     "(gentext) Font %d file \"%s\" is not a vplot font.",
	     ttxfont, fname);
	ttxfont = ERRFONT;
	return;
    }


/*
 * It begins promisingly... but is the opening string followed by
 * the binary integer FONTCHECK? Perhaps the file needs to be
 * byte-swapped?
 */
    file_ok = YES;
    need_swab = NO;

    nread = read (fd, (char *) &fontcheck, sizeof (int));

    if (nread != sizeof (int))
	file_ok = NO;
    else
    {
	if (fontcheck != FONTCHECK)
	{
	    swab_font ((char *) &fontcheck, (int) sizeof (int));
	    if (fontcheck != FONTCHECK)
		file_ok = NO;
	    else
	    {
		need_swab = YES;
#ifdef DEBUG
		ERR (WARN, name,
		     "(gentext) Font %d file \"%s\" needs byte-swapping.",
		     ttxfont, fname);
#endif
	    }
	}
    }

    if (file_ok == NO)
    {
    font_garbled:
/*
 * I guess if I were really tidy I would free up the memory allocated
 * to a font found to be truncated on disk. This way the partial
 * font will remain in memory but marked as "unloaded" since it was
 * found to be partially damaged. Well, I guess this way a manic
 * debugger can examine it in memory...
 */
	close (fd);
	ERR (WARN, name,
	     "(gentext) Font %d file \"%s\" is garbled or truncated.",
	     ttxfont, fname);
	ttxfont = ERRFONT;
	return;
    }

/*
 *  Suck the fonts into memory, byte-swapping as we go if indicated.
 *  The start of the file contains the length and then the
 *  offsets from the beginning to the 7 structures defining the font.
 */
    nread = read (fd, (char *) &length, sizeof (int));
    if (nread != sizeof (int))
	goto font_garbled;
    if (need_swab)
	swab_font ((char *) &length, (int) sizeof (int));

    newfont = sf_charalloc (length);

/* Offsets to the 7 structures defining the font... */
    nread = read (fd, (char *) offs, 7 * sizeof (int));
    if (nread != 7 * sizeof (int))
	goto font_garbled;
    if (need_swab)
	swab_font ((char *) offs, 7 * (int) sizeof (int));

/* The font data itself */
    nread = read (fd, (char *) newfont, length);
    if (nread != length)
	goto font_garbled;
    if (need_swab)
	swab_font (newfont, length);

/* Done! */
    close (fd);

/* The vital parameters (dimensions, bounds, and lengths) */
    font[ttxfont].dim = (int *) (newfont + offs[0]);
/* Pointers to the addresses of the glyphs themselves */
    font[ttxfont].saddr = (int *) (newfont + offs[1]);
/* Left widths */
    font[ttxfont].swidthl = (int *) (newfont + offs[2]);
/* Right widths */
    font[ttxfont].swidthr = (int *) (newfont + offs[3]);
/* Vertical symbol hot-spot position */
    font[ttxfont].symbol = (int *) (newfont + offs[4]);
/* The actual glyph strokes */
    font[ttxfont].svec = (unsigned int *) (newfont + offs[5]);
/* Ligature data */
    font[ttxfont].lig = (int *) (newfont + offs[6]);

/* Whether or not this font has been loaded into memory or not yet */
    font[ttxfont].load = YES;
}
