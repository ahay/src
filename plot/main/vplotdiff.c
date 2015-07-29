/* Vplot diff - see if 2 vplot files represent "identical" plots. */
/*
  Copyright (C) 2008 Joseph A. Dellinger
  
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
 * Author Joe Dellinger, 23 May 2008  -- 8 June 2008.
 * Vplotdiff is extensively modified from pldb.c in Madagascar,
 * which was in turn slightly modified from the pldb.c in
 * the 1994 revision of SEPlib. pldb.c was originally written
 * by Rob Clayton, then modified by Joe Dellinger and others.
 *
 * NOTE: The line numbers referenced in the output are with
 * respect to line numbers in the output of pldb. So if pldb
 * is changed, this program must be changed as well!
 *
 * PROBLEMS: Some of the changes vppen might make vplotdiff is
 * not smart enough to recognize as "not really changing the file".
 * In particular, vppen will omit a zero-length text string, which
 * vplotdiff will object to as a change, even though no difference
 * will be visible in the plot.
 * Also, whenever a pattern or color table is needed, vplotdiff
 * checks ALL the patterns and the COMPLETE color table, not just
 * the ones that are used. I think this is a feature and not a bug,
 * but others might disagree.
 * Should a text command leave the current pen position "undefined",
 * or set to "an unknown but defined location, i.e. the end of the
 * text"? I have it the latter, and it is legal to allow one text
 * command to follow another without a move in between.
 * -- Joe Dellinger, 7 June 2008.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include <rsf.h>
#include <rsfplot.h>
#include "../../pens/include/params.h"

/* These values must be well outside what is possible for a short integer */
#define VPLOTDIFF_NOT_INITIALIZED 99999
#define VPLOTDIFF_TEXT_INITIALIZED 99990

/* After this many errors stop complaining */
#define ENOUGH_ERRORS	10
#define SLOP_FACTOR 1.15 
#define MAX_TXT 40

#define RASTER_TOL 10 /* raster tolerance */
#define VECTOR_TOL 2  /* vector tolerance */

enum {NTC_NONEED, NTC_NEED, NTC_TEXT, NTC_EOF};

/*
 * All the attributes that define the current vplot
 * plotting state excepting those that are arrays.
 */
struct vplot_current_state
{
    int setstyle;
    int setstyle_count;
    int origin1;
    int origin2;
    int origin_count;
    int move1;
    int move2;
    int move_count;
    int txalign1;
    int txalign2;
    int txalign_count;
    int txfontprec1;
    int txfontprec2;
    int txfontprec3;
    int txfontprec_count;
    int overlay;
    int overlay_count;
    int color;
    int color_count;
    int fat;
    int fat_count;
    int window1;
    int window2;
    int window3;
    int window4;
    int window_count;
};
/* Plot attributes that are arrays */
struct color_table_current_state
{
    int largest;
    int count[NPAT];
    int red[NPAT];
    int green[NPAT];
    int blue[NPAT];
};
struct dash_state
{
    int count;
    int dashon;
    int dashes[MAXDASH * 2];
};
struct pattern_state
{
    int largest;
    int count[NPAT + 1];
    int nmul[NPAT + 1];
    int nx[NPAT + 1];
    int ny[NPAT + 1];
    int dim[NPAT + 1];
    int *patbits[NPAT + 1];
};
struct warn_state
{
    int setstyle;
    int origin;
    int window;
    int coltab;
    int move;
    int txalign;
    int txfontprec;
    int overlay;
    int color;
    int fat;
    int dash;
    int pattern;
};

static void text (const int debug, char *name, int count, FILE * stream);
/* Exactly like vp_getint, except reads from a stream instead of stdin. */
static int vp_getint2 (char *name, int count, FILE * stream);
static int vp_getint3 (FILE * stream);
/* Handle commands that set attributes but don't actually plot anything. */
static int check_vplot1 (const int debug,
			 char *name, int c, FILE * stream,
			 int *count, int *group,
			 struct vplot_current_state *state,
			 const struct vplot_current_state reset_state,
			 struct color_table_current_state *cstate,
			 struct dash_state *dstate,
			 struct pattern_state *pstate,
			 struct warn_state *warn,
			 struct warn_state *needtocheck, char *command);
/* Handle commands that plot things */
static int check_vplot2simple (int c1,
			       const int debug1,
			       char *name1, FILE * stream1, int *count1,
			       struct vplot_current_state *state1,
			       const int debug2,
			       char *name2, FILE * stream2, int *count2,
			       struct vplot_current_state *state2);
static int check_vplot2text (const int debug1, char *name1, int c1,
			     FILE * stream1, int *count1, const int debug2,
			     char *name2, int c2, FILE * stream2,
			     int *count2);
static int check_vplot2ras (const int debug1, char *name1, int c1,
			    FILE * stream1, int *count1,
			    struct color_table_current_state *cstate1,
			    const int debug2, char *name2, int c2,
			    FILE * stream2, int *count2,
			    struct color_table_current_state *cstate2);
/* For printing comprehensible error messages */
static void vplot_debug (int c, FILE * stream);
/* Checks that the vplot state matches */
static int check_state (const char *command1, const char *command2,
			const char *string1, const char *string2,
			const char *count1,
			const struct vplot_current_state *state1,
			const struct color_table_current_state *cstate1,
			const struct dash_state *dstate1,
			const struct pattern_state *pstate1,
			const char *count2,
			const struct vplot_current_state *state2,
			const struct color_table_current_state *cstate2,
			const struct dash_state *dstate2,
			const struct pattern_state *pstate2,
			struct warn_state *warn,
			struct warn_state *needtocheck1,
			struct warn_state *needtocheck2);

static const char *documentation[] = {
    "",
    "",
    "NAME",
    "		 sfvplotdiff -- (VPLOT DIFFerences) ",
    "SYNOPSIS",
    "		sfvplotdiff file1.vpl file2.vpl",
    "-----------------------------------------------------------------------",
    "	Reads two input vplot files; writes comments to standard error.",
    "   Return codes: ",
    "	 -1: There was a syntax error in one of the input vplot files.",
    "	     (Any vplot syntax error causes an immediate exit.)",
    "	  0: The two files represent equivalent plots.",
    "	  1: The two files differ in a vplot attribute.",
    "	  2: The two files differ in a plotting command.",
    "           (Such errors will cause an immediate exit.)",
    "",
    "",
    "",
    "SEE ALSO",
    "	 sfpldb, sfplas",
    "	 man vplot"
};

int
main (int argc, char *argv[])
{
    int c1, c2, i, ii, ndoc;
    FILE *stream1;
    FILE *stream2;
/* Starts at 1 because of comment from pldb added to starts of files. */
    int count1 = 1, count2 = 1;
    char count_string1[40], count_string2[40];
    char command1[40], command2[40];
    int error_count1 = 0;
    int group1 = 0, group2 = 0;

    struct pattern_state pstate1;
    struct pattern_state pstate2;
    struct dash_state dstate1;
    struct dash_state dstate2;
    struct color_table_current_state *cstate1;
    struct color_table_current_state *cstate2;
    struct vplot_current_state state1;
    struct vplot_current_state state2;
    const struct vplot_current_state reset_state = {
	/* style */ VPLOTDIFF_NOT_INITIALIZED, 0,
	/* origin */ VPLOTDIFF_NOT_INITIALIZED, VPLOTDIFF_NOT_INITIALIZED, 0,
	/* move */ VPLOTDIFF_NOT_INITIALIZED, VPLOTDIFF_NOT_INITIALIZED, 0,
	/* txalign */ 0, 0, 0,
	/* txfontprec */ 3, 2, 0, 0,
	/* overlay */ VPLOTDIFF_NOT_INITIALIZED, 0,
	/* color */ 7, 0,
	/* fat */ 0, 0,
	/* window */ -32760, -24570, 32760, 24570, 0
    };
    struct warn_state warn = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
/* setstyle, origin, window, and coltab start out and forever stay set. */
    struct warn_state needtocheck1 = { NTC_NEED, NTC_NEED, NTC_NEED, NTC_NEED,
				       NTC_NONEED, NTC_NONEED, NTC_NONEED, NTC_NONEED,
				       NTC_NONEED, NTC_NONEED,
				       NTC_NONEED, NTC_NONEED };
    struct warn_state needtocheck2 = { NTC_NEED, NTC_NEED, NTC_NEED, NTC_NEED,
				       NTC_NONEED, NTC_NONEED, NTC_NONEED, NTC_NONEED,
				       NTC_NONEED, NTC_NONEED,
				       NTC_NONEED, NTC_NONEED };

/*
 * If these are turned on, then on stdout you get a pldb-format
 * dump of the corresponding input file.
 */
    const int debug1 = 0, debug2 = 0;

    if (1 == argc && isatty (fileno (stdin)))
    {				/* no input - do selfdoc */
	ndoc = sizeof (documentation) / sizeof (documentation[0]);
	for (i = 0; i < ndoc; i++)
	{
	    printf ("%s\n", documentation[i]);
	}
	exit (0);
    }

    sf_init(argc,argv);

    if (argc < 3) sf_error("Need 2 input files.");

    if (NULL == (stream1 = fopen (argv[1], "r")))
	sf_error("Cannot open 1st file %s:", argv[1]);

    if (NULL == (stream2 = fopen (argv[2], "r")))
	sf_error("Cannot open 2nd file %s:", argv[2]);

/*
 *  Initialize plot state
 */
    state1 = reset_state;
    state2 = reset_state;

    cstate1 = (struct color_table_current_state *) sf_alloc(1,sizeof(*cstate1));
    cstate2 = (struct color_table_current_state *) sf_alloc(1,sizeof(*cstate2));

    cstate1->largest = -1;
    for (ii = 0; ii <= MAX_COL; ii++)
    {
	cstate1->count[ii] = 0;

	cstate1->red[ii] = VPLOTDIFF_NOT_INITIALIZED;
	cstate1->green[ii] = VPLOTDIFF_NOT_INITIALIZED;
	cstate1->blue[ii] = VPLOTDIFF_NOT_INITIALIZED;
    }

    cstate1->red[0] = 0;
    cstate1->green[0] = 0;
    cstate1->blue[0] = 0;

    cstate1->red[1] = 0;
    cstate1->green[1] = 0;
    cstate1->blue[1] = MAX_GUN;

    cstate1->red[2] = MAX_GUN;
    cstate1->green[2] = 0;
    cstate1->blue[2] = 0;

    cstate1->red[3] = MAX_GUN;
    cstate1->green[3] = 0;
    cstate1->blue[3] = MAX_GUN;

    cstate1->red[4] = 0;
    cstate1->green[4] = MAX_GUN;
    cstate1->blue[4] = 0;

    cstate1->red[5] = 0;
    cstate1->green[5] = MAX_GUN;
    cstate1->blue[5] = MAX_GUN;

    cstate1->red[6] = MAX_GUN;
    cstate1->green[6] = MAX_GUN;
    cstate1->blue[6] = 0;

    cstate1->red[7] = MAX_GUN;
    cstate1->green[7] = MAX_GUN;
    cstate1->blue[7] = MAX_GUN;

    cstate2->largest = -1;
    for (ii = 0; ii <= MAX_COL; ii++)
    {
	cstate2->count[ii] = 0;

	cstate2->red[ii] = cstate1->red[ii];
	cstate2->green[ii] = cstate1->green[ii];
	cstate2->blue[ii] = cstate1->blue[ii];
    }

    dstate1.dashon = 0;
    dstate2.dashon = 0;
    dstate1.count = 0;
    dstate2.count = 0;

    pstate1.largest = -1;
    pstate2.largest = -1;
    for (ii = 1; ii <= NPAT; ii++)
    {
	pstate1.count[ii] = 0;
	pstate2.count[ii] = 0;

	pstate1.nmul[ii] = 0;
	pstate1.nx[ii] = 0;
	pstate1.ny[ii] = 0;
	pstate1.dim[ii] = 0;
	pstate1.patbits[ii] = NULL;
	pstate2.nmul[ii] = 0;
	pstate2.nx[ii] = 0;
	pstate2.ny[ii] = 0;
	pstate2.dim[ii] = 0;
	pstate2.patbits[ii] = NULL;
    }

    strcpy (command1, "");
    strcpy (command2, "");

/*
 * Handle erases at the start of the file as a special case.
 */
    c1 = getc (stream1);
    c2 = getc (stream2);

    if (c1 == EOF) sf_error("%s is empty!", argv[1]);
    if (c2 == EOF) sf_error("%s is empty!", argv[2]);

/* Now we can be sure we won't try to push back an EOF. */
    if (c1 == VP_ERASE && c2 != VP_ERASE)
    {
	ungetc (c2, stream2);
/* Skip the extra initial erase in file 1 */
	count1++;
	if (debug1)
	    printf ("%c\n", c1);
    }
    else if (c2 == VP_ERASE && c1 != VP_ERASE)
    {
	ungetc (c1, stream1);
/* Skip the extra initial erase in file 2 */
	count2++;
	if (debug2)
	    printf ("%c\n", c2);
    }
    else
    {
/* No special trickery required; get on with it. */
	ungetc (c1, stream1);
	ungetc (c2, stream2);
    }

/*
 * Main loop
 *
 * Read file 1 until we get to a plotting command, then
 * Read file 2 until we get to a plotting command.
 * Check the vplot states are identical at that point.
 * If yes, then process next plotting command in parallel.
 * Check that that command matches.
 * "EOF" is treated like "another plotting command".
 */

    while (1)
    {
	/* Read file 1 until we get to a plotting command */
	while (1)
	{
	    c1 = getc (stream1);
	    if (0 != check_vplot1 (debug1, argv[1],
				   c1, stream1, &count1, &group1, &state1,
				   reset_state, cstate1, &dstate1,
				   &pstate1, &warn, &needtocheck1, command1))
		break;
	}
	/* Read file 2 until we get to a plotting command */
	while (1)
	{
	    c2 = getc (stream2);
	    if (0 != check_vplot1 (debug2, argv[2],
				   c2, stream2, &count2, &group2, &state2,
				   reset_state, cstate2, &dstate2,
				   &pstate2, &warn, &needtocheck2, command2))
		break;
	}

	/* Get the line numbers of the plotting commands (may be "EOF") */
	if (c1 == EOF)
	    strcpy (count_string1, "EOF");
	else
	    sprintf (count_string1, "%d", count1 + 1);

	if (c2 == EOF)
	    strcpy (count_string2, "EOF");
	else
	    sprintf (count_string2, "%d", count2 + 1);

	/*
	 * Both streams are now poised at plotting commands.
	 * Do their vplot states match?
	 * The needtocheck structures allow us to only check what's necessary.
	 */
	if (0 !=
	    check_state (command1, command2, argv[1], argv[2],
			 count_string1, &state1, cstate1, &dstate1,
			 &pstate1, count_string2, &state2, cstate2,
			 &dstate2, &pstate2, &warn, &needtocheck1, &needtocheck2))
	    error_count1++;

	if (error_count1 > ENOUGH_ERRORS) sf_error("maximum error count reached.");

/*
 * Both hit the end at the same time. We are done,
 * perhaps even successfully. This is the only
 * way to successfully exit.
 */
	if (c1 == EOF && c2 == EOF)
	    break;

/*
 * Check that plotting commands match
 */
	switch (c1)
	{
	    /* All the simple commands that can only be done one way */
	    case EOF:
	    case VP_DRAW:
	    case VP_PLINE:
	    case VP_PMARK:
	    case VP_AREA:
	    case VP_OLDAREA:
	    case VP_BACKGROUND:
	    case VP_ERASE:
		if (c1 != c2)
		{
		    sf_warning("command mismatch.");
		    fprintf (stderr, "\tFile %s, line %s, command: ",
			     argv[1], count_string1);
		    vplot_debug (c1, stream1);
		    fprintf (stderr, "\tFile %s, line %s, command: ",
			     argv[2], count_string2);
		    vplot_debug (c2, stream2);
		    exit (2);
		}
		error_count1 += check_vplot2simple (c1,
						    debug1, argv[1], stream1,
						    &count1, &state1,
						    debug2, argv[2], stream2,
						    &count2, &state2);
		break;
		/* Three different ways to make text */
	    case VP_GTEXT:
	    case VP_OLDTEXT:
	    case VP_TEXT:
		if (c2 != VP_GTEXT && c2 != VP_TEXT && c2 != VP_OLDTEXT)
		{
		    sf_warning ("command mismatch.");
		    fprintf (stderr, "\tFile %s, line %s, command: ",
			     argv[1], count_string1);
		    vplot_debug (c1, stream1);
		    fprintf (stderr, "\tFile %s, line %s, command: ",
			     argv[2], count_string2);
		    vplot_debug (c2, stream2);
		    exit (2);
		}
		error_count1 +=
		    check_vplot2text (debug1, argv[1], c1, stream1, &count1,
				      debug2, argv[2], c2, stream2, &count2);

                /*
		 * OK if the next command is another text command,
		 * not OK otherwise.
		 */
		state1.move1 = VPLOTDIFF_TEXT_INITIALIZED;
		state1.move2 = VPLOTDIFF_TEXT_INITIALIZED;
		state1.move_count = count1;
		state2.move1 = VPLOTDIFF_TEXT_INITIALIZED;
		state2.move2 = VPLOTDIFF_TEXT_INITIALIZED;
		state2.move_count = count2;
		break;
		/* Two different ways to make raster */
	    case VP_BYTE_RASTER:
	    case VP_BIT_RASTER:
		if (c2 != VP_BYTE_RASTER && c2 != VP_BIT_RASTER)
		{
		    sf_warning ("command mismatch.");
		    fprintf (stderr, "\tFile %s, line %s, command: ",
			     argv[1], count_string1);
		    vplot_debug (c1, stream1);
		    fprintf (stderr, "\tFile %s, line %s, command: ",
			     argv[2], count_string2);
		    vplot_debug (c2, stream2);
		    exit (2);
		}
		error_count1 +=
		    check_vplot2ras (debug1, argv[1], c1, stream1, &count1,
				     cstate1,
				     debug2, argv[2], c2, stream2, &count2, 
				     cstate2);

		state1.move1 = VPLOTDIFF_NOT_INITIALIZED;
		state1.move2 = VPLOTDIFF_NOT_INITIALIZED;
		state1.move_count = count1;
		state2.move1 = VPLOTDIFF_NOT_INITIALIZED;
		state2.move2 = VPLOTDIFF_NOT_INITIALIZED;
		state2.move_count = count2;
		break;
	    default:
		sf_error("%s, line %s: unknown command %c (%o octal)",
			 argv[1], count_string1, c1, c1);
		break;
	}

	if (error_count1 > ENOUGH_ERRORS)
	    sf_error ("maximum error count reached.");
    }

    if (error_count1 > 0)
	exit (1);
    else
	exit (0);
}

/*
 * Some commands have a text string following them.
 * This function reads those text strings from a stream,
 * checking for unexpected EOF's while doing so.
 */
static void
text (const int debug, char *name, int count, FILE * stream)
{
    int c;

    while ('\0' != (c = getc (stream)))
    {
	if (c == EOF)
	    sf_error("%s, line %d: Unexpected EOF in text!",name, count);
	    
	if (debug)
	    putchar (c);
    }
    if (debug)
	putchar ('\n');
}

/*
 * Extract and decode an integer from a stream.
 * This routine differs from the usual vp_getint
 * in that it does not assume the input is standard
 * input, and it checks for unexpected EOF's.
 */
int
vp_getint2 (char *name, int count, FILE * stream)
{
    int w1, w2;
    unsigned short w;

    w1 = getc (stream);
    if (w1 == EOF)
	sf_error ("%s, line %d: Unexpected EOF in integer!",name, count);

    w2 = getc (stream);
    if (w2 == EOF)
	sf_error("%s, line %d: Unexpected EOF in integer!",name, count);

    w = ((unsigned short) w1) + (((unsigned short) w2) << 8);

    return ((int) ((short) w));
}

/* Called from vplot_debug */
int
vp_getint3 (FILE * stream)
{
    int w1, w2;
    unsigned short w;

    w1 = getc (stream);
    if (w1 == EOF) sf_error("EOF");

    w2 = getc (stream);
    if (w2 == EOF) sf_error("EOF");

    w = ((unsigned short) w1) + (((unsigned short) w2) << 8);

    return ((int) ((short) w));
}

/*
 * This function handles vplot commands that set state variables
 * but don't actually draw anything. There is only one command
 * that both changes state variables and "plots" something:
 * the "erase" command.
 *
 * This routine will also catch any illegal commands.
 */
int
check_vplot1 (const int debug,
	      char *name, int c, FILE * stream, int *count, int *group,
	      struct vplot_current_state *state,
	      const struct vplot_current_state reset_state,
	      struct color_table_current_state *cstate,
	      struct dash_state *dstate, struct pattern_state *pstate,
	      struct warn_state *warn, struct warn_state *needtocheck,
	      char *command)
{
    int a, ix, iy, iz, npts, col_tab_no, i, n0;
    int ii, nmul, nx, ny;
    int ipat, col, j, bit;
    int x, y, red, green, blue, fat;
    int xmin, ymin, xmax, ymax, off, rep;

/*
 * Setstyle, origin, window, and coltab are just left set as they effect
 * everything.
 * For the others, until we find out otherwise, assume they don't matter.
 */
    needtocheck->move = NTC_NONEED;
    needtocheck->txalign = NTC_NONEED;
    needtocheck->txfontprec = NTC_NONEED;
    needtocheck->overlay = NTC_NONEED;
    needtocheck->color = NTC_NONEED;
    needtocheck->fat = NTC_NONEED;
    needtocheck->dash = NTC_NONEED;
    needtocheck->pattern = NTC_NONEED;
    sprintf (command, " not set");

/*
 * In this routine we only process commands
 * that don't actually plot anything.
 * For the plotting commands, we just do
 * few simple sanity checks and note what
 * matters for that command.
 *
 * "EOF" acts like a plotting command here.
 */
    switch (c)
    {
/*
 * First, plotting commands.
 */
	case EOF:
	    /* Make sure no groups left open at End Of File */
	    if (*group != 0) 
		sf_error("%s: %d unclosed group(s) at EOF.",
			 name, *group);

/*
 * At the end of the file, check EVERYTHING.
 * Someone may concatenate this vplot file onto the start of another,
 * and who knows what they may use in that file. Potentially the entire
 * state may matter.
 * It is OK to leave the pen position undefined at the end of a file, so
 * long as it is done CONSISTENTLY between the two files. The three
 * possibilities are "it was left undefined", "it was left pointing at
 * the end of a text string", and "it was left set to some particular
 * coordinate position".
 */
	    needtocheck->move = NTC_EOF;
	    needtocheck->txalign = NTC_NEED;
	    needtocheck->txfontprec = NTC_NEED;
	    needtocheck->overlay = NTC_NEED;
	    needtocheck->color = NTC_NEED;
	    needtocheck->fat = NTC_NEED;
	    needtocheck->dash = NTC_NEED;
	    needtocheck->pattern = NTC_NEED;
	    strcpy (command, "");
	    return 1;
	    break;
	case VP_DRAW:
	    /*
	     * The draw command has its own special check
	     * that the current pen position matches.
	     * So turn OFF the generic check!
	     * Since this is so counterintuitive I am making
	     * it explicit here.
	     */
	    needtocheck->move = NTC_NONEED;

	    needtocheck->color = NTC_NEED;
	    needtocheck->fat = NTC_NEED;
	    needtocheck->dash = NTC_NEED;
	    sprintf (command, " (%c)", c);
	    return 1;
	    break;
	case VP_PLINE:
	    needtocheck->color = NTC_NEED;
	    needtocheck->fat = NTC_NEED;
	    needtocheck->dash = NTC_NEED;
	    sprintf (command, " (%c)", c);
	    return 1;
	    break;
	case VP_PMARK:
	    needtocheck->txfontprec = NTC_NEED;
	    needtocheck->color = NTC_NEED;
	    needtocheck->fat = NTC_NEED;
	    sprintf (command, " (%c)", c);
	    return 1;
	    break;
	case VP_GTEXT:
	case VP_OLDTEXT:
	case VP_TEXT:
	    /*
	     * This special value indicates the command asking for
	     * the check is a text command, which has special rules.
	     */
	    needtocheck->move = NTC_TEXT;
	    needtocheck->txalign = NTC_NEED;
	    needtocheck->txfontprec = NTC_NEED;
	    needtocheck->color = NTC_NEED;
	    needtocheck->fat = NTC_NEED;
	    sprintf (command, " (%c)", c);
	    return 1;
	    break;
	case VP_BIT_RASTER:
	case VP_BYTE_RASTER:
	    needtocheck->overlay = NTC_NEED;
	    sprintf (command, " (%c)", c);
	    return 1;
	    break;
	case VP_OLDAREA:
	    needtocheck->overlay = NTC_NEED;
	    needtocheck->color = NTC_NEED;
	    sprintf (command, " (%c)", c);
	    return 1;
	    break;
	case VP_AREA:
	    needtocheck->overlay = NTC_NEED;
	    needtocheck->color = NTC_NEED;
	    needtocheck->pattern = NTC_NEED;
	    sprintf (command, " (%c)", c);
	    return 1;
	    break;
	case VP_BACKGROUND:
	    /*
	     * This command depends on what color 0 is set to,
	     * but as we always check the entire color table
	     * whenever we plot anything anyway, no additional
	     * checking is required. Just flag it as a plotting
	     * command by returning 1.
	     * 
	     * Should we later decide not to always check coltab,
	     * then we'll need here:
	     * needtocheck->coltab = NTC_NEED;
	     */
	    return 1;
	    break;
/*
 * An erase is a plotting command that also changes the plot state.
 */
	case VP_ERASE:
	    (*count)++;
	    /* Make sure no groups left open at an erase */
	    if (*group != 0)
		sf_error("%s: %d unclosed group(s) at erase at line %d.",
			 name, *group, *count);

	    /* An erase command resets much of the plot state */
	    *state = reset_state;
	    state->setstyle_count = *count;
	    warn->setstyle = 1;
	    state->origin_count = *count;
	    warn->origin = 1;
	    state->move_count = *count;
	    warn->move = 1;
	    state->txalign_count = *count;
	    warn->txalign = 1;
	    state->txfontprec_count = *count;
	    warn->txfontprec = 1;
	    state->overlay_count = *count;
	    warn->overlay = 1;
	    state->color_count = *count;
	    warn->color = 1;
	    state->fat_count = *count;
	    warn->fat = 1;
	    state->window_count = *count;
	    warn->window = 1;
	    dstate->dashon = 0;
	    dstate->count = *count;
	    warn->dash = 1;
	    /* But, it also counts as a command that does something! */
	    (*count)--;
	    sprintf (command, " (%c)", c);
	    return 1;
	    break;
/*
 * Next, commands that we can (mostly) just ignore as they neither
 * change the vplot state nor plot anything.
 */
	case VP_MESSAGE:
	    if (debug)
		printf ("%c\n", c);
	    (*count)++;
	    (*count)++;
	    text (debug, name, *count, stream);
	    break;
	case VP_BREAK:
	case VP_PURGE:
	case VP_NOOP:
	    if (debug)
		printf ("%c\n", c);
	    (*count)++;
	    break;
	case VP_BEGIN_GROUP:
	    if (debug)
		printf ("%c\n", c);
	    (*count)++;
	    (*count)++;
	    text (debug, name, *count, stream);
	    (*group)++;
	    break;
	case VP_END_GROUP:
	    if (debug)
		printf ("%c\n", c);
	    (*count)++;
	    (*group)--;
	    if (*group < 0)
		sf_error("%s: No group to end at line %d!", name, *count);
	    break;
/*
 * Finally, commands that change the vplot state.
 */
	case VP_SETSTYLE:
	    a = getc (stream);
	    if (EOF == a)
		sf_error("%s: Unexpected EOF in setstyle!", name);
	    a = tolower (a);
	    if (debug)
		printf ("%c %c\n", c, a);
	    (*count)++;
	    state->setstyle = a;
	    state->setstyle_count = *count;
	    warn->setstyle = 1;
	    /* The setstyle command also resets the origin */
	    state->origin1 = reset_state.origin1;
	    state->origin2 = reset_state.origin2;
	    state->origin_count = *count;
	    warn->origin = 1;
	    break;
	case VP_ORIGIN:
	    (*count)++;
	    x = vp_getint2 (name, *count, stream);
	    y = vp_getint2 (name, *count, stream);
	    if (debug)
		printf ("%c %d %d\n", c, x, y);
	    state->origin1 = x;
	    state->origin2 = y;
	    state->origin_count = *count;
	    warn->origin = 1;
	    break;
	case VP_MOVE:
	    (*count)++;
	    x = vp_getint2 (name, *count, stream);
	    y = vp_getint2 (name, *count, stream);
	    if (debug)
		printf ("%c %d %d\n", c, x, y);
	    state->move1 = x;
	    state->move2 = y;
	    state->move_count = *count;
	    warn->move = 1;
	    break;
	case VP_TXALIGN:
	    (*count)++;
	    ix = vp_getint2 (name, *count, stream);
	    iy = vp_getint2 (name, *count, stream);
	    if (debug)
		printf ("%c %d %d\n", c, ix, iy);
	    state->txalign1 = ix;
	    state->txalign2 = iy;
	    state->txalign_count = *count;
	    warn->txalign = 1;
	    break;
	case VP_TXFONTPREC:
	    (*count)++;
	    ix = vp_getint2 (name, *count, stream);
	    iy = vp_getint2 (name, *count, stream);
	    iz = vp_getint2 (name, *count, stream);
	    if (debug)
		printf ("%c %d %d %d\n", c, ix, iy, iz);
	    state->txfontprec1 = ix;
	    state->txfontprec2 = iy;
	    state->txfontprec3 = iz;
	    state->txfontprec_count = *count;
	    warn->txfontprec = 1;
	    break;
	case VP_SETDASH:
	    (*count)++;
	    npts = vp_getint2 (name, *count, stream);
	    if (npts > MAXDASH)
		sf_error("%s: Dash pattern of %d parts is too complex.",
			 name, npts);
	    if (debug)
		printf ("%c %d\n", c, npts);
	    dstate->dashon = npts;
	    dstate->count = *count;
	    for (ii = 0; ii < npts; ii++)
	    {
		(*count)++;
		x = vp_getint2 (name, *count, stream);
		y = vp_getint2 (name, *count, stream);
		if (debug)
		    printf ("%d %d\n", x, y);
		dstate->dashes[2 * ii] = x;
		dstate->dashes[2 * ii + 1] = y;
	    }
	    warn->dash = 1;
	    break;
	case VP_OVERLAY:
	    (*count)++;
	    ix = vp_getint2 (name, *count, stream);
	    if (debug)
		printf ("%c %d\n", c, ix);
	    state->overlay = ix;
	    state->overlay_count = *count;
	    warn->overlay = 1;
	    break;
	case VP_COLOR:
	    (*count)++;
	    ix = vp_getint2 (name, *count, stream);
	    if (debug)
		printf ("%c %d\n", c, ix);
	    state->color = ix;
	    state->color_count = *count;
	    warn->color = 1;
	    break;
	case VP_SET_COLOR_TABLE:
	    (*count)++;
	    col_tab_no = vp_getint2 (name, *count, stream);
	    red = vp_getint2 (name, *count, stream);
	    green = vp_getint2 (name, *count, stream);
	    blue = vp_getint2 (name, *count, stream);
	    if (debug)
		printf ("%c %d %d %d %d\n", c, col_tab_no, red, green, blue);
	    if (col_tab_no > MAX_COL)
		sf_error("%s line %d: Illegal color table number %d",
			 name, *count, col_tab_no);

	    if (col_tab_no > cstate->largest)
	    {
		cstate->largest = col_tab_no;
	    }

	    if ((red   < 0) || (red   > MAX_GUN) ||
		(green < 0) || (green > MAX_GUN) || 
		(blue  < 0) || (blue  > MAX_GUN))
		sf_error("%s line %d: color table %d=(%d,%d,%d) out of range",
			 name, *count, col_tab_no, red, green, blue);

	    cstate->red[col_tab_no] = red;
	    cstate->green[col_tab_no] = green;
	    cstate->blue[col_tab_no] = blue;
	    cstate->count[col_tab_no] = *count;
	    warn->coltab = 1;
	    break;
	case VP_FAT:
	    (*count)++;
	    fat = vp_getint2 (name, *count, stream);
	    if (debug)
		printf ("%c %d\n", c, fat);
	    state->fat = fat;
	    state->fat_count = *count;
	    warn->fat = 1;
	    break;
	case VP_WINDOW:
	    (*count)++;
	    xmin = vp_getint2 (name, *count, stream);
	    ymin = vp_getint2 (name, *count, stream);
	    xmax = vp_getint2 (name, *count, stream);
	    ymax = vp_getint2 (name, *count, stream);
	    if (debug)
		printf ("%c %d %d %d %d\n", c, xmin, ymin, xmax, ymax);
	    state->window1 = xmin;
	    state->window2 = ymin;
	    state->window3 = xmax;
	    state->window4 = ymax;
	    state->window_count = *count;
	    warn->window = 1;
	    break;
	case VP_PATLOAD:
	    (*count)++;
	    nmul = vp_getint2 (name, *count, stream);
	    nx = vp_getint2 (name, *count, stream);
	    ny = vp_getint2 (name, *count, stream);
	    ipat = vp_getint2 (name, *count, stream);
	    if (debug)
		printf ("%c %d %d %d %d\n", c, nmul, nx, ny, ipat);

	    if (ipat < 1 || ipat > NPAT)
		sf_error("%s: ipat = %d out of bounds.",
			 name, ipat);

	    if (ipat > pstate->largest)
	    {
		pstate->largest = ipat;
	    }
	    pstate->nmul[ipat] = nmul;
	    pstate->nx[ipat] = nx;
	    pstate->ny[ipat] = ny;
	    pstate->count[ipat] = *count;
	    if (pstate->patbits[ipat] != NULL)
	    {
		free ((char *) pstate->patbits[ipat]);
		pstate->patbits[ipat] = NULL;
	    }

	    if (-1 == nx)
	    {
		pstate->dim[ipat] = ny * 2 * 4;
		n0 = (pstate->dim[ipat] > 0) ? pstate->dim[ipat] : 1;
		if ((pstate->patbits[ipat] = (int *) malloc ((unsigned) n0 *
							     sizeof (int))) == NULL)
		    sf_error("%s: Cannot allocate memory for pattern %d.",
			     name, ipat);
		
		for (i = 0; i < ny * 2; i++)
		{
		    (*count)++;
		    fat = vp_getint2 (name, *count, stream);
		    col = vp_getint2 (name, *count, stream);
		    off = vp_getint2 (name, *count, stream);
		    rep = vp_getint2 (name, *count, stream);
		    if (debug)
			printf ("%d %d %d %d\n", fat, col, off, rep);
		    pstate->patbits[ipat][4 * i + 0] = fat;
		    pstate->patbits[ipat][4 * i + 1] = col;
		    pstate->patbits[ipat][4 * i + 2] = off;
		    pstate->patbits[ipat][4 * i + 3] = rep;
		}
	    }
	    else
	    {
		pstate->dim[ipat] = nx * ny;
		n0 = (pstate->dim[ipat] > 0) ? pstate->dim[ipat] : 1;
		pstate->patbits[ipat] = sf_intalloc ((unsigned) n0);

		ii = 0;
		for (i = 0; i < nx; i++)
		{
		    (*count)++;
		    for (j = 0; j < ny; j++)
		    {
			bit = vp_getint2 (name, *count, stream);
			if (bit >= 0 && bit <= 9)
			{
			    if (debug)
				putchar (bit + '0');
			}
			else
			{
			    if (debug)
				putchar (bit + 'A' - 10);
			}
			pstate->patbits[ipat][ii++] = bit;
		    }
		    if (debug)
			putchar ('\n');
		}
	    }
	    warn->pattern = 1;
	    break;
	default:
	    sf_error("%s, line %d: unknown command %c (%o octal)",
		     name, (*count) + 1, c, c);
	    break;
    }

    return 0;
}

/*
 * Check that two "text" commands are equivalent.
 */
int
check_vplot2text (const int debug1,
		  char *name1, int c1, FILE * stream1, int *count1,
		  const int debug2,
		  char *name2, int c2, FILE * stream2, int *count2)
{
    int errorcount = 0, texterror = 0;
    int x1, y1, xup1, yup1;
    int x2, y2, xup2, yup2;
    int key1, size1, orient1;
    int key2, size2, orient2;
    char command1[80], command2[80];
    float xfpath1=0., yfpath1=0., xfup1=0., yfup1=0.;
    float xfpath2=0., yfpath2=0., xfup2=0., yfup2=0.;
    float fsize, forient;
    char textstring1[MAX_TXT], textstring2[MAX_TXT];
    int ii;
    
    switch (c1)
    {
	case VP_GTEXT:
	    (*count1)++;
	    x1 = vp_getint2 (name1, *count1, stream1);
	    y1 = vp_getint2 (name1, *count1, stream1);
	    xup1 = vp_getint2 (name1, *count1, stream1);
	    yup1 = vp_getint2 (name1, *count1, stream1);
	    sprintf (command1, "%c %d %d %d %d", c1, x1, y1, xup1, yup1);
	    if (debug1)
		printf ("%c %d %d %d %d\n", c1, x1, y1, xup1, yup1);
	    
	    xfpath1 = x1;
	    yfpath1 = y1;
	    xfup1 = xup1;
	    yfup1 = yup1;
	    break;
	case VP_OLDTEXT:
	    (*count1)++;
	    key1 = vp_getint2 (name1, *count1, stream1);
	    size1 = (key1 & 037);
	    orient1 = (key1 & 0140) >> 5;
	    orient1 = 90 * orient1;
	    c1 = VP_TEXT;
	    goto vp_text1;
	    break;
	case VP_TEXT:
	    (*count1)++;
	    size1 = vp_getint2 (name1, *count1, stream1);
	    orient1 = vp_getint2 (name1, *count1, stream1);
	vp_text1:
	    sprintf (command1, "%c %d %d", c1, size1, orient1);
	    if (debug1)
		printf ("%c %d %d\n", c1, size1, orient1);
	    
	    /*
	     * Convert to GKS-style path - up representation.
	     */
	    fsize = (float) TEXTVECSCALE *(float) RPERIN *
		(float) size1 / (float) TXPERIN;
	    forient = (3.141592654 / 180.) * (float) orient1;
	    xfpath1 = cos (forient) * fsize;
	    yfpath1 = sin (forient) * fsize;
	    xfup1 = -sin (forient) * fsize;
	    yfup1 = cos (forient) * fsize;
	    break;
    }
    
    switch (c2)
    {
	case VP_GTEXT:
	    (*count2)++;
	    x2 = vp_getint2 (name2, *count2, stream2);
	    y2 = vp_getint2 (name2, *count2, stream2);
	    xup2 = vp_getint2 (name2, *count2, stream2);
	    yup2 = vp_getint2 (name2, *count2, stream2);
	    sprintf (command2, "%c %d %d %d %d", c2, x2, y2, xup2, yup2);
	    if (debug2)
		printf ("%c %d %d %d %d\n", c2, x2, y2, xup2, yup2);
	    
	    xfpath2 = x2;
	    yfpath2 = y2;
	    xfup2 = xup2;
	    yfup2 = yup2;
	    break;
	case VP_OLDTEXT:
	    (*count2)++;
	    key2 = vp_getint2 (name2, *count2, stream2);
	    size2 = (key2 & 037);
	    orient2 = (key2 & 0140) >> 5;
	    orient2 = 90 * orient2;
	    c2 = VP_TEXT;
	    goto vp_text2;
	    break;
	case VP_TEXT:
	    (*count2)++;
	    size2 = vp_getint2 (name2, *count2, stream2);
	    orient2 = vp_getint2 (name2, *count2, stream2);
	vp_text2:
	    sprintf (command2, "%c %d %d", c2, size2, orient2);
	    if (debug2)
		printf ("%c %d %d\n", c2, size2, orient2);
	    
	    /*
	     * Convert to GKS-style path - up representation.
	     */
	    fsize = (float) TEXTVECSCALE *(float) RPERIN *
		(float) size2 / (float) TXPERIN;
	    forient = (3.141592654 / 180.) * (float) orient2;
	    xfpath2 = cos (forient) * fsize;
	    yfpath2 = sin (forient) * fsize;
	    xfup2 = -sin (forient) * fsize;
	    yfup2 = cos (forient) * fsize;
	    break;
    }
    
    if ((fabs (xfpath1 - xfpath2) > (float) TEXTVECSCALE * SLOP_FACTOR) ||
	(fabs (yfpath1 - yfpath2) > (float) TEXTVECSCALE * SLOP_FACTOR) ||
	(fabs (xfup1 - xfup2) > (float) TEXTVECSCALE * SLOP_FACTOR) ||
	(fabs (yfup1 - yfup2) > (float) TEXTVECSCALE * SLOP_FACTOR))
    {
	sf_warning ("Coordinate mismatch: text");
	fprintf (stderr,
		 "\tFile %s, line %d: command: %s\n", name1, *count1, command1);
	fprintf (stderr,
		 "\tFile %s, line %d: command: %s\n", name2, *count2, command2);
	
	errorcount++;
    }
    
/* Do the text strings match? */
    (*count1)++;
    (*count2)++;
    texterror = 0;
    ii = 0;
    c1 = -2;
    c2 = -2;
    textstring1[MAX_TXT - 1] = '\0';
    textstring2[MAX_TXT - 1] = '\0';
    
    while (1)
    {
	if (c1 != '\0')
	{
	    c1 = getc (stream1);
	    if (c1 == EOF)
		sf_error("%s, line %d: Unexpected EOF in text!",
			 name1, *count1);
	}
	
	if (c2 != '\0')
	{
	    c2 = getc (stream2);
	    if (c2 == EOF)
		sf_error("%s, line %d: Unexpected EOF in text!",
			 name2, *count2);
	}
	
	if (ii < MAX_TXT - 1)
	{
	    textstring1[ii] = c1;
	    textstring2[ii] = c2;
	}
	
	if (c1 == '\0' && c2 == '\0')
	{
	    break;
	}
	
	if (c1 != c2)
	{
	    texterror = 1;
	}
	
	if (debug1 && c1 != '\0')
	    putchar (c1);
	if (debug2 && c2 != '\0')
	    putchar (c2);
	
	ii++;
    }
    
    if (debug1)
	putchar ('\n');
    if (debug2)
	putchar ('\n');
    
    if (texterror > 0)
    {
	sf_warning ("Text string mismatch");
	fprintf (stderr, "\tFile %s, line %d: %s\n", name1, *count1, textstring1);
	fprintf (stderr, "\tFile %s, line %d: %s\n", name2, *count2, textstring2);
	
	errorcount++;
    }
    
    return errorcount;
}

/*
 * Check that two raster commands are equivalent.
 */
int
check_vplot2ras (const int debug1,
		 char *name1, int c1, FILE * stream1, int *count1,
		 struct color_table_current_state *cstate1,
		 const int debug2,
		 char *name2, int c2, FILE * stream2, int *count2,
		 struct color_table_current_state *cstate2)
{
    int byte1, i1;
    int ii1, num_byte1;
    int ras_orient1, ras_offset1, xpix1, ypix1, num_rep1, num_pat1, pos1;
    int ibyte1 = 0;
    int xll1, yll1;
    int xur1, yur1;
    int byte2, i2;
    int ii2, num_byte2;
    int ras_orient2, ras_offset2, xpix2, ypix2, num_rep2, num_pat2, pos2;
    int ibyte2 = 0;
    int xll2, yll2;
    int xur2, yur2;
    int errorcount = 0, rasmismatch = 0;
    unsigned char *red1, *green1, *blue1;
    int *line1;
    int icount, icolno;
    int iline, irep;
    
/*
 * The Madagascar version of libvplot has no "blast" option
 * (IE, "compressed" raster). Looks like Sergey took it out.
 * The pen filters still know how to read it, though.
 * I'm supporting it here to maintain compatibility with SEPlib
 * and my own historical vplot files!
 */
    (*count1)++;
    ras_orient1 = vp_getint2 (name1, *count1, stream1);
    ras_offset1 = vp_getint2 (name1, *count1, stream1);
    if (debug1)
	printf ("%c %d %d\n", c1, ras_orient1, ras_offset1);
    (*count2)++;
    ras_orient2 = vp_getint2 (name2, *count2, stream2);
    ras_offset2 = vp_getint2 (name2, *count2, stream2);
    if (debug2)
	printf ("%c %d %d\n", c2, ras_orient2, ras_offset2);
    
    if (ras_orient1 != ras_orient2)
    {
	sf_warning ("Command mismatch: Raster orientation");
	fprintf (stderr, "\tFile %s, line %d: %c %d %d\n", name1, *count1,
		 c1, ras_orient1, ras_offset1);
	fprintf (stderr, "\tFile %s, line %d: %c %d %d\n", name2, *count2,
		 c2, ras_orient2, ras_offset2);
	
	exit (2);
    }
    
    (*count1)++;
    xll1 = vp_getint2 (name1, *count1, stream1);
    yll1 = vp_getint2 (name1, *count1, stream1);
    if (debug1)
	printf ("%d %d\n", xll1, yll1);
    (*count2)++;
    xll2 = vp_getint2 (name2, *count2, stream2);
    yll2 = vp_getint2 (name2, *count2, stream2);
    if (debug2)
	printf ("%d %d\n", xll2, yll2);
    (*count1)++;
    xur1 = vp_getint2 (name1, *count1, stream1);
    yur1 = vp_getint2 (name1, *count1, stream1);
    if (debug1)
	printf ("%d %d\n", xur1, yur1);
    (*count2)++;
    xur2 = vp_getint2 (name2, *count2, stream2);
    yur2 = vp_getint2 (name2, *count2, stream2);
    if (debug2)
	printf ("%d %d\n", xur2, yur2);
    
    if ((abs (xll1 - xll2) > 1) || (abs (yll1 - yll2) > 1) ||
	(abs (xur1 - xur2) > 1) || (abs (yur1 - yur2) > 1))
    {
	sf_warning ("Coordinate mismatch: Raster");
	fprintf (stderr,
		 "\tFile %s, line %d: (%d,%d,%d,%d)\n",
		 name1, *count1 - 1, xll1, yll1, xur1, yur1);
	fprintf (stderr,
		 "\tFile %s, line %d: (%d,%d,%d,%d)\n",
		 name2, *count2 - 1, xll2, yll2, xur2, yur2);
	
	errorcount++;
    }
    
    (*count1)++;
    xpix1 = vp_getint2 (name1, *count1, stream1);
    ypix1 = vp_getint2 (name1, *count1, stream1);
    if (debug1)
	printf ("%d %d\n", xpix1, ypix1);
    (*count2)++;
    xpix2 = vp_getint2 (name2, *count2, stream2);
    ypix2 = vp_getint2 (name2, *count2, stream2);
    if (debug2)
	printf ("%d %d\n", xpix2, ypix2);
    
    if ((xpix1 != xpix2) || (ypix1 != ypix2))
    {
	sf_warning ("Command mismatch: Raster size");
	fprintf (stderr,
		 "\tFile %s, line %d: %d %d\n", name1, *count1, xpix1, ypix1);
	fprintf (stderr,
		 "\tFile %s, line %d: %d %d\n", name2, *count2, xpix2, ypix2);
	
	exit (2);
    }
    
/*
 * Read the two rasters into memory as RGB triplets.
 */

    red1   = sf_ucharalloc (xpix1 * ypix1);
    green1 = sf_ucharalloc (xpix1 * ypix1);
    blue1 =  sf_ucharalloc (xpix1 * ypix1);
    line1 = sf_intalloc (xpix1 * ypix1);
    
    for (i1 = 0; i1 < ypix1; i1 += num_rep1)
    {
	(*count1)++;
	num_rep1 = vp_getint2 (name1, *count1, stream1);
	if (debug1)
	    printf ("     %d  ------ lines %d to %d\n",
		    num_rep1, i1, i1 + num_rep1 - 1);
	for (pos1 = 0; pos1 < xpix1; pos1 += num_byte1 * num_pat1)
	{
	    (*count1)++;
	    num_pat1 = vp_getint2 (name1, *count1, stream1);
	    num_byte1 = vp_getint2 (name1, *count1, stream1);
	    if (num_pat1 < 0 || num_byte1 < 0)
		sf_error("%s, line %d: Error in raster field",
			 name1, *count1);

	    if (debug1)
		printf ("%d %d  -- start byte %d in y_line %d\n",
			num_pat1, num_byte1, pos1, i1);

	    if (VP_BYTE_RASTER == c1)
	    {
		for (ii1 = 0; ii1 < num_byte1; ii1++)
		{
		    (*count1)++;
		    byte1 = getc (stream1);
		    if (EOF == byte1)
			sf_error("%s, line %d: "
				 "Unexpected EOF in byte raster!",
				 name1, *count1);

		    if (debug1)
			printf ("%d\n", byte1);
		    
		    icolno = ras_offset1 + byte1;
		    if (icolno < 0 || icolno > MAX_COL)
			sf_error("%s, line %d: "
				 "Attempt to use out of range color number in byte raster!",
				 name1, *count1);

		    if (cstate1->red[icolno] == VPLOTDIFF_NOT_INITIALIZED)
			sf_error("%s, line %d: "
				 "Attempt to use undefined color in byte raster!",
				 name1, *count1);
		    
		    for (iline = 0; iline < num_rep1; iline++)
		    {
			for (irep = 0; irep < num_pat1; irep++)
			{
			    icount =
				(i1 + iline) * xpix1 + irep * num_byte1 + pos1 + ii1;
			    
			    if (icount < 0 || icount >= xpix1 * ypix1)
				sf_error("%s, line %d: "
					 "Raster line not length promised!",
					 name1, *count1);

			    red1[icount] = cstate1->red[icolno];
			    green1[icount] = cstate1->green[icolno];
			    blue1[icount] = cstate1->blue[icolno];
			    line1[icount] = *count1;
			}
		    }
		}
	    }
	    else
	    {			/* VP_BIT_RASTER */
		for (ii1 = 0; ii1 < num_byte1; ii1++)
		{
		    (*count1)++;
		    if (0 == ii1 % 8)
		    {
			ibyte1 = getc (stream1);
			if (EOF == ibyte1)
			    sf_error("%s, line %d: "
				     "Unexpected EOF in bit raster!",
				     name1, *count1);
		    }
		    else
		    {
			ibyte1 <<= 1;
		    }
		    if (debug1)
			printf ("%d\n", ((ibyte1 >> 7) & 001));
		    
		    icolno = ras_offset1 * ((ibyte1 >> 7) & 001);
		    if (icolno < 0 || icolno > MAX_COL)
			sf_error("%s, line %d: "
				 "Attempt to use out of range color number in bit raster!",
				 name1, *count1);

		    if (cstate1->red[icolno] == VPLOTDIFF_NOT_INITIALIZED)
			sf_error("%s, line %d: "
				 "Attempt to use undefined color in bit raster!",
				 name1, *count1);
		    
		    for (iline = 0; iline < num_rep1; iline++)
		    {
			for (irep = 0; irep < num_pat1; irep++)
			{
			    icount =
				(i1 + iline) * xpix1 + irep * num_byte1 + pos1 + ii1;
		    
			    if (icount < 0 || icount >= xpix1 * ypix1)
				sf_error("%s, line %d: "
					 "Raster line not length promised!",
					 name1, *count1);
			    
			    red1[icount] = cstate1->red[icolno];
			    green1[icount] = cstate1->green[icolno];
			    blue1[icount] = cstate1->blue[icolno];
			    line1[icount] = *count1;
			}
		    }
		}
	    }
	}
    }
    
    for (i2 = 0; i2 < ypix2; i2 += num_rep2)
    {
	(*count2)++;
	num_rep2 = vp_getint2 (name2, *count2, stream2);
	if (debug2)
	    printf ("     %d  ------ lines %d to %d\n",
		    num_rep2, i2, i2 + num_rep2 - 1);
	for (pos2 = 0; pos2 < xpix2; pos2 += num_byte2 * num_pat2)
	{
	    (*count2)++;
	    num_pat2 = vp_getint2 (name2, *count2, stream2);
	    num_byte2 = vp_getint2 (name2, *count2, stream2);
	    if (num_pat2 < 0 || num_byte2 < 0)
		sf_error("%s, line %d: Error in raster field",
			 name2, *count2);

	    if (debug2)
		printf ("%d %d  -- start byte %d in y_line %d\n",
			num_pat2, num_byte2, pos2, i2);
	    
	    if (VP_BYTE_RASTER == c2)
	    {
		for (ii2 = 0; ii2 < num_byte2; ii2++)
		{
		    (*count2)++;
		    byte2 = getc (stream2);
		    if (EOF == byte2)
			sf_error("%s, line %d: Unexpected EOF in byte raster!",
				 name2, *count2);

		    if (debug2)
			printf ("%d\n", byte2);
		    
		    icolno = ras_offset2 + byte2;
		    if (icolno < 0 || icolno > MAX_COL)
			sf_error("%s, line %d: Attempt to use out of range color number in byte raster!",
				 name2, *count2);

		    if (cstate2->red[icolno] == VPLOTDIFF_NOT_INITIALIZED)
			sf_error("%s, line %d: Attempt to use undefined color in byte raster!",
				 name2, *count2);
		    
		    for (iline = 0; iline < num_rep2; iline++)
		    {
			for (irep = 0; irep < num_pat2; irep++)
			{
			    icount =
				(i2 + iline) * xpix2 + irep * num_byte2 + pos2 + ii2;
			    
			    if (icount < 0 || icount >= xpix2 * ypix2)
				sf_error("%s, line %d: Raster line not length promised!",
					 name2, *count2);

			    if ((abs
				 ((int) red1[icount] -
				  cstate2->red[icolno]) > RASTER_TOL)
				||
				(abs
				 ((int) green1[icount] -
				  cstate2->green[icolno]) > RASTER_TOL)
				||
				(abs
				 ((int) blue1[icount] - 
				  cstate2->blue[icolno]) > RASTER_TOL))
			    {
				if (rasmismatch == ENOUGH_ERRORS)
				{
				    sf_warning("Further raster color mismatches suppressed.");
				}
				else if (rasmismatch < ENOUGH_ERRORS)
				{
				    sf_warning ("Byte raster color mismatch");
				    fprintf (stderr,
					     "\tFile %s, line %d: color=(%d,%d,%d).\n",
					     name1,
					     line1[icount],
					     (int)
					     red1[icount],
					     (int) green1[icount], (int) blue1[icount]);
				    fprintf (stderr,
					     "\tFile %s, line %d: color=(%d,%d,%d).\n",
					     name2, *count2,
					     cstate2->
					     red[icolno],
					     cstate2->
					     green[icolno], cstate2->blue[icolno]);
				}
				
				rasmismatch++;
			    }
			}
		    }
		}
	    }
	    else
	    {			/* VP_BIT_RASTER */
		for (ii2 = 0; ii2 < num_byte2; ii2++)
		{
		    (*count2)++;
		    if (0 == ii2 % 8)
		    {
			ibyte2 = getc (stream2);
			if (EOF == ibyte2)
			    sf_error("%s, line %d: Unexpected EOF in bit raster!",
				     name2, *count2);
		    }
		    else
		    {
			ibyte2 <<= 1;
		    }
		    if (debug2)
			printf ("%d\n", ((ibyte2 >> 7) & 001));
		    
		    icolno = ras_offset2 * ((ibyte2 >> 7) & 001);
		    if (icolno < 0 || icolno > MAX_COL)
			sf_error("%s, line %d: Attempt to use out of range color number in bit raster!",
				 name2, *count2);

		    if (cstate2->red[icolno] == VPLOTDIFF_NOT_INITIALIZED)
			sf_error("%s, line %d: Attempt to use undefined color in bit raster!",
				 name2, *count2);
		    
		    for (iline = 0; iline < num_rep2; iline++)
		    {
			for (irep = 0; irep < num_pat2; irep++)
			{
			    icount =
				(i2 + iline) * xpix2 + irep * num_byte2 + pos2 + ii2;
			    
			    if (icount < 0 || icount >= xpix2 * ypix2)
				sf_error("s, line %d: Raster line not length promised!",
					 name2, *count2);

			    if ((abs
				 ((int) red1[icount] -
				  cstate2->red[icolno]) > 1)
				||
				(abs
				 ((int) green1[icount] -
				  cstate2->green[icolno]) > 1)
				||
				(abs
				 ((int) blue1[icount] - cstate2->blue[icolno]) > 1))
			    {
				if (rasmismatch == ENOUGH_ERRORS)
				{
				    sf_warning ("Further raster color mismatches suppressed.");
				}
				else if (rasmismatch < ENOUGH_ERRORS)
				{
				    sf_warning ("Bit raster color mismatch");
				    fprintf (stderr,
					     "\tFile %s, line %d: color=(%d,%d,%d).\n",
					     name1,
					     line1[icount],
					     (int) red1[icount],
					     (int) green1[icount], 
					     (int) blue1[icount]);
				    fprintf (stderr,
					     "\tFile %s, line %d: color=(%d,%d,%d).\n",
					     name2, *count2,
					     cstate2->red[icolno],
					     cstate2->green[icolno], 
					     cstate2->blue[icolno]);
				}
				
				rasmismatch++;
			    }
			}
		    }
		}
	    }
	}
    }
    
    free ((char *) red1);
    free ((char *) green1);
    free ((char *) blue1);
    free ((char *) line1);
    
    if (rasmismatch > 0)
	errorcount++;
    
    return errorcount;
}

/*
 * When we are about to draw something, check that the vplot states
 * in the two files match. By "state" I mean all the variables that
 * a vplot command that does not actually draw anything can set. These
 * variables then can effect how later commands are interpreted.
 * Most, but not all, of these are reset on an "erase" command.
 * Some, but not all, are defaulted to "sensible" values at the start
 * of the file and upon an erase.
 */
static int
check_state (const char *command1, const char *command2,
	     const char *file1, const char *file2,
	     const char *count1,
	     const struct vplot_current_state *state1,
	     const struct color_table_current_state *cstate1,
	     const struct dash_state *dstate1,
	     const struct pattern_state *pstate1,
	     const char *count2,
	     const struct vplot_current_state *state2,
	     const struct color_table_current_state *cstate2,
	     const struct dash_state *dstate2,
	     const struct pattern_state *pstate2,
	     struct warn_state *warn,
	     struct warn_state *needtocheck1, struct warn_state *needtocheck2)
{
    int ii, jj;
    char string1[120], string2[120];
    char string1a[80], string2a[80];
    char string1b[80], string2b[80];
    char string1c[80], string2c[80];
    int got_error = 0;
    int enough_already;
    
    if (warn->setstyle == 1)
    {
	if (state1->setstyle != state2->setstyle)
	{
	    if (state1->setstyle == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string1, "defaulted");
	    else
		sprintf (string1, "%c", state1->setstyle);
	    
	    if (state2->setstyle == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string2, "defaulted");
	    else
		sprintf (string2, "%c", state2->setstyle);
	    
	    sf_warning ("State mismatch: setstyle");
	    fprintf (stderr,
		     "\tFile %s, line %s%s: style=%s from line %d.\n",
		     file1, count1, command1, string1, state1->setstyle_count);
	    fprintf (stderr,
		     "\tFile %s, line %s%s: style=%s from line %d.\n",
		     file2, count2, command2, string2, state2->setstyle_count);
	    
	    got_error = 1;
	    warn->setstyle = 0;
	}
    }
    
    if (warn->origin == 1)
    {
	if ((abs (state1->origin1 - state2->origin1) > 1) ||
	    (abs (state1->origin2 - state2->origin2) > 1))
	{
	    if (state1->origin1 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string1, "defaulted");
	    else
		sprintf (string1, "%d", state1->origin1);
	    
	    if (state1->origin2 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string1a, "defaulted");
	    else
		sprintf (string1a, "%d", state1->origin2);
	    
	    if (state2->origin1 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string2, "defaulted");
	    else
		sprintf (string2, "%d", state2->origin1);
	    
	    if (state2->origin2 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string2a, "defaulted");
	    else
		sprintf (string2a, "%d", state2->origin2);
	    
	    sf_warning ("State mismatch: origin");
	    fprintf (stderr,
		     "\tFile %s, line %s%s: origin=(%s,%s) from line %d.\n",
		     file1, count1, command1, string1, string1a,
		     state1->origin_count);
	    fprintf (stderr,
		     "\tFile %s, line %s%s: origin=(%s,%s) from line %d.\n", file2,
		     count2, command2, string2, string2a, state2->origin_count);
	    
	    got_error = 1;
	    warn->origin = 0;
	}
    }

    if (warn->move == 1 && (needtocheck1->move != NTC_NONEED || 
			    needtocheck2->move != NTC_NONEED))
    {
	/*
	 * If needtocheck is set to NTC_TEXT it indicates it is a
	 * text command asking for the check.
	 * If the move value is VPLOTDIFF_TEXT_INITIALIZED
	 * it indicates the pen position was left unknown by
	 * a text command.
	 * It is OK to follow one text command immediately by another;
	 * the text just continues from where it left off.
	 *
	 * If needtocheck is set to NTC_EOF it indicates it is an
	 * EOF command asking for the check.
	 * At the EOF it is OK if the pen position is undefined,
	 * so long as it's undefined CONSISTENTLY.
	 */
	if (needtocheck1->move != NTC_EOF && needtocheck1->move != NTC_NONEED)
	{
	    if (((state1->move1 == VPLOTDIFF_NOT_INITIALIZED)
		 || (state1->move2 == VPLOTDIFF_NOT_INITIALIZED))
		|| (needtocheck1->move == NTC_NEED &&
		    ((state1->move1 == VPLOTDIFF_TEXT_INITIALIZED)
		     || (state1->move2 == VPLOTDIFF_TEXT_INITIALIZED))))
	    {
		sf_warning ("Coordinate problem: pen position undefined");
		fprintf (stderr,
			 "\tFile %s, line %s%s: position undefined from line %d.\n",
			 file1, count1, command1, state1->move_count);

		got_error = 1;
		warn->move = 0;
	    }
	}
	if (needtocheck2->move != NTC_EOF && needtocheck2->move != NTC_NONEED)
	{
	    if (((state2->move1 == VPLOTDIFF_NOT_INITIALIZED)
		 || (state2->move2 == VPLOTDIFF_NOT_INITIALIZED))
		|| (needtocheck2->move == NTC_NEED &&
		    ((state2->move1 == VPLOTDIFF_TEXT_INITIALIZED)
		     || (state2->move2 == VPLOTDIFF_TEXT_INITIALIZED))))
	    {
		sf_warning ("Coordinate problem: pen position undefined");
		fprintf (stderr,
			 "\tFile %s, line %s%s: position undefined from line %d.\n",
			 file2, count2, command2, state2->move_count);

		got_error = 1;
		warn->move = 0;
	    }
	}

	if ((abs (state1->move1 - state2->move1) > VECTOR_TOL) ||
	    (abs (state1->move2 - state2->move2) > VECTOR_TOL))
	{
	    if (state1->move1 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string1, "undefined");
	    else if (state1->move1 == VPLOTDIFF_TEXT_INITIALIZED)
		strcpy (string1, "end of text");
	    else
		sprintf (string1, "%d", state1->move1);

	    if (state1->move2 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string1a, "undefined");
	    else if (state1->move2 == VPLOTDIFF_TEXT_INITIALIZED)
		strcpy (string1a, "end of text");
	    else
		sprintf (string1a, "%d", state1->move2);

	    if (state2->move1 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string2, "undefined");
	    else if (state2->move1 == VPLOTDIFF_TEXT_INITIALIZED)
		strcpy (string2, "end of text");
	    else
		sprintf (string2, "%d", state2->move1);

	    if (state2->move2 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string2a, "undefined");
	    else if (state2->move2 == VPLOTDIFF_TEXT_INITIALIZED)
		strcpy (string2a, "end of text");
	    else
		sprintf (string2a, "%d", state2->move2);

	    sf_warning ("State mismatch: move");
	    fprintf (stderr,
		     "\tFile %s, line %s%s: move=(%s,%s) from line %d.\n",
		     file1, count1, command1, string1, string1a,
		     state1->move_count);
	    fprintf (stderr,
		     "\tFile %s, line %s%s: move=(%s,%s) from line %d.\n",
		     file2, count2, command2, string2, string2a,
		     state2->move_count);

	    got_error = 1;
	    warn->move = 0;
	}
    }

    if (warn->txalign == 1 &&
	(needtocheck1->txalign != NTC_NONEED || 
	 needtocheck2->txalign != NTC_NONEED))
    {
	if ((state1->txalign1 != state2->txalign1) ||
	    (state1->txalign2 != state2->txalign2))
	{
	    if (state1->txalign1 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string1, "defaulted");
	    else
		sprintf (string1, "%d", state1->txalign1);

	    if (state1->txalign2 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string1a, "defaulted");
	    else
		sprintf (string1a, "%d", state1->txalign2);

	    if (state2->txalign1 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string2, "defaulted");
	    else
		sprintf (string2, "%d", state2->txalign1);

	    if (state2->txalign2 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string2a, "defaulted");
	    else
		sprintf (string2a, "%d", state2->txalign2);

	    sf_warning ("State mismatch: txalign");
	    fprintf (stderr,
		     "\tFile %s, line %s%s: txalign=(%s,%s) from line %d.\n",
		     file1, count1, command1, string1, string1a,
		     state1->txalign_count);
	    fprintf (stderr,
		     "\tFile %s, line %s%s: txalign=(%s,%s) from line %d.\n",
		     file2, count2, command2, string2, string2a,
		     state2->txalign_count);

	    got_error = 1;
	    warn->txalign = 0;
	}
    }

    if (warn->txfontprec == 1 &&
	(needtocheck1->txfontprec != NTC_NONEED || 
	 needtocheck2->txfontprec != NTC_NONEED))
    {
	if ((state1->txfontprec1 != state2->txfontprec1) ||
	    (state1->txfontprec2 != state2->txfontprec2) ||
	    (state1->txfontprec3 != state2->txfontprec3))
	{
	    if (state1->txfontprec1 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string1, "defaulted");
	    else
		sprintf (string1, "%d", state1->txfontprec1);

	    if (state1->txfontprec2 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string1a, "defaulted");
	    else
		sprintf (string1a, "%d", state1->txfontprec2);

	    if (state1->txfontprec3 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string1b, "defaulted");
	    else
		sprintf (string1b, "%d", state1->txfontprec3);

	    if (state2->txfontprec1 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string2, "defaulted");
	    else
		sprintf (string2, "%d", state2->txfontprec1);

	    if (state2->txfontprec2 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string2a, "defaulted");
	    else
		sprintf (string2a, "%d", state2->txfontprec2);

	    if (state2->txfontprec3 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string2b, "defaulted");
	    else
		sprintf (string2b, "%d", state2->txfontprec3);

	    sf_warning ("State mismatch: txfontprec");
	    fprintf (stderr,
		     "\tFile %s, line %s%s: txfontprec=(%s,%s,%s) from line %d.\n",
		     file1, count1, command1, string1, string1a, string1b,
		     state1->txfontprec_count);
	    fprintf (stderr,
		     "\tFile %s, line %s%s: txfontprec=(%s,%s,%s) from line %d.\n",
		     file2, count2, command2, string2, string2a, string2b,
		     state2->txfontprec_count);

	    got_error = 1;
	    warn->txfontprec = 0;
	}
    }

    if (warn->overlay == 1 &&
	(needtocheck1->overlay != NTC_NONEED || 
	 needtocheck2->overlay != NTC_NONEED))
    {
	if ((state1->overlay != state2->overlay))
	{
	    if (state1->overlay == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string1, "defaulted");
	    else
		sprintf (string1, "%d", state1->overlay);

	    if (state2->overlay == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string2, "defaulted");
	    else
		sprintf (string2, "%d", state2->overlay);

	    sf_warning ("State mismatch: overlay");
	    fprintf (stderr,
		     "\tFile %s, line %s%s: overlay=%s from line %d.\n",
		     file1, count1, command1, string1, state1->overlay_count);
	    fprintf (stderr,
		     "\tFile %s, line %s%s: overlay=%s from line %d.\n",
		     file2, count2, command2, string2, state2->overlay_count);

	    got_error = 1;
	    warn->overlay = 0;
	}
    }

    if (warn->color == 1 &&
	(needtocheck1->color != NTC_NONEED || 
	 needtocheck2->color != NTC_NONEED))
    {
	if ((state1->color != state2->color))
	{
	    if (state1->color == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string1, "defaulted");
	    else
		sprintf (string1, "%d", state1->color);

	    if (state2->color == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string2, "defaulted");
	    else
		sprintf (string2, "%d", state2->color);

	    sf_warning ("State mismatch: color");
	    fprintf (stderr,
		     "\tFile %s, line %s%s: color=%s from line %d.\n",
		     file1, count1, command1, string1, state1->color_count);
	    fprintf (stderr,
		     "\tFile %s, line %s%s: color=%s from line %d.\n",
		     file2, count2, command2, string2, state2->color_count);

	    got_error = 1;
	    warn->color = 0;
	}
    }

    if (warn->fat == 1 && (needtocheck1->fat != NTC_NONEED || 
			   needtocheck2->fat != NTC_NONEED))
    {
	if ((state1->fat != state2->fat))
	{
	    if (state1->fat == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string1, "defaulted");
	    else
		sprintf (string1, "%d", state1->fat);

	    if (state2->fat == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string2, "defaulted");
	    else
		sprintf (string2, "%d", state2->fat);

	    sf_warning ("State mismatch: fat");
	    fprintf (stderr,
		     "\tFile %s, line %s%s: fat=%s from line %d.\n",
		     file1, count1, command1, string1, state1->fat_count);
	    fprintf (stderr,
		     "\tFile %s, line %s%s: fat=%s from line %d.\n",
		     file2, count2, command2, string2, state2->fat_count);

	    got_error = 1;
	    warn->fat = 0;
	}
    }

    if (warn->window == 1)
    {
	if ((abs (state1->window1 - state2->window1) > 1) ||
	    (abs (state1->window2 - state2->window2) > 1) ||
	    (abs (state1->window3 - state2->window3) > 1) ||
	    (abs (state1->window4 - state2->window4) > 1))
	{
	    if (state1->window1 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string1, "defaulted");
	    else
		sprintf (string1, "%d", state1->window1);

	    if (state1->window2 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string1a, "defaulted");
	    else
		sprintf (string1a, "%d", state1->window2);

	    if (state1->window3 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string1b, "defaulted");
	    else
		sprintf (string1b, "%d", state1->window3);

	    if (state1->window4 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string1c, "defaulted");
	    else
		sprintf (string1c, "%d", state1->window4);

	    if (state2->window1 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string2, "defaulted");
	    else
		sprintf (string2, "%d", state2->window1);

	    if (state2->window2 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string2a, "defaulted");
	    else
		sprintf (string2a, "%d", state2->window2);

	    if (state2->window3 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string2b, "defaulted");
	    else
		sprintf (string2b, "%d", state2->window3);

	    if (state2->window4 == VPLOTDIFF_NOT_INITIALIZED)
		strcpy (string2c, "defaulted");
	    else
		sprintf (string2c, "%d", state2->window4);

	    sf_warning ("State mismatch: window");
	    fprintf (stderr,
		     "\tFile %s, line %s%s: window=(%s,%s,%s,%s) from line %d.\n",
		     file1, count1, command1, string1, string1a, string1b,
		     string1c, state1->window_count);
	    fprintf (stderr,
		     "\tFile %s, line %s%s: window=(%s,%s,%s,%s) from line %d.\n",
		     file2, count2, command2, string2, string2a, string2b,
		     string2c, state2->window_count);

	    got_error = 1;
	    warn->window = 0;
	}
    }

    if (warn->dash == 1 && (needtocheck1->dash != NTC_NONEED || 
			    needtocheck2->dash != NTC_NONEED))
    {
	warn->dash = 0;

	if (dstate1->dashon != dstate2->dashon)
	{
	    if (0==dstate1->dashon) 
	    {
		got_error = 0;

		for (ii = 0; ii < dstate2->dashon; ii++)
		{
		    if (0 != dstate2->dashes[2 * ii] || 0 != dstate2->dashes[2 * ii + 1])
		    {
			got_error = 1;
			break;
		    }
		}
	    }
	    else if (0==dstate2->dashon) 
	    {
		got_error = 0;

		for (ii = 0; ii < dstate1->dashon; ii++)
		{
		    if (0 != dstate1->dashes[2 * ii] || 0 != dstate1->dashes[2 * ii + 1])
		    {
			got_error = 1;
			break;
		    }
		}
	    }
	    else
	    {
		got_error = 1;
	    }

	    if (0 != got_error) {
		sf_warning ("State mismatch: dash");
		fprintf (stderr,
			 "\tFile %s, line %s%s: pattern of %d dashes from line %d.\n",
			 file1, count1, command1, dstate1->dashon, dstate1->count);
		fprintf (stderr,
			 "\tFile %s, line %s%s: pattern of %d dashes from line %d.\n",
			 file2, count2, command2, dstate2->dashon, dstate2->count);
	    }
	}
	else
	{
	    for (ii = 0; ii < dstate1->dashon; ii++)
	    {
		if ((abs
		     (dstate1->dashes[2 * ii] -
		      dstate2->dashes[2 * ii]) > 1)
		    ||
		    (abs
		     (dstate1->dashes[2 * ii + 1] -
		      dstate2->dashes[2 * ii + 1]) > 1))
		{
		    sf_warning ("State mismatch: dash");
		    fprintf (stderr,
			     "\tFile %s, line %s%s: discrepancy in dash pattern from line %d.\n",
			     file1, count1, command1, dstate1->count);
		    fprintf (stderr,
			     "\tFile %s, line %s%s: discrepancy in dash pattern from line %d.\n",
			     file2, count2, command2, dstate2->count);

		    got_error = 1;
		    break;
		}
	    }
	}
    }

    if (warn->pattern == 1 &&
	(needtocheck1->pattern != NTC_NONEED || 
	 needtocheck2->pattern != NTC_NONEED))
    {
/*
 * LOTS of ways patterns can fail to match!
 */

/* Largest pattern is different; that's a mismatch right there */
	if (pstate1->largest != pstate2->largest)
	{
	    sf_warning ("State mismatch: pattern");
	    if (pstate1->largest < 0)
	    {
		fprintf (stderr,
			 "\tFile %s, line %s%s: no patterns set.\n", file1, count1,
			 command1);
	    }
	    else
	    {
		fprintf (stderr,
			 "\tFile %s, line %s%s: patterns set from line(s)",
			 file1, count1, command1);
		for (ii = 0; ii <= pstate1->largest; ii++)
		{
		    if (pstate1->patbits[ii] != NULL)
		    {
			fprintf (stderr, " %d", pstate1->count[ii]);
		    }
		}
		fprintf (stderr, ".\n");
	    }
	    if (pstate2->largest < 0)
	    {
		fprintf (stderr,
			 "\tFile %s, line %s%s: no patterns set.\n", file2, count2,
			 command2);
	    }
	    else
	    {
		fprintf (stderr,
			 "\tFile %s, line %s%s: patterns set from line(s)",
			 file2, count2, command2);
		for (ii = 0; ii <= pstate2->largest; ii++)
		{
		    if (pstate2->patbits[ii] != NULL)
		    {
			fprintf (stderr, " %d", pstate2->count[ii]);
		    }
		}
		fprintf (stderr, ".\n");
	    }

	    got_error = 1;
	    warn->pattern = 0;
	}
	else
	{
/*
 * Largest pattern is the same in both files.
 * Check them against each other pattern by pattern.
 */
	    for (ii = 0; ii <= pstate1->largest; ii++)
	    {
		if ((pstate1->patbits[ii] == NULL)
		    && (pstate2->patbits[ii] == NULL))
		{
/* Both patterns numbered ii are empty --- keep scanning. */
		    continue;
		}
		else if ((pstate1->patbits[ii] != NULL)
			 && (pstate2->patbits[ii] != NULL))
		{
/* Both patterns numbered ii have something. Do they match dimensions? */
		    if ((pstate1->dim[ii] != pstate2->dim[ii]) ||
			(pstate1->nmul[ii] != pstate2->nmul[ii]) ||
			(pstate1->nx[ii] != pstate2->nx[ii]) ||
			(pstate1->ny[ii] != pstate2->ny[ii]))
		    {
			sf_warning ("State mismatch: pattern");
			fprintf (stderr,
				 "\tFile %s, line %s%s: pattern %d from line %d.\n",
				 file1, count1, command1, ii, pstate1->count[ii]);
			fprintf (stderr,
				 "\tFile %s, line %s%s: pattern %d from line %d.\n",
				 file2, count2, command2, ii, pstate2->count[ii]);

			got_error = 1;
			warn->pattern = 0;
		    }
		    else
		    {
/* Dimensions match. Do the patterns match exactly? */

			for (jj = 0; jj < pstate1->dim[ii]; jj++)
			{
			    if (pstate1->patbits[ii][jj] != pstate2->patbits[ii][jj])
			    {
				sf_warning ("State mismatch: pattern");
				fprintf (stderr,
					 "\tFile %s, line %s%s: pattern %d from line %d.\n",
					 file1, count1, command1, ii,
					 pstate1->count[ii]);
				fprintf (stderr,
					 "\tFile %s, line %s%s: pattern %d from line %d.\n",
					 file2, count2, command2, ii,
					 pstate2->count[ii]);

				got_error = 1;
				warn->pattern = 0;
				break;
			    }
			}
		    }
		}
		else if (pstate1->patbits[ii] != NULL)
		{
/* Pattern 1 set, but the other is not. */
		    sf_warning ("State mismatch: pattern");
		    fprintf (stderr,
			     "\tFile %s, line %s%s: pattern %d from line %d set.\n",
			     file1, count1, command1, ii, pstate1->count[ii]);
		    fprintf (stderr,
			     "\tFile %s, line %s%s: pattern %d not set.\n",
			     file2, count2, command2, ii);

		    got_error = 1;
		    warn->pattern = 0;
		}
		else
		{
/* Pattern 2 set, but the other is not. */
		    sf_warning ("State mismatch: pattern");
		    fprintf (stderr,
			     "\tFile %s, line %s%s: pattern %d not set.\n",
			     file1, count1, command1, ii);
		    fprintf (stderr,
			     "\tFile %s, line %s%s: pattern %d from line %d set.\n",
			     file2, count2, command2, ii, pstate2->count[ii]);

		    got_error = 1;
		    warn->pattern = 0;
		}
	    }
	}
    }

/*
 * Currently needtocheck.coltab is initialized to NTC_NEED and never changed.
 * We'll leave the needtocheck check in for coltab, however,
 * in case someone later comes up with a case where
 * it doesn't need to be checked.
 */
    if (warn->coltab == 1 &&
	(needtocheck1->coltab != NTC_NONEED || 
	 needtocheck2->coltab != NTC_NONEED))
    {
	enough_already = 0;
	for (ii = 0; ii <= MAX_COL; ii++)
	{
	    /*
	     * No need to check any further, we're off the end
	     * of what's set on both.
	     */
	    if (ii > cstate1->largest && ii > cstate2->largest)
		break;

	    if ((abs (cstate1->red[ii] - cstate2->red[ii]) > 1) ||
		(abs (cstate1->green[ii] - cstate2->green[ii]) > 1) ||
		(abs (cstate1->blue[ii] - cstate2->blue[ii]) > 1))
	    {
		if (cstate1->count[ii] == 0)
		{
		    if (ii <= 7)
			sprintf (string1, "=(%d,%d,%d) (defaulted)",
				 cstate1->red[ii], cstate1->green[ii],
				 cstate1->blue[ii]);
		    else
			strcpy (string1, " never set");
		}
		else
		    sprintf (string1, "=(%d,%d,%d) from line %d",
			     cstate1->red[ii], cstate1->green[ii],
			     cstate1->blue[ii], cstate1->count[ii]);

		if (cstate2->count[ii] == 0)
		{
		    if (ii <= 7)
			sprintf (string2, "=(%d,%d,%d) (defaulted)",
				 cstate2->red[ii], cstate2->green[ii],
				 cstate2->blue[ii]);
		    else
			strcpy (string2, " never set");
		}
		else
		    sprintf (string2, "=(%d,%d,%d) from line %d",
			     cstate2->red[ii], cstate2->green[ii],
			     cstate2->blue[ii], cstate2->count[ii]);

		sf_warning ("State mismatch: color table");
		fprintf (stderr,
			 "\tFile %s, line %s%s: color table entry %d%s.\n",
			 file1, count1, command1, ii, string1);
		fprintf (stderr,
			 "\tFile %s, line %s%s: color table entry %d%s.\n",
			 file2, count2, command2, ii, string2);

		got_error = 1;
		warn->coltab = 0;

		enough_already++;
		if (enough_already >= ENOUGH_ERRORS)
		{
		    sf_warning ("Further color table mismatches suppressed.");
		    break;
		}
	    }
	}
    }

    return got_error;
}

/*
 * Used for printing comprehensible error messages.
 * (IE, giving error messages that match the output of pldb.)
 */
void
vplot_debug (int c, FILE * stream)
{
    int npts, mtype, key;
    int orient;
    int ras_orient, ras_offset;
    int x, y, xcor, ycor, msize, size;

    switch (c)
    {
	case EOF:
	    fprintf (stderr, "EOF\n");
	    break;
	case VP_DRAW:
	    x = vp_getint3 (stream);
	    y = vp_getint3 (stream);
	    fprintf (stderr, "%c %d %d\n", c, x, y);
	    break;
	case VP_PLINE:
	    npts = vp_getint3 (stream);
	    fprintf (stderr, "%c %d\n", c, npts);
	    break;
	case VP_PMARK:
	    npts = vp_getint3 (stream);
	    mtype = vp_getint3 (stream);
	    msize = vp_getint3 (stream);
	    fprintf (stderr, "%c %d %d %d\n", c, npts, mtype, msize);
	    break;
	case VP_GTEXT:
	    x = vp_getint3 (stream);
	    y = vp_getint3 (stream);
	    xcor = vp_getint3 (stream);
	    ycor = vp_getint3 (stream);
	    fprintf (stderr, "%c %d %d %d %d\n", c, x, y, xcor, ycor);
	    break;
	case VP_OLDTEXT:
	    key = vp_getint3 (stream);
	    size = (key & 037);
	    orient = (key & 0140) >> 5;
	    fprintf (stderr, "%c %d %d\n", VP_TEXT, size, 90 * orient);
	    break;
	case VP_TEXT:
	    size = vp_getint3 (stream);
	    orient = vp_getint3 (stream);
	    fprintf (stderr, "%c %d %d\n", c, size, orient);
	    break;
	case VP_BIT_RASTER:
	case VP_BYTE_RASTER:
	    ras_orient = vp_getint3 (stream);
	    ras_offset = vp_getint3 (stream);
	    fprintf (stderr, "%c %d %d\n", c, ras_orient, ras_offset);
	    break;
	case VP_AREA:
	    npts = vp_getint3 (stream);
	    fprintf (stderr, "%c %d\n", c, npts);
	    break;
	case VP_OLDAREA:
	    npts = vp_getint3 (stream);
	    fprintf (stderr, "%c %d\n", c, npts);
	    break;
	case VP_BACKGROUND:
	case VP_ERASE:
	    fprintf (stderr, "%c\n", c);
	    break;
	default:
	    fprintf (stderr, "unknown?! (This shouldn't happen.)\n");
	    break;
    }

    return;
}

/*
 * Check that two "simple" vplot commands
 * are equivalent. By "simple" I mean there is
 * only one way to accomplish that command!
 */
int
check_vplot2simple (int c,
		    const int debug1, char *name1, FILE * stream1,
		    int *count1,
		    struct vplot_current_state *state1,
		    const int debug2, char *name2, FILE * stream2,
		    int *count2, struct vplot_current_state *state2)
{
    int npts1, mtype1, i1;
    int ii1, xmask1, ymask1;
    int x1=0, y1=0, fat1, msize1;
    int npts2, mtype2;
    int xmask2, ymask2;
    int x2=0, y2=0, fat2, msize2;
    int errorcount = 0;

/* Only handle commands that plot something */

    switch (c)
    {
	case VP_DRAW:
	    (*count1)++;
	    x1 = vp_getint2 (name1, *count1, stream1);
	    y1 = vp_getint2 (name1, *count1, stream1);
	    if (debug1)
		printf ("%c %d %d\n", c, x1, y1);

	    (*count2)++;
	    x2 = vp_getint2 (name2, *count2, stream2);
	    y2 = vp_getint2 (name2, *count2, stream2);
	    if (debug2)
		printf ("%c %d %d\n", c, x2, y2);

/*
 * Is the "current pen position" set?
 */
	    if (((state1->move1 == VPLOTDIFF_NOT_INITIALIZED)
		|| (state1->move2 == VPLOTDIFF_NOT_INITIALIZED)) ||
		((state1->move1 == VPLOTDIFF_TEXT_INITIALIZED)
		 || (state1->move2 == VPLOTDIFF_TEXT_INITIALIZED)))
	    {
		sf_warning ("Coordinate problem: draw from unknown location");
		fprintf (stderr,
			 "\tFile %s, lines %d and %d: draw from (undefined,undefined) to (%d,%d)\n",
			 name1, state1->move_count, *count1, x1, y1);

		errorcount++;
	    }

	    if (((state2->move1 == VPLOTDIFF_NOT_INITIALIZED)
		|| (state2->move2 == VPLOTDIFF_NOT_INITIALIZED)) ||
		((state2->move1 == VPLOTDIFF_TEXT_INITIALIZED)
		 || (state2->move2 == VPLOTDIFF_TEXT_INITIALIZED)))
	    {
		sf_warning ("Coordinate problem: draw from unknown location");
		fprintf (stderr,
			 "\tFile %s, lines %d and %d: draw from (undefined,undefined) to (%d,%d)\n",
			 name2, state2->move_count, *count2, x2, y2);

		errorcount++;
	    }

	    if (errorcount > 0)
	    {
		if ((abs (x1 - x2) > VECTOR_TOL) || 
		    (abs (y1 - y2) > VECTOR_TOL))
		{
		    sf_warning ("Coordinate mismatch: draw");
		    fprintf (stderr,
			     "\tFile %s, line %d: draw to (%d,%d)\n",
			     name1, *count1, x1, y1);
		    fprintf (stderr,
			     "\tFile %s, line %d: draw to (%d,%d)\n",
			     name2, *count2, x2, y2);

		    errorcount++;
		}
	    }
	    else
	    {
/*
 * If we are closer swapped around, then check it that way.
 */
		if ((abs (x1 - x2) + abs (y1 - y2) +
		     abs (state1->move1 - state2->move1) + abs (state1->move2 -
								state2->move2)) <=
		    (abs (x1 - state2->move1) + abs (y1 - state2->move2) +
		     abs (state1->move1 - x2) + abs (state1->move2 - y2)))
		{
		    if ((abs (x1 - x2) > VECTOR_TOL) || 
			(abs (y1 - y2) > VECTOR_TOL) ||
			(abs (state1->move1 - state2->move1) > VECTOR_TOL) ||
			(abs (state1->move2 - state2->move2) > VECTOR_TOL))
		    {
			sf_warning ("Coordinate mismatch: draw");
			fprintf (stderr,
				 "\tFile %s, lines %d and %d: draw from (%d,%d) to (%d,%d)\n",
				 name1, state1->move_count, *count1,
				 state1->move1, state1->move2, x1, y1);
			fprintf (stderr,
				 "\tFile %s, lines %d and %d: draw from (%d,%d) to (%d,%d)\n",
				 name2, state2->move_count, *count2,
				 state2->move1, state2->move2, x2, y2);

			errorcount++;
		    }
		}
		else
		{
		    if ((abs (x1 - state2->move1) > VECTOR_TOL) || 
			(abs (y1 - state2->move2) > VECTOR_TOL)
			|| (abs (state1->move1 - x2) > VECTOR_TOL)
			|| (abs (state1->move2 - y2) > VECTOR_TOL))
		    {
			sf_warning ("Coordinate mismatch: draw");
			fprintf (stderr,
				 "\tFile %s, lines %d and %d: (%d,%d) to (%d,%d)\n",
				 name1, state1->move_count, *count1,
				 state1->move1, state1->move2, x1, y1);
			fprintf (stderr,
				 "\tFile %s, lines %d and %d: (%d,%d) to (%d,%d) (swapped)\n",
				 name2, *count2, state2->move_count,
				 x2, y2, state2->move1, state2->move2);

			errorcount++;
		    }
		}
	    }

	    state1->move1 = x1;
	    state1->move2 = y1;
	    state1->move_count = *count1;
	    state2->move1 = x2;
	    state2->move2 = y2;
	    state2->move_count = *count2;
	    break;
	case VP_PLINE:
	    (*count1)++;
	    npts1 = vp_getint2 (name1, *count1, stream1);
	    if (debug1)
		printf ("%c %d\n", c, npts1);

	    (*count2)++;
	    npts2 = vp_getint2 (name2, *count2, stream2);
	    if (debug2)
		printf ("%c %d\n", c, npts2);

	    if (npts1 != npts2)
	    {
		sf_warning ("Command mismatch: pline");
		fprintf (stderr,
			 "\tFile %s, line %d: %c %d\n", name1, *count1, c, npts1);
		fprintf (stderr,
			 "\tFile %s, line %d: %c %d\n", name2, *count2, c, npts2);

		exit (2);
	    }


	    for (ii1 = 0; ii1 < npts1; ii1++)
	    {
		(*count1)++;
		x1 = vp_getint2 (name1, *count1, stream1);
		y1 = vp_getint2 (name1, *count1, stream1);
		if (debug1)
		    printf ("%d %d\n", x1, y1);

		(*count2)++;
		x2 = vp_getint2 (name2, *count2, stream2);
		y2 = vp_getint2 (name2, *count2, stream2);
		if (debug2)
		    printf ("%d %d\n", x2, y2);

		if ((abs (x1 - x2) > VECTOR_TOL) || 
		    (abs (y1 - y2) > VECTOR_TOL))
		{
		    sf_warning ("Coordinate mismatch: pline");
		    fprintf (stderr,
			     "\tFile %s, line %d: %d %d\n", name1, *count1, x1, y1);
		    fprintf (stderr,
			     "\tFile %s, line %d: %d %d\n", name2, *count2, x2, y2);

		    errorcount++;
		}
	    }

	    if (npts1 > 0)
	    {
		state1->move1 = x1;
		state1->move2 = y1;
		state1->move_count = *count1;
		state2->move1 = x2;
		state2->move2 = y2;
		state2->move_count = *count2;
	    }
	    break;
	case VP_PMARK:
	    (*count1)++;
	    npts1 = vp_getint2 (name1, *count1, stream1);
	    mtype1 = vp_getint2 (name1, *count1, stream1);
	    msize1 = vp_getint2 (name1, *count1, stream1);
	    if (debug1)
		printf ("%c %d %d %d\n", c, npts1, mtype1, msize1);

	    (*count2)++;
	    npts2 = vp_getint2 (name2, *count2, stream2);
	    mtype2 = vp_getint2 (name2, *count2, stream2);
	    msize2 = vp_getint2 (name2, *count2, stream2);
	    if (debug2)
		printf ("%c %d %d %d\n", c, npts2, mtype2, msize2);

	    if (npts1 != npts2 || mtype1 != mtype2 || msize1 != msize2)
	    {
		sf_warning ("Command mismatch: pmark");
		fprintf (stderr,
			 "\tFile %s, line %d: %c %d %d %d\n",
			 name1, *count1, c, npts1, mtype1, msize1);
		fprintf (stderr,
			 "\tFile %s, line %d: %c %d %d %d\n",
			 name2, *count2, c, npts2, mtype2, msize2);

		exit (2);
	    }

	    while (npts1--)
	    {
		(*count1)++;
		x1 = vp_getint2 (name1, *count1, stream1);
		y1 = vp_getint2 (name1, *count1, stream1);
		if (debug1)
		    printf ("%d %d\n", x1, y1);

		(*count2)++;
		x2 = vp_getint2 (name2, *count2, stream2);
		y2 = vp_getint2 (name2, *count2, stream2);
		if (debug2)
		    printf ("%d %d\n", x2, y2);

		if ((abs (x1 - x2) > VECTOR_TOL) || 
		    (abs (y1 - y2) > VECTOR_TOL))
		{
		    sf_warning ("Coordinate mismatch: pmark");
		    fprintf (stderr,
			     "\tFile %s, line %d: %d %d\n", name1, *count1, x1, y1);
		    fprintf (stderr,
			     "\tFile %s, line %d: %d %d\n", name2, *count2, x2, y2);

		    errorcount++;
		}
	    }

	    state1->move1 = VPLOTDIFF_NOT_INITIALIZED;
	    state1->move2 = VPLOTDIFF_NOT_INITIALIZED;
	    state1->move_count = *count1;
	    state2->move1 = VPLOTDIFF_NOT_INITIALIZED;
	    state2->move2 = VPLOTDIFF_NOT_INITIALIZED;
	    state2->move_count = *count2;
	    break;
	case VP_AREA:
	    (*count1)++;
	    npts1 = vp_getint2 (name1, *count1, stream1);
	    if (debug1)
		printf ("%c %d\n", c, npts1);

	    (*count2)++;
	    npts2 = vp_getint2 (name2, *count2, stream2);
	    if (debug2)
		printf ("%c %d\n", c, npts2);

	    if (npts1 != npts2)
	    {
		sf_warning ("Command mismatch: area");
		fprintf (stderr,
			 "\tFile %s, line %d: %c %d\n", name1, *count1, c, npts1);
		fprintf (stderr,
			 "\tFile %s, line %d: %c %d\n", name2, *count2, c, npts2);

		exit (2);
	    }

	    for (i1 = 0; i1 < npts1; i1++)
	    {
		(*count1)++;
		x1 = vp_getint2 (name1, *count1, stream1);
		y1 = vp_getint2 (name1, *count1, stream1);
		if (debug1)
		    printf ("%d %d\n", x1, y1);

		(*count2)++;
		x2 = vp_getint2 (name2, *count2, stream2);
		y2 = vp_getint2 (name2, *count2, stream2);
		if (debug2)
		    printf ("%d %d\n", x2, y2);

		if ((abs (x1 - x2) > VECTOR_TOL) || 
		    (abs (y1 - y2) > VECTOR_TOL))
		{
		    sf_warning ("Coordinate mismatch: area");
		    fprintf (stderr,
			     "\tFile %s, line %d: %d %d\n", name1, *count1, x1, y1);
		    fprintf (stderr,
			     "\tFile %s, line %d: %d %d\n", name2, *count2, x2, y2);

		    errorcount++;
		}
	    }
	    state1->move1 = VPLOTDIFF_NOT_INITIALIZED;
	    state1->move2 = VPLOTDIFF_NOT_INITIALIZED;
	    state1->move_count = *count1;
	    state2->move1 = VPLOTDIFF_NOT_INITIALIZED;
	    state2->move2 = VPLOTDIFF_NOT_INITIALIZED;
	    state2->move_count = *count2;
	    break;
	case VP_OLDAREA:
	    (*count1)++;
	    npts1 = vp_getint2 (name1, *count1, stream1);
	    if (debug1)
		printf ("%c %d\n", c, npts1);

	    (*count2)++;
	    npts2 = vp_getint2 (name2, *count2, stream2);
	    if (debug2)
		printf ("%c %d\n", c, npts2);

	    if (abs (npts1 - npts2) > VECTOR_TOL)
	    {
		sf_warning ("Command mismatch: oldarea");
		fprintf (stderr,
			 "\tFile %s, line %d: %c %d\n", name1, *count1, c, npts1);
		fprintf (stderr,
			 "\tFile %s, line %d: %c %d\n", name2, *count2, c, npts2);

		exit (2);
	    }

	    (*count1)++;
	    fat1 = vp_getint2 (name1, *count1, stream1);
	    xmask1 = vp_getint2 (name1, *count1, stream1);
	    ymask1 = vp_getint2 (name1, *count1, stream1);
	    if (debug1)
		printf ("%d %d %d\n", fat1, xmask1, ymask1);

	    (*count2)++;
	    fat2 = vp_getint2 (name2, *count2, stream2);
	    xmask2 = vp_getint2 (name2, *count2, stream2);
	    ymask2 = vp_getint2 (name2, *count2, stream2);
	    if (debug2)
		printf ("%d %d %d\n", fat2, xmask2, ymask2);

	    if (fat1 != fat2 || xmask1 != xmask2 || ymask1 != ymask2)
	    {
		sf_warning ("Command mismatch: oldarea");
		fprintf (stderr,
			 "\tFile %s, line %d: %d %d %d\n",
			 name1, *count1, fat1, xmask1, ymask1);
		fprintf (stderr,
			 "\tFile %s, line %d: %d %d %d\n",
			 name2, *count2, fat2, xmask2, ymask2);

		exit (2);
	    }

	    for (i1 = 0; i1 < npts1; i1++)
	    {
		(*count1)++;
		x1 = vp_getint2 (name1, *count1, stream1);
		y1 = vp_getint2 (name1, *count1, stream1);
		if (debug1)
		    printf ("%d %d\n", x1, y1);

		(*count2)++;
		x2 = vp_getint2 (name2, *count2, stream2);
		y2 = vp_getint2 (name2, *count2, stream2);
		if (debug2)
		    printf ("%d %d\n", x2, y2);

		if ((abs (x1 - x2) > VECTOR_TOL) || 
		    (abs (y1 - y2) > VECTOR_TOL))
		{
		    sf_warning ("Coordinate mismatch: oldarea");
		    fprintf (stderr,
			     "\tFile %s, line %d: %d %d\n", name1, *count1, x1, y1);
		    fprintf (stderr,
			     "\tFile %s, line %d: %d %d\n", name2, *count2, x2, y2);

		    errorcount++;
		}
	    }
	    state1->move1 = VPLOTDIFF_NOT_INITIALIZED;
	    state1->move2 = VPLOTDIFF_NOT_INITIALIZED;
	    state1->move_count = *count1;
	    state2->move1 = VPLOTDIFF_NOT_INITIALIZED;
	    state2->move2 = VPLOTDIFF_NOT_INITIALIZED;
	    state2->move_count = *count2;
	    break;
	case VP_BACKGROUND:
	case VP_ERASE:
	    if (debug1)
		printf ("%c\n", c);
	    (*count1)++;

	    if (debug2)
		printf ("%c\n", c);
	    (*count2)++;
	    break;
	default:
	    sf_warning ("%s, line %d: unknown command %c (%o octal)",
			name1, *count1, c, c);
	    sf_error   ("%s, line %d: unknown command %c (%o octal)",
			name2, *count2, c, c);
	    break;
    }

    return errorcount;
}

