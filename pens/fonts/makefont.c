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
 *  source file:   ./filters/include/vplotfonts/makefont.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger Oct 18 1987
 *	Made arrays allocated dynamically
 * Joe Dellinger May 4 1990
 *	Don't try to mix shorts and ints, just make everything ints.
 *	Some machines have pointer alignment problems. Change the
 *	magic FONTCHECK number so that the old fonts won't blow up Vplot
 *	if they're accidentally used.
 */

/*
 * This program converts easily understandable "Vplot font" files into
 * an "include" form which is suitable for being included into gentext.c,
 * and also simultaneously into a binary ("bin") form that can be loaded at
 * runtime by gentext.c. (The Makefile in this directory shows how the
 * several resulting files must be combined.)
 *
 * makefont font_name < Vplot_font > Vplot_font.include
 *
 * If no font_name is given, it will be called "XXXXX".
 *
 * A "Vplot font" file has the following form:
 * Line 1:   start  end
 * These are the ASCII values of the beginning and ending characters
 * which will be in this font.
 * (Beginning and ending as in lowest and highest numbered in the ASCII
 * collating sequence... it is OK to not have every character in that interval)
 * Line 2: letter line space
 * Space between letters, space between lines, width of a space.
 * Line 3: top cap half base bottom
 * Vertical coordinates for "top" (Top of highest normal character in font,
 * typically), "cap" (Top of typical capital letter), "half" (What you align
 * on for vertical centering), "base" (Bottom of most letters, except those
 * with descenders), "bottom" (Bottom of lowest normal descender)
 * N Lines: Lig val1 val2 val3 ...
 * The first number is the integer ascii value of a glyph which is a ligature
 * of the ascii characters with value val1, val2, ... (up to 6)
 * A blank line signifies the end of the ligatures.
 * Start of a character line: e left right char symb
 * 'e' is simply the Vplot character for "erase", and is fixed.
 * "char" tells what ASCII code to assign this glyph to... (need be in no
 * particular order!) If "char" is outside the normal ASCII range it can also
 * be represented \num , where num is the "ASCII" number of the glyph. num
 * is NOT restricted to be less than 256! These characters can be accessed
 * via the \vnum command of gentext.
 * left, right, give the left and right half-widths
 * of the character (it is assumed that all characters are centered at
 * coordinate zero left and right), symb gives the vertical coordinate for
 * centering if this character is used as a symbol. If "left", "right",
 * or "symb" are unknown, these may be replaced instead by "5555" and
 * this program will calculate them for you. (Also leftright=0 set inside
 * this program will force such calculation of all left and right widths,
 * as will symbol=0 force such calculation of all vertical centering
 * coordinates. Instead of 5555, "symb" can also simply be omitted.)
 * Move line: m xcoordinate ycoordinate
 * Start of actual glyph coordinates... you must always start with a move.
 * Coordinates must be LESS THAN 64 in absolute value.
 * Draw line: d xcoordinate ycoordinate
 * Typically the next line will be a draw.
 * Start polygon line: A npts
 * It is also possible to put in a polygon.
 * Polygon line: tab xcoordinate ycoordinate
 * there should be following a Start Polygon line npts lines of coordinates,
 * each beginning with a tab or a space.
 * Continue in this way until the next 'e' signals the start of the next
 * character.
 *
 * Files in this format may be viewed using plas:
 * plas < Vplot_font | ?pen xcenter=0 ycenter=0 scale=150 ...
 *
 * Refer to the pen.vplot_font file given here as an example.
 *
 * - Joe Dellinger, Stanford University Dept of Earth Science
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#define FONTCHECK -1990 /* Magic number that identifies font files */

#define EOCBIT 0x8000	/* END OF CHARACTER MARK OR POLYGON BIT */
#define DRAWBIT 0x4000	/* DRAW BIT */
#define XBIT 0x0040
#define YBIT 0x2000
#define UNDEFINED	-1

/* Too big a number to be a possible width or height, which is bounded by 64 */
#define MAGIC		5555
#define MAX_LENGTH	100000	/* Max number of lines in input Vplot file */
#define MAX_FONT	32768	/* Max number of symbols in a font */
#define MAX_STRING	10	/* Max length of char number */

int             length[7] =
{
    0, 0, 0, 0, 0, 0, 0
};
int             length2[7] =
{
    0, 0, 0, 0, 0, 0, 0
};

static void write_int(int fd, int* c)
{
    if (sizeof(int) != write (fd, (char *) c, sizeof (int))) exit(2);
}

int main (int argc, char *argv[])
{
    char           *mde, *cp;
    int            *x, *y, *z, *xp, *yp, *zp;
    char           *ch, *chp;
    int            *addr;
    int             addr_index;
    int             xout=0, xmax=0, xmin=0, ymin=0, ymax=0, nc = 0;
    int            *lwidth, *rwidth;
    int            *symb;
    int             pad, leftright, symbol;
    char            string[133];
    int             start, end;
    int             first = 1;
    int             letter, line, space, top, cap, half, base, bottom;
    int             left=0, right=0, vsymb=0, count=0;
    char            name[100];
    int             i, j, iaddr;
    int             fd;
    unsigned int    uint;
    int             sint;
    int             integer;
    int             lengtht;
    int             lig[7];

    mde = (char *) malloc (MAX_LENGTH * sizeof (char));
    x = (int *) malloc (MAX_LENGTH * sizeof (int));
    y = (int *) malloc (MAX_LENGTH * sizeof (int));
    z = (int *) malloc (MAX_LENGTH * sizeof (int));
    ch = (char *) malloc (MAX_LENGTH * MAX_STRING * sizeof (char));
    addr = (int *) malloc (MAX_FONT * sizeof (int));
    lwidth = (int *) malloc (MAX_FONT * sizeof (int));
    rwidth = (int *) malloc (MAX_FONT * sizeof (int));
    symb = (int *) malloc (MAX_FONT * sizeof (int));

    for (i = 0; i < MAX_FONT; i++)
    {
	addr[i] = UNDEFINED;
	lwidth[i] = UNDEFINED;
	rwidth[i] = UNDEFINED;
    }
    for (i = 0; i < MAX_LENGTH; i++)
    {
	x[i] = MAGIC;
	y[i] = MAGIC;
	z[i] = MAGIC;
    }

    xp = x;
    yp = y;
    zp = z;
    cp = mde;
    chp = ch;

/*
 * Change these next 3 lines in order to add padding or to always ignore
 * the given width and height information
 */
    pad = 0;
    leftright = 1;
    symbol = 1;


    if (argc == 1)
    {
	strcpy (name, "XXXXX");
    }
    else
    {
	strcpy (name, argv[1]);
    }
    strcat (name, "_");

    if (NULL == fgets (string,132,stdin)) exit(1);
    sscanf (string, "%d %d", &start, &end);
    if (NULL == fgets (string,132,stdin)) exit(1);
    sscanf (string, "%d %d %d", &letter, &line, &space);
    space -= 2 * letter;
    if (NULL == fgets (string,132,stdin)) exit(1);
    sscanf (string, "%d %d %d %d %d", &top, &cap, &half, &base, &bottom);

    sprintf (string, "%slig", name);
    fd = creat (string, 0777);
    printf ("int %slig[] =\n{\n", name);

    while (1)
    {
	if (NULL == fgets (string,132,stdin)) exit(1);
	/* At most 6 characters in a ligature! */
	lig[0] = 0;
	i = sscanf (string, "%d %d %d %d %d %d %d",
		    lig, lig + 1, lig + 2, lig + 3, lig + 4, lig + 5, lig + 6);
	if (i <= 1)
	{
	    printf ("%d,  ", 0);
	    integer = 0;
	    write_int (fd, &integer);
	    length[6] += sizeof (int);

	    printf ("\n};\n\n");
	    close (fd);
	    break;
	}
	else
	{
	    printf ("%d,  ", i - 1);
	    integer = i - 1;
	    write_int (fd, &integer);
	    length[6] += sizeof (int);

	    printf ("%d, ", lig[0]);
	    integer = lig[0];
	    write_int (fd, &integer);
	    length[6] += sizeof (int);

	    for (j = 1; j < i; j++)
	    {
		printf ("%d,", lig[j]);
		integer = lig[j];
		write_int (fd, &integer);
		length[6] += sizeof (int);
	    }
	    printf ("\n");
	}
    }

    while (fgets (string,132,stdin))
    {
	*chp = '\0';
	*(chp + 1) = '\0';
	*(chp + 2) = '\0';
	sscanf (string, "%c %d %d %s %d", cp++, xp++, yp++, chp, zp++);
	chp += MAX_STRING;
	nc++;
    }
    *cp = 'e';
    *chp = '\0';
    nc++;

    printf ("unsigned int %svec[] =\n{\n", name);
    sprintf (string, "%ssvec", name);
    fd = creat (string, 0777);

    for (i = 0, iaddr = 0; i < nc; ++i)
    {
	switch (mde[i])
	{
	    case ('e'):
		/* Things to do when ending a character */
		if (!first)
		{
		    xout = EOCBIT;
		    if (!leftright || left == MAGIC || right == MAGIC)
		    {
			lwidth[addr_index] = -xmin + pad / 2.;
			rwidth[addr_index] = xmax + pad / 2.;
		    }
		    else
		    {
			lwidth[addr_index] = left;
			rwidth[addr_index] = right;
		    }
		    if (!symbol || vsymb == MAGIC)
		    {
			symb[addr_index] = (ymax + ymin) / 2;
		    }
		    else
		    {
			symb[addr_index] = vsymb;
		    }
		}
		else
		    first = 0;

		/* Things to do when beginning a character */
		if (ch[i * MAX_STRING + 0] != '\0')
		{
		    count = 0;
		    ymax = -10000;
		    ymin = 10000;
		    xmax = -10000;
		    xmin = 10000;
		    if (ch[i * MAX_STRING + 0] != '\\' || ch[i * MAX_STRING + 1] == '\0')
		    {
			addr_index = ch[i * MAX_STRING + 0];
		    }
		    else
		    {
			sscanf (&ch[i * MAX_STRING + 0], "\\%d ", &addr_index);
		    }
		    if (addr_index < start || addr_index > end)
		    {
			fprintf (stderr, "Character value %d out of promised range %d to %d!\n", addr_index, start, end);
			exit (1);
		    }
		    addr_index -= start;
		    addr[addr_index] = iaddr;
		    left = x[i];
		    right = y[i];
		    vsymb = z[i];
		}
		break;
	    case 'A':
		count = x[i];
		if (count < 3)
		{
		    fprintf (stderr, "Polygon must have more than 2 vertices!\n");
		    exit (2);
		}
		/* Don't increment iaddr, no output */
		continue;
		break;
	    case '\t':
		mde[i] = ' ';
	    case ' ':
		count--;
	    case ('d'):
	    case ('m'):
		xout = x[i] < 0 ? -(x[i]) : x[i];
		if (x[i] < 0)
		    xout |= XBIT;
		xout |= ((y[i] < 0 ? -(y[i]) : y[i]) << 7);
		if (y[i] < 0)
		    xout |= YBIT;
		if (mde[i] == 'd')
		    xout |= DRAWBIT;
		else
		    if (mde[i] == ' ')
		    {
			xout |= DRAWBIT;
			if (count > 0)
			    xout |= EOCBIT;
		    }

		xmax = x[i] < xmax ? xmax : x[i];
		xmin = x[i] > xmin ? xmin : x[i];
		ymax = y[i] < ymax ? ymax : y[i];
		ymin = y[i] > ymin ? ymin : y[i];
		break;
	    default:
		/* ignore this line */
		continue;
		break;
	}

	/*
	 * Only increment pointer for those things which output something,
	 * and also the first time. 
	 */
	iaddr++;
	if (i > 0)
	{
	    printf ("%d,", xout);
	    uint = xout;
	    if (sizeof(unsigned int) != 
		write (fd, (char *) &uint, sizeof (unsigned int))) exit(2);
	    length[5] += sizeof (unsigned int);
	    if (xout == EOCBIT)
		printf ("\n");
	}
    }				/* end for loop  */
    close (fd);

/*
 * Check info
 */
    sprintf (string, "%scheck", name);
    fd = creat (string, 0777);
    if (20 != write (fd, (char *) "Vplot Binary fonT  \n", 20)) exit(2);
    integer = FONTCHECK;
    write_int (fd, &integer);
    close (fd);

/*
 * Now write out all the Binary information.
 * We will write out each different sort of thing into a different
 * file, and then use "cat" to put them all together.
 * Keep track of how long each piece is as we write it out
 * so that we construct the "table of contents" to the concatenated file.
 */

    sprintf (string, "%saddr", name);
    fd = creat (string, 0777);

    printf ("};\n\n\nint %saddr[] =\n{\n", name);
    for (i = 0; i <= end - start; ++i)
    {
	printf ("%d,", addr[i]);
	integer = addr[i];
	write_int (fd, &integer);
	length[1] += sizeof (int);
    }
    close (fd);

    sprintf (string, "%swidthl", name);
    fd = creat (string, 0777);

    printf ("\n};\n\n\nint %swidthl[] =\n{\n", name);
    for (i = 0; i <= end - start; ++i)
    {
	printf ("%d,", lwidth[i]);
	sint = lwidth[i];
	write_int (fd, &sint);
	length[2] += sizeof (int);
    }
    close (fd);

    sprintf (string, "%swidthr", name);
    fd = creat (string, 0777);

    printf ("\n};\n\n\nint %swidthr[] =\n{\n", name);
    for (i = 0; i <= end - start; ++i)
    {
	printf ("%d,", rwidth[i]);
	sint = rwidth[i];
	write_int (fd, &sint);
	length[3] += sizeof (int);
    }
    close (fd);

    sprintf (string, "%ssymbol", name);
    fd = creat (string, 0777);

    printf ("\n};\n\n\nint %ssymbol[] =\n{\n", name);
    for (i = 0; i <= end - start; ++i)
    {
	printf ("%d,", symb[i]);
	sint = symb[i];
	write_int (fd, &sint);
	length[4] += sizeof (int);
    }
    close (fd);

    printf ("\n};\n\n\n");

    printf ("int %sdim[] = {%d, %d, %d, %d, %d, %d, %d, %d, %d, %d};\n",
	    name, bottom, base, half, cap, top, letter, line, space, start, end);
    sprintf (string, "%sdim", name);
    fd = creat (string, 0777);
    sint = bottom;
    write_int (fd, &sint);
    sint = base;
    write_int (fd, &sint);
    sint = half;
    write_int (fd, &sint);
    sint = cap;
    write_int (fd, &sint);
    sint = top;
    write_int (fd, &sint);
    sint = letter;
    write_int (fd, &sint);
    sint = line;
    write_int (fd, &sint);
    sint = space;
    write_int (fd, &sint);
    sint = start;
    write_int (fd, &sint);
    sint = end;
    write_int (fd, &sint);
    length[0] += 10 * sizeof (int);
    close (fd);

/*
 * Length is an array which says how long each piece is;
 * Length2 says at what offset each begins.
 * Lengtht gives the total length of the whole file (not including the
 * table of contents itself, which goes at the beginning and is of
 * fixed length).
 */

    lengtht = length[0];
    for (i = 1; i < 7; i++)
    {
	length2[i] = length2[i - 1] + length[i - 1];
	lengtht += length[i];
    }

    sprintf (string, "%sheader", name);
    fd = creat (string, 0777);
    write_int (fd, &lengtht);
    if (7 * sizeof(int) != write (fd, (char *) length2, 7 * sizeof (int)))
	exit(3);
    close (fd);
    exit(0);
}
