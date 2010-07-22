/* Plot Assembler - convert ascii to vplot. */
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
/* Modified from the original version by Joe Dellinger */

#include <stdio.h>
#include <math.h>
#include <string.h>

#include <rsfplot.h>

#define MAXLINE 2016

#define GETLINE(a) if (NULL== fgets (a, MAXLINE, stdin)) \
	sf_error("fgets error")

static void text (void);
static void setunits(char c);

static float fatscale, txvecscale, txscale, scale, colscale, hscale;
static char  *documentation[] = {
    "",
    "",
    "NAME",
    "		sfplas 	(PLot ASsembler) -- deformat vector plot commands",
    "SYNOPSIS",
    "		sfplas < plot.asc > plot.vpl",
    "------------------------------------------------------------------------",
    "	Reads stdin; writes stdout.  Complement to 'pldb'.",
    "	Converts vplot-language files in ascii (human-readable) form",
    "	into the hybrid ascii/binary form that the 'pen' filters understand.",
    "",
    "   A backslash at the end of a line of text allows multi-line strings.",
    "Options",
    "      -v: Use vplot units (the default)",
    "      -i: Use inches",
    "      -c: Use centimeters",
    "  Note that these options apply to ALL geometric attributes,",
    "  including line width and character height scale factors.",
    "",
    "The default may be set within the file by having as the first line",
    "#plas: X",
    "where X is one of V, C, or I for Vplot, Cm, or Inches",
    "SEE ALSO",
    "	sfpldb."
};

int main (int argc, char* argv[])
{
    char units='v';
    int npat, line_count, ix, iy, iz, npts, mtype, orient, maskx, masky, ic;
    int i, j, nmul, ipat, col, ras_orient, ras_offset, col_tab_no, nx, ny;
    int xpix, ypix, num_rep, count, num_pat, num_byte, byte, ibyte=0, ndoc;
    float x, y, xcor, ycor, msize, size, fat, off, rep, red, green, blue;
    float xmin, ymin, xmax, ymax, xvplot, yvplot;
    char line[MAXLINE], c, a, *ptr;

    if (isatty(fileno(stdin))) { /* no input - do selfdoc */
	ndoc = sizeof(documentation)/sizeof(documentation[0]);
	for (i = 0; i < ndoc; i++) {
	    printf ("%s\n", documentation[i]);
	}
	exit(1);
    }

    vp_init();

    for (ic=1; ic < argc; ic++) {
	if ('-' == argv[ic][0]) units=argv[ic][1];
    }
    setunits (units);

    npat=0;
    for (line_count=0; fgets (line, MAXLINE, stdin) != NULL; line_count++) {
	c = line[0];

	switch (c) {
	    case '#':
		if (0 == line_count &&
		    0 == strncmp (line, "#plas: ", 7)) setunits (line[7]);
		break;
	    case VP_MESSAGE:
		putchar (c);
		text ();
		break;
	    case VP_SETSTYLE:
		sscanf (line, "%*c %c", &a);
		putchar (c);
		putchar (a);
		break;
	    case VP_ERASE:
	    case VP_PURGE:
	    case VP_BREAK:
	    case VP_NOOP:
	    case VP_BACKGROUND:
		putchar (c);
		break;
	    case VP_TXALIGN:
		putchar (c);
		sscanf (line, "%*c %d %d", &ix, &iy);
		vp_putint (ix);
		vp_putint (iy);
		break;
	    case VP_ORIGIN:
	    case VP_MOVE:
	    case VP_DRAW:
		putchar (c);
		sscanf (line, "%*c %f %f", &x, &y);
		vp_putfloat0 (x * scale);
		vp_putfloat0 (y * scale);
		break;
	    case VP_TXFONTPREC:
		putchar (c);
		sscanf (line, "%*c %d %d %d", &ix, &iy, &iz);
		vp_putint (ix);
		vp_putint (iy);
		vp_putint (iz);
		break;
	    case VP_PLINE:
	    case VP_SETDASH:
		putchar (c);
		sscanf (line, "%*c %d", &npts);
		vp_putint (npts);
		while (npts--) {
		    GETLINE(line);
		    sscanf (line, "%f %f", &x, &y);
		    vp_putfloat0 (x * scale);
		    vp_putfloat0 (y * scale);
		}
		break;
	    case VP_PMARK:
		putchar (c);
		sscanf (line, "%*c %d %d %f", &npts, &mtype, &msize);
		vp_putint (npts);
		vp_putint (mtype);
		vp_putfloat0 (msize * txscale);
		while (npts--) {
		    GETLINE(line);
		    sscanf (line, "%f %f", &x, &y);
		    vp_putfloat0 (x * scale);
		    vp_putfloat0 (y * scale);
		}
		break;
	    case VP_BEGIN_GROUP:
		putchar (c);
		text ();
		break;
	    case VP_END_GROUP:
		putchar (c);
		break;
	    case VP_OLDTEXT:
	    case VP_TEXT:
		putchar (VP_TEXT);
		sscanf (line, "%*c %f %d", &size, &orient);
		vp_putfloat0 (size * txscale);
		if (VP_OLDTEXT == c) orient *= 90;
		vp_putint (orient);
		text ();
		break;
	    case VP_GTEXT:
		putchar (c);
		sscanf (line, "%*c %f %f %f %f", &x, &y, &xcor, &ycor);
		vp_putfloat0 (x * scale * txvecscale);
		vp_putfloat0 (y * scale * txvecscale);
		vp_putfloat0 (xcor * scale * txvecscale);
		vp_putfloat0 (ycor * scale * txvecscale);
		text ();
		break;
	    case VP_COLOR:
	    case VP_OVERLAY:
		putchar (c);
		sscanf (line, "%*c %d", &ix);
		vp_putint (ix);
		break;
	    case VP_SET_COLOR_TABLE:
		putchar (c);
		sscanf (line, "%*c %d %f %f %f",
			&col_tab_no, &red, &green, &blue);
		vp_putint (col_tab_no);
		vp_putfloat0 (red * colscale);
		vp_putfloat0 (green * colscale);
		vp_putfloat0 (blue * colscale);
		break;
	    case VP_FAT:
		putchar (c);
		sscanf (line, "%*c %f", &fat);
		vp_putfloat0 (fat * fatscale);
		break;
	    case VP_WINDOW:
		putchar (c);
		sscanf (line, "%*c %f %f %f %f", &xmin, &ymin, &xmax, &ymax);
		vp_putfloat0 (scale * xmin);
		vp_putfloat0 (scale * ymin);
		vp_putfloat0 (scale * xmax);
		vp_putfloat0 (scale * ymax);
		break;
	    case VP_OLDAREA:
		putchar (c);
		sscanf (line, "%*c %d", &npts);
		GETLINE(line);
		sscanf (line, "%f %d %d", &fat, &maskx, &masky);
		vp_putint (npts);
		vp_putfloat0 (fat * fatscale);
		vp_putint (maskx);
		vp_putint (masky);
		for (i = 0; i < npts; i++) {
		    GETLINE(line);
		    sscanf (line, "%f %f", &x, &y);
		    vp_putfloat0 (x * scale);
		    vp_putfloat0 (y * scale);
		}
		break;
	    case VP_AREA:
		putchar (c);
		sscanf (line, "%*c %d", &npts);
		vp_putint (npts);
		for (i = 0; i < npts; i++) {
		    GETLINE(line);
		    sscanf (line, "%f %f", &x, &y);
		    vp_putfloat0 (scale * x);
		    vp_putfloat0 (scale * y);
		}
		break;
	    case VP_PATLOAD:
		npat++;
		ipat = 0;
		sscanf (line, "%*c %d %d %d %d", &nmul, &nx, &ny, &ipat);
		putchar (c);
		vp_putint (nmul);
		vp_putint (nx);
		vp_putint (ny);
		if (0 == ipat) ipat = npat;
		vp_putint (ipat);
		if (-1 != nx) {
		    for (i = 0; i < nx; i++) {
			GETLINE(line);
			for (j = 0, ptr = line; j < ny; j++, ptr++) {
			    if ('\0' == *ptr || '\n' == *ptr)
				fprintf (stderr, "null/nl");
			    if (*ptr >= '0' && *ptr <= '9') {
				ipat = (*ptr) - '0';
			    } else {
				ipat = (*ptr) - 'A' + 10;
			    }
			    vp_putint (ipat);
			}
		    }
		} else {
		    for (i = 0; i < ny * 2; i++) {
			GETLINE(line);
			sscanf (line, "%f %d %f %f", &fat, &col, &off, &rep);
			vp_putfloat0 (fat * fatscale);
			vp_putint (col);
			vp_putfloat0 (off * hscale);
			vp_putfloat0 (rep * hscale);
		    }
		} break;
	    case VP_BYTE_RASTER:
	    case VP_BIT_RASTER:
		sscanf (line, "%*c %d %d", &ras_orient, &ras_offset);
		putchar (c);
		vp_putint (ras_orient);
		vp_putint (ras_offset);
		GETLINE(line);
		sscanf (line, "%f %f", &xcor, &ycor);
		vp_putfloat0 (xcor * scale);
		vp_putfloat0 (ycor * scale);
		GETLINE(line);
		sscanf (line, "%f %f", &xvplot, &yvplot);
		vp_putfloat0 (scale * xvplot);
		vp_putfloat0 (scale * yvplot);
		GETLINE(line);
		sscanf (line, "%d %d", &xpix, &ypix);
		vp_putint (xpix);
		vp_putint (ypix);

		for (j = 0; j < ypix; j += num_rep) {
		    GETLINE(line);
		    sscanf (line, "%d", &num_rep);
		    if (num_rep < 1)
			fprintf (stderr, "Bad Raster repetition factor\n");
		    vp_putint (num_rep);
		    for (count=0; count < xpix; count += num_byte * num_pat) {
			GETLINE(line);
			sscanf (line, "%d %d", &num_pat, &num_byte);
			vp_putint (num_pat);
			vp_putint (num_byte);
			for (i = 0; i < num_byte; i++) {
			    GETLINE(line);
			    sscanf (line, "%d", &byte);
			    if (VP_BYTE_RASTER == c) {
				putchar ((char) byte);
			    } else {
				if (i % 8 == 0)
				    ibyte = 0;
				ibyte |= ((byte != 0) << (7 - (i % 8)));
				if (i % 8 == 7 || i == num_byte - 1)
				    putchar ((char) ibyte);
			    }
			}
		    }
		}
		break;
	    default:
		break;  /* Treat unknown characters as comments */
	}
    }

    return 0;
}

static void text (void)
{
    char c;

    while ('\n' != (c = getchar ())) {
	if ('\\' != c) {
	    putchar (c);
	} else {
	    c = getchar ();
	    if ('\n' == c) continue;

	    putchar ('\\');
	    putchar (c);
	}
    }
    putchar ('\0');
}

static void setunits(char c) {
    switch (c) {
	case 'I':
	case 'i':
	    fatscale = FATPERIN;
	    txvecscale = TEXTVECSCALE;
	    txscale = TXPERIN;
	    scale = RPERIN;
	    colscale = MAX_GUN;
	    hscale = RPERIN;
	    break;
	case 'C':
	case 'c':
	    fatscale = FATPERIN / 2.54;
	    txvecscale = TEXTVECSCALE;
	    txscale = TXPERIN / 2.54;
	    scale = RPERIN / 2.54;
	    colscale = MAX_GUN;
	    hscale = RPERIN / 2.54;
	    break;
	case 'V':
	case 'v':
	default:
	    fatscale = 1.;
	    txvecscale = 1.;
	    txscale = 1.;
	    scale = 1.;
	    colscale = 1.;
	    hscale = 1;
	    break;
    }
}
