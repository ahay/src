/* Plot Debugger - convert vplot to ascii.

Takes: < plot.vpl > plot.asc
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <rsfplot.h>

static void text (void);

int main (int argc, char* argv[])
{
    int c, a, ix, iy, iz, npts, mtype, key, col_tab_no, byte, i;
    int ii, num_byte, nmul, nx, ny, xmask, ymask, orient;
    int ras_orient, ras_offset, xpix, ypix, num_rep, num_pat, pos;
    int ibyte=0, ipat, col, j, bit;
    float x, y, xcor, ycor, red, green, blue, fat, msize, size;
    float xmin, ymin, xmax, ymax, xvplot, yvplot, off, rep;
    float fatscale, txvecscale, txscale, scale, colscale, hscale;

    fatscale = 1.;
    txvecscale = 1.;
    txscale = 1.;
    scale = 1.;
    colscale = 1.;
    hscale = 1.;

    printf ("#plas: Vplot units used in this file\n");

    while (EOF != (c = getchar ())) {
	switch (c) {
	    case VP_MESSAGE:
		printf ("%c\n", c);
		text ();
	    break;
	    case VP_SETSTYLE:
		a = getchar ();
		printf ("%c %c\n", c, a);
		break;
	    case VP_ERASE:
	    case VP_BREAK:
	    case VP_PURGE:
	    case VP_NOOP:
		printf ("%c\n", c);
		break;
	    case VP_ORIGIN:
	    case VP_MOVE:
	    case VP_DRAW:
		x = scale * vp_getint();
		y = scale * vp_getint();
		printf ("%c %g %g\n", c, x, y);
		break;
	    case VP_TXALIGN:
		ix = vp_getint();
		iy = vp_getint();
	    printf ("%c %d %d\n", c, ix, iy);
	    break;
	    case VP_TXFONTPREC:
		ix = vp_getint();
		iy = vp_getint();
		iz = vp_getint();
		printf ("%c %d %d %d\n", c, ix, iy, iz);
		break;
	    case VP_PLINE:
	    case VP_SETDASH:
		npts = vp_getint();
		printf ("%c %d\n", c, npts);
		while (npts--) {
		    x = scale * vp_getint();
		    y = scale * vp_getint();
		    printf ("%g %g\n", x, y);
		}
		break;
	    case VP_PMARK:
		npts = vp_getint();
		mtype = vp_getint();
		msize = txscale * vp_getint();
		printf ("%c %d %d %g\n", c, npts, mtype, msize);
		while (npts--) {
		    x = scale * vp_getint();
		    y = scale * vp_getint();
		    printf ("%g %g\n", x, y);
		}
		break;
	    case VP_BEGIN_GROUP:
		printf ("%c\n", c);
		text ();
		break;
	    case VP_END_GROUP:
		printf ("%c\n", c);
		break;
	    case VP_GTEXT:
		x = scale * txvecscale * vp_getint();
		y = scale * txvecscale * vp_getint();
		xcor = scale * txvecscale * vp_getint();
		ycor = scale * txvecscale * vp_getint();
		printf ("%c %g %g %g %g\n", c, x, y, xcor, ycor);
		text ();
		break;
	    case VP_OLDTEXT:
		key = vp_getint();
		size = txscale * (key & 037);
		orient = (key & 0140) >> 5;
		printf ("%c %g %d\n", VP_TEXT, size, 90 * orient);
		text ();
		break;
	    case VP_TEXT:
		size = txscale * vp_getint();
		orient = vp_getint();
		printf ("%c %g %d\n", c, size, orient);
		text ();
		break;
	    case VP_OVERLAY:
	    case VP_COLOR:
		ix = vp_getint();
		printf ("%c %d\n", c, ix);
		break;
	    case VP_SET_COLOR_TABLE:
		col_tab_no = vp_getint();
		red = colscale * vp_getint();
		green = colscale * vp_getint();
		blue = colscale * vp_getint();
		printf ("%c %d %g %g %g\n", c, col_tab_no, red, green, blue);
		break;
	    case VP_FAT:
		fat = fatscale * vp_getint();
		printf ("%c %g\n", c, fat);
		break;
	    case VP_WINDOW:
		xmin = scale * vp_getint();
		ymin = scale * vp_getint();
		xmax = scale * vp_getint();
		ymax = scale * vp_getint();
		printf ("%c %g %g %g %g\n", c, xmin, ymin, xmax, ymax);
		break;
	    case VP_BIT_RASTER:
	    case VP_BYTE_RASTER:
		ras_orient = vp_getint();
		ras_offset = vp_getint();
		printf ("%c %d %d\n", c, ras_orient, ras_offset);
		xcor = scale * vp_getint();
		ycor = scale * vp_getint();
		printf ("%g %g\n", xcor, ycor);
		xvplot = scale * vp_getint();
		yvplot = scale * vp_getint();
		printf ("%g %g\n", xvplot, yvplot);
		xpix = vp_getint();
		ypix = vp_getint();
		printf ("%d %d\n", xpix, ypix);
		for (i = 0; i < ypix; i += num_rep) {
		    num_rep = vp_getint();
		    printf ("     %d  ------ lines %d to %d\n",
			    num_rep, i, i + num_rep - 1);
		    for (pos = 0; pos < xpix; pos += num_byte * num_pat) {
			num_pat = vp_getint();
			num_byte = vp_getint();
			printf ("%d %d  -- start byte %d in y_line %d\n",
				num_pat, num_byte, pos, i);
			if (num_pat < 0 || num_byte < 0)
			{
			    fprintf (stderr, "Error in raster field\n");
			    exit (1);
			}

			if (VP_BYTE_RASTER == c) {
			    for (ii = 0; ii < num_byte; ii++) {
				byte = getchar ();
				if (EOF == byte) {
				    fprintf (stderr, "Error in raster field\n");
				    exit (1);
				}
				printf ("%d\n", byte);
			    }
			} else { /* VP_BIT_RASTER */
			    for (ii = 0; ii < num_byte; ii++) {
				if (0 == ii % 8) {
				    ibyte = getchar ();
				} else {
				    ibyte <<= 1;
				}
				printf ("%d\n", ((ibyte >> 7) & 001));
			    }
			}
		    }
		}
		break;
	    case VP_AREA:
		npts = vp_getint();
		printf ("%c %d\n", c, npts);
		for (i = 0; i < npts; i++) {
		    x = scale * vp_getint();
		    y = scale * vp_getint();
		    printf ("%g %g\n", x, y);
		}
		break;
	    case VP_PATLOAD:
		nmul = vp_getint();
		nx = vp_getint();
		ny = vp_getint();
		ipat = vp_getint();
		printf ("%c %d %d %d %d\n",
			c, nmul, nx, ny, ipat);
		if (-1 == nx) {
		    for (i = 0; i < ny * 2; i++) {
			fat = fatscale * vp_getint();
			col = vp_getint();
			off = hscale * vp_getint();
			rep = hscale * vp_getint();
			printf ("%g %d %g %g\n",
				fat, col, off, rep);
		    }
		} else {
		    for (i = 0; i < nx; i++) {
			for (j = 0; j < ny; j++) {
			    bit = vp_getint();
			    if (bit >= 0 && bit <= 9) {
				putchar (bit + '0');
			    } else {
				putchar (bit + 'A' - 10);
			    }
			}
			putchar ('\n');
		    }
		}
		break;
	    case VP_OLDAREA:
		npts = vp_getint();
		printf ("%c %d\n", c, npts);
		fat = fatscale * vp_getint();
		xmask = vp_getint();
		ymask = vp_getint();
		printf ("%g %d %d\n", fat, xmask, ymask);
		for (i = 0; i < npts; i++)
		{
		    x = scale * vp_getint();
		    y = scale * vp_getint();
		    printf ("%g %g\n", x, y);
		}
		break;
	    default:
		fprintf (stderr,"******unknown command %c %o\n", c, c);
		exit(0);
		break;
	}
    }
    
    exit (0);
}

static void text (void)
{
    int c;

    while ('\0' != (c = getchar ())) {
	putchar (c);
    }
    putchar ('\n');
}

/* 	$Id: pldb.c,v 1.5 2003/10/01 23:41:18 fomels Exp $	 */

