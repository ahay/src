#include <stdio.h>

#include <rsf.h>

#include "device.h"

static int brake = BREAK_BREAK;
static bool first_time = true;
static bool mono = false;

static void init_colors(void)
{
    /* Set up the default color table */
    
    for (ii = 0; ii < 8; ii++) {
	color_set[ii][STATUS] = SET;
	color_set[ii][_RED] = MAX_GUN * ((ii & 2) / 2);
	color_set[ii][_GREEN] = MAX_GUN * ((ii & 4) / 4);
	color_set[ii][_BLUE] = MAX_GUN * ((ii & 1) / 1);
	color_set[ii][_GREY] = 
	    greycorr ((int) ((MAX_GUN * ii) / 7));
    }

    if (mono) num_col = 0;
    
    num_col_8 = (num_col > 8) ? num_col : 8;
    
    if (mono) {
	color_set[0][MAP] = 0;
	for (ii = 1; ii <= MAX_COL; ii++) {
	    color_set[ii][MAP] = 7;
	}
	for (ii = num_col_8; ii < 256; ii++) {
	    color_set[ii][_GREY] = 
		color_set[((ii - 8) % 7) + 1][_GREY];
	}
	
	/* Put a grey scale in the upper half 
	   of the color table */
	for (ii = 256; ii <= MAX_COL; ii++) {
	    color_set[ii][_GREY] = greycorr (ii - 256);
	}
    } else { /* not mono */
	for (ii = 0; ii < num_col_8; ii++) {
	    color_set[ii][MAP] = ii;
	}
	for (ii = num_col_8; ii <= MAX_COL; ii++) {
	    color_set[ii][MAP] = ((ii - 8) % 7) + 1;
	}
    }

    setstyle (new_style);
}

void vp_dovplot (void)
{
    int c, ii;

    while (((c = getchar ()) == VP_ERASE) || 
	   (c==VP_SET_COLOR_TABLE) || 
	   ((c == VP_BREAK) && (brake == BREAK_ERASE)))
    {
	switch (c) {
	    case VP_SET_COLOR_TABLE:	/* set color table entry */
		if(first_time) init_colors();
		set_table() ;
		break;
	    case VP_ERASE:
	    case VP_BREAK:
		starterase = 1;
		break;
	}
    }
}

#if 0 

    /*
     * Files that consist only of erase characters and nothing else are
     * ignored.
     */
    if (c == EOF)
	return;
    ungetc ((char) c, pltin);

    if (first_time)
    {
	dev.reset ();
/*
 * Device is now officially open for orders.
 */
	if (dev_xmax <= dev_xmin ||
	    dev_ymax <= dev_ymin ||
	    pixels_per_inch == 0. ||
	    aspect_ratio == 0. ||
	    num_col == -1)
	    ERR (FATAL, name, "Critical variables left unset by device!");

/*
 * Set maximum clipping window for dev.attributes(SET_WINDOW)
 */
	xwmax_last = dev_xmax;
	xwmin_last = dev_xmin;
	ywmax_last = dev_ymax;
	ywmin_last = dev_ymin;

/*
 * Set up color maps
 */
	init_colors ();

	ever_called = 1;
    }
    else
    {
	dev.close (CLOSE_FLUSH);
	if (epause > 0)
	{
	    sleep ((unsigned) epause);
	}
	else if (epause < 0)
	{
	    dev.close (CLOSE_FLUSH);
	    message (MESG_ERASE);
	    message (MESG_ON);
	    message (MESG_HOME);
	    message (MESG_READY);
	    message (MESG_HIGHLIGHT_ON);
	    message (MESG_TEXT, "Type Return to Continue...  ");
	    message (MESG_DONE);
/*	    dev.interact (INT_PAUSE, controltty, string);*/
	    ii = dev.interact (INT_PAUSE, controltty, string);
	    if (ii == DOVPLOT_EXIT)
	    {
		dev.close (CLOSE_FLUSH);
		return;
	    }
	    message (MESG_HIGHLIGHT_OFF);
	    if (!allowecho)
	    {
		message (MESG_READY);
		message (MESG_TEXT, CRLF);
	    }
	    message (MESG_DONE);
	    message (MESG_OFF);
	    message (MESG_ERASE);
	}
	/*
	 * Inquire point back from device
	 */
	if (interact[0] != '\0')
	{
	    getapoint ();
	}
    }

    /*
     * Erase the screen to background color
     */
    if (((erase & FORCE_INITIAL) && (first_time || (erase & DO_LITERALS)))
	|| (starterase && (erase & DO_LITERALS)))
    {

	if (first_time){
	    dev.erase (ERASE_START);
}
	else
	{
	    dev.erase (ERASE_MIDDLE);
	    nplots++;
	}
	dev.close (CLOSE_FLUSH);
    }

    first_time = NO;

/*
 * Reset fatness, cur_color, etc.
 */
    new_style = default_style;
    setstyle (new_style);
    reset ();

/*
 * Make SURE the color is what it's supposed to be, just to be safe.
 */
    dev.attributes (SET_COLOR, cur_color, 0, 0, 0);
    need_devcolor = NO;

    message (MESG_OFF);
    dev.close (CLOSE_FLUSH);

/*
 * Start a group that will contain this frame (always group 0).
 */
    ii = ftell (pltin);
    sprintf (group_name, "%s.%d", pltname, nplots);
    dev.attributes (BEGIN_GROUP, group_number, ii, 0, 0);
    group_number++;

    if (framewindows)
	outline_window ();

/*
 * Finally, here's the main loop that does the actual processing of vplot.
 * Notice that some tricky commands (VP_DRAW, VP_SET_COLOR_TABLE) look ahead
 * at the next command to try to anticipate often-occuring situations,
 * so things get a little complicated here.
 */
    while ((c = getc (pltin)) != EOF)
    {
	switch (c)		/* command list */
	{
	case VP_SETSTYLE:	/* set the style */
/*    fprintf(stderr,"in dovplot 2\n");*/
	    c = getc (pltin);
	    if ((c == 'r') || (c == 'R') || (c == 'm') || (c == 'M'))
		new_style = ROTATED;
	    else if ((c == 'o') || (c == 'O'))
		new_style = OLD;
	    else if ((c == 'a') || (c == 'A'))
		new_style = ABSOLUTE;
	    else
		new_style = STANDARD;
	    setstyle (new_style);

	    if (framewindows)
		outline_window ();

/*    fprintf(stderr,"in dovplot 3\n");*/
	    break;
	case VP_MOVE:		/* move */
/*    fprintf(stderr,"in dovplot 4\n");*/
	    /*
	     * Reset position in dash pattern.
	     */
	    dashpos = 0.;

	    GETXY (xold, yold);
/*    fprintf(stderr,"in dovplot 5\n");*/
	    break;
	case VP_DRAW:		/* draw */
/*    fprintf(stderr,"in dovplot 6\n");*/
	    GETXY (xnew, ynew);
	    update_color ();
	    while ((c = getc (pltin)) == VP_DRAW)
	    {
		GETXY (xnewer, ynewer);
		/*
		 * Is it the same point?
		 */
		if (xnewer == xnew && ynewer == ynew)
		    continue;
		/*
		 * Is it colinear and horizontal or vertical?
		 */
		if ((ynewer == ynew && ynew == yold &&
		     ((xnewer > xnew) == (xnew > xold))) ||
		    (xnewer == xnew && xnew == xold &&
		     ((ynewer > ynew) == (ynew > yold))))
		{
		    ynew = ynewer;
		    xnew = xnewer;
		    continue;
		}
		dev.vector (xold, yold, xnew, ynew, fat, dashon);
		xold = xnew;
		yold = ynew;
		xnew = xnewer;
		ynew = ynewer;
	    }
	    dev.vector (xold, yold, xnew, ynew, fat, dashon);
	    xold = xnew;
	    yold = ynew;
	    if (c == EOF)
		goto End_of_file;
	    ungetc ((char) c, pltin);
/*    fprintf(stderr,"in dovplot 7\n");*/
	    break;
	case VP_PLINE:		/* polyline */
/*    fprintf(stderr,"in dovplot 8\n");*/
	    /*
	     * Reset position in dash pattern.
	     */
	    dashpos = 0.;

	    npts = geth (pltin);
	    if (npts == 0)
		break;
	    GETXY (xold, yold);
	    npts--;
	    if (npts == 0)
		break;
	    GETXY (xnew, ynew);
	    npts--;
	    update_color ();
	    while (npts > 0)
	    {
		GETXY (xnewer, ynewer);
		npts--;
		/*
		 * Is it the same point?
		 */
		if (xnewer == xnew && ynewer == ynew)
		    continue;
		/*
		 * Is it colinear and horizontal or vertical?
		 */
		if ((ynewer == ynew && ynew == yold &&
		     ((xnewer > xnew) == (xnew > xold))) ||
		    (xnewer == xnew && xnew == xold &&
		     ((ynewer > ynew) == (ynew > yold))))
		{
		    ynew = ynewer;
		    xnew = xnewer;
		    continue;
		}
		dev.vector (xold, yold, xnew, ynew, fat, dashon);
		xold = xnew;
		yold = ynew;
		xnew = xnewer;
		ynew = ynewer;
	    }
	    dev.vector (xold, yold, xnew, ynew, fat, dashon);
	    xold = xnew;
	    yold = ynew;
/*    fprintf(stderr,"in dovplot 9\n");*/
	    break;
	case VP_PMARK:		/* polymarker */
/*    fprintf(stderr,"in dovplot 10\n");*/
	    npts = geth (pltin);/* how many markers ? */
	    type = geth (pltin);/* what symbol? (any positive integer) */
	    size = geth (pltin);/* How big? */
	    size = size * mkscale * vdevscale * RPERIN / TXPERIN;

	    if (npts == 0)
		break;
	    /* allocate space for the points */
	    marker_vec = (int *) malloc ((unsigned) (npts * 2 * sizeof (int)));
	    mvec = marker_vec;
	    if (mvec == NULL)
		ERR (FATAL, name, "Can't malloc memory for markers!");

	    /* read the locations, and transform them to device coordinates */
	    for (ii = 0; ii < npts; ii++)
	    {
		GETXY (xnewer, ynewer);
		*mvec = xnewer;
		++mvec;
		*mvec = ynewer;
		++mvec;
	    }
	    update_color ();
	    /* call the device routine to display the markers */
	    dev.marker (npts, type, size, marker_vec);
	    /* release the storage used for the points */
	    free ((char *) marker_vec);
/*    fprintf(stderr,"in dovplot a\n");*/
	    break;
	case VP_ORIGIN:	/* set origin */
/*    fprintf(stderr,"in dovplot b\n");*/
	    xorigin = geth (pltin);
	    yorigin = geth (pltin);
/*    fprintf(stderr,"in dovplot c\n");*/
	    break;
	case VP_BEGIN_GROUP:
/*    fprintf(stderr,"in dovplot d\n");*/
	    ii = ftell (pltin) - 1;
	    getvpstring ();
	    strncpy (group_name, txbuffer, MAXFLEN);
	    dev.attributes (BEGIN_GROUP, group_number, ii, 0, 0);
	    group_number++;
/*    fprintf(stderr,"in dovplot e\n");*/
	    break;
	case VP_END_GROUP:
/*    fprintf(stderr,"in dovplot f\n");*/
	    group_number--;
	    if (group_number < 1)
	    {
		ERR (WARN, name,
		     "group invalidly nested");
		group_number = 1;
	    }
	    else
	    {
		ii = dev.attributes (END_GROUP, group_number, USER_EGROUP, 0, 0);
		if (ii == DOVPLOT_EXIT)
		{
		    dev.close (CLOSE_FLUSH);
		    return;
		}
		if (ii != DOVPLOT_CONT)
		    ERR (WARN, name, "dev.attributes(END_GROUP,...) returned junk.");
	    }
/*    fprintf(stderr,"in dovplot g\n");*/
	    break;
	case VP_GTEXT:		/* GKS-like text */
/*    fprintf(stderr,"in dovplot h\n");*/
	    xtext0 = 0;
	    ytext0 = 0;
	    vptodevxy_text (xtext0, ytext0, &xtext0, &ytext0);
	    GETXY_TEXT (xtext1, ytext1);
	    GETXY_TEXT (xtext2, ytext2);
    g_text:
	    xtext1 -= xtext0;
	    xtext2 -= xtext0;
	    ytext1 -= ytext0;
	    ytext2 -= ytext0;

	    savefat = fat;
	    savefatmult = fatmult;
	    fatmult *= txscale;
	    if (ifat >= 0)
	    {
		fat = fatmult * (float) (ifat + fatbase);
	    }
	    else
	    {
		fat = -1;
	    }
	    update_color ();
	    getvpstring ();
/*
 * Fonts less than NUMGENFONT reserved for gentext fonts:
 * up to the device to enforce that rule, though.
 */
	    dev.text (txbuffer,
		      txscale * (float) xtext1 / TEXTVECSCALE,
		      txscale * (float) ytext1 / TEXTVECSCALE,
		      txscale * (float) xtext2 / TEXTVECSCALE,
		      txscale * (float) ytext2 / TEXTVECSCALE);
	    fat = savefat;
	    fatmult = savefatmult;
/*    fprintf(stderr,"in dovplot 11\n");*/
	    break;
	case VP_TEXT:		/* text */
/*    fprintf(stderr,"in dovplot 12\n");*/
	    size = geth (pltin);
	    size = size * (float) RPERIN / (float) TXPERIN;
	    orient = (int) geth (pltin);
    new_text:
	    xtext0 = 0;
	    ytext0 = 0;
	    /* Character path direction */
	    xtext1 =
	     ROUND (TEXTVECSCALE * size * cos (orient * 3.14159 / 180.));
	    ytext1 =
	     ROUND (TEXTVECSCALE * size * sin (orient * 3.14159 / 180.));
	    /* Character up vector direction */
	    orient += 90;
	    xtext2 =
	     ROUND (TEXTVECSCALE * size * cos (orient * 3.14159 / 180.));
	    ytext2 =
	     ROUND (TEXTVECSCALE * size * sin (orient * 3.14159 / 180.));
	    vptodevxy_text (xtext0, ytext0, &xtext0, &ytext0);
	    vptodevxy_text (xtext1, ytext1, &xtext1, &ytext1);
	    vptodevxy_text (xtext2, ytext2, &xtext2, &ytext2);
/*    fprintf(stderr,"in dovplot i\n");*/
	    goto g_text;
	    break;
	case VP_OLDTEXT:	/* archaic format text */
/*    fprintf(stderr,"in dovplot 13\n");*/
	    if ((key = geth (pltin)) < 0)
		ERR (FATAL, name, "invalid text key");
	    size = (key & 037);
	    size = size * (float) RPERIN / (float) TXPERIN;
	    orient = (int) (((key & 0140) >> 5) * 90);
/*    fprintf(stderr,"in dovplot k\n");*/
	    goto new_text;
	    break;
	case VP_OLDAREA:	/* polygon */
/*    fprintf(stderr,"in dovplot j\n");*/
	    npts = geth (pltin);
	    afat = geth (pltin);
	    if (afat >= 0)
	    {
		afat = fatmult * (fatbase + afat);
	    }
	    nx_temp = geth (pltin);
	    ny_temp = geth (pltin);

/*
 * If a monochrome device, then fill with dot pattern.
 * If a color device, fill with solid color.
 */
	    nx = nx_temp * patternmult;
	    if (nx_temp == 1 || (nx_temp > 0 && nx == 0))
		nx = 1;
	    ny = ny_temp * patternmult;
	    if (ny_temp == 1 || (ny_temp > 0 && ny == 0))
		ny = 1;

	    nx_orig = nx;
	    ny_orig = ny;

	    if (!mono)
	    {
		if (nx_temp == 0 || ny_temp == 0)
		{
		    nx = 0;
		    ny = 0;
		}
		else
		{
		    nx = 1;
		    ny = 1;
		}
	    }
	    /*
	     * Create a temporary pattern
	     */
	    ipat = 0;
	    if (nx * ny > 0)
	    {
		if ((ptr = (int *) calloc ((unsigned) nx * ny, sizeof (int))) == NULL)
		    ERR (FATAL, name, "cannot alloc memory to load pattern");
		pat[ipat].patbits = ptr;
		ptr[(nx * ny) - 1] = cur_color;
	    }
	    else
		pat[ipat].patbits = NULL;

	    pat[ipat].xdim = nx;
	    pat[ipat].ydim = ny;
	    pat[ipat].xdim_orig = nx_orig;
	    pat[ipat].ydim_orig = ny_orig;
	    npts = getpolygon (npts);
	    update_color ();

/* Do the polygon */
	    if (npts > 2 && shade && nx > 0 && ny > 0)
		dev.area (npts, vxbuffer);

/* And then its border */
	    if (afat >= 0)
		vecoutline (vxbuffer);

	    if (nx * ny > 0)
		free ((char *) ptr);
/*    fprintf(stderr,"in dovplot l\n");*/
	    break;
	case VP_AREA:		/* polygon fill */
/*    fprintf(stderr,"in dovplot m\n");*/
	    npts = geth (pltin);
	    ipat = pat_color;

	    if (pat[ipat].patbits == NULL && shade)
	    {
		/*
		 * Create a default pattern (solid fill with this color) If
		 * you don't like this default for a particular device, ie, a
		 * black and white one, then the device itself can set up
		 * defaults for colors 0 through 7 in dev.open
		 */
		nx = 1;
		ny = 1;

		if ((ptr = (int *) calloc ((unsigned) nx * ny, sizeof (int))) == NULL)
		    ERR (FATAL, name, "cannot alloc memory to load pattern");
		pat[ipat].patbits = ptr;
		ptr[(nx * ny) - 1] = cur_color;
		pat[ipat].xdim = nx;
		pat[ipat].ydim = ny;
		pat[ipat].xdim_orig = nx;
		pat[ipat].ydim_orig = ny;
	    }

	    npts = getpolygon (npts);

	    if (!shade)
	    {
		/* At least draw the boundary to show where it is */
		afat = 0;
		update_color ();
		vecoutline (vxbuffer);
	    }
	    else
	    {
		/* See whether raster or hatch area */
		if (pat[ipat].xdim >= 0)
		{
		    /* raster */
		    if (npts > 2 && pat[ipat].ydim > 0 && pat[ipat].xdim > 0)
		    {
			update_color ();
			dev.area (npts, vxbuffer);
		    }
		}
		else
		{
		    /* hatch */
		    numhatch = -pat[ipat].xdim;
		    angle = 3.14159 * pat[ipat].ydim / 180.;
		    ptr = pat[ipat].patbits;
		    for (i = 0; i < numhatch * 2; i++)
		    {
			hafat[i] = (*ptr++);
			hacol[i] = (*ptr++);
			haoff[i] = (*ptr++);
			hasiz[i] = (*ptr++);
		    }
		    if (npts > 2)
		    {
			/*
			 * Hatch patterns don't rotate. The polygon does, the
			 * fill pattern doesn't.
			 */
			genhatch (npts, numhatch, angle, hafat, hacol, haoff, hasiz, vxbuffer);
		    }
		}
	    }
/*    fprintf(stderr,"in dovplot n\n");*/
	    break;
	case VP_FAT:		/* fat */
/*    fprintf(stderr,"in dovplot o\n");*/
	    /*
	     * negative fat always means to not draw the line at all.
	     */
	    ifat = geth (pltin);
	    if (ifat >= 0)
	    {
		fat = fatmult * (float) (ifat + fatbase);
	    }
	    else
	    {
		fat = -1;
	    }
	    dev.attributes (NEW_FAT, fat, 0, 0, 0);
/*    fprintf(stderr,"in dovplot p\n");*/
	    break;
	case VP_COLOR:		/* change color */
/*    fprintf(stderr,"in dovplot q\n");*/
	    pat_color = geth (pltin);
	    if (pat_color > MAX_COL || pat_color < 0)
	    {
		ERR (WARN, name, "bad color number %d (max %d, min 0)",
		     pat_color, MAX_COL);
		pat_color = DEFAULT_COLOR;
	    }
	    next_color = COLOR_MAP (pat_color);
	    if (next_color != cur_color)
	    {
		cur_color = next_color;
		need_devcolor = YES;
	    }
	    /*
	     * Pattern zero reserved for the OLD_AREA command, so increment
	     * the rest by one to make room.
	     */
	    pat_color++;
/*    fprintf(stderr,"in dovplot r\n");*/
	    break;
	case VP_SET_COLOR_TABLE:	/* set color table entry */
     set_table();
	    break;
	case VP_PURGE:		/* purge pltout buffers */
/*    fprintf(stderr,"in dovplot u\n");*/
	    dev.close (CLOSE_FLUSH);
/*    fprintf(stderr,"in dovplot v\n");*/
	    break;
	case VP_BREAK:		/* break */
/*    fprintf(stderr,"in dovplot w\n");*/
	    if (brake == BREAK_IGNORE)
/*    fprintf(stderr,"in dovplot x\n");*/
		break;
	    /* NOTE break IS IN AN "if": WE USUALLY FALL THROUGH HERE! */
	case VP_ERASE:		/* erase (and break falls through here too) */
/*    fprintf(stderr,"in dovplot y\n");*/
	    dev.close (CLOSE_FLUSH);

	    /*
	     * Erases and breaks can't occur inside groups.
	     */
	    group_number--;
	    if (group_number != 0)
	    {
		ERR (WARN, name,
		     "group contains erase or break");
		group_number = 0;
	    }
	    else
	    {
		if (erase & DO_LITERALS)
		{
		    if (c == VP_ERASE || brake)
			ii = dev.attributes (END_GROUP, group_number, ERASE_EGROUP, 0, 0);
		    else
			ii = dev.attributes (END_GROUP, group_number, BREAK_EGROUP, 0, 0);
		}
		else
		    ii = dev.attributes (END_GROUP, group_number, IGNORE_EGROUP, 0, 0);

		if (ii == DOVPLOT_EXIT)
		    return;
		if (ii != DOVPLOT_CONT)
		    ERR (WARN, name, "dev.attributes(END_GROUP,...) returned junk.");
	    }

	    if (epause < 0)
	    {
		message (MESG_ERASE);
		message (MESG_ON);
		message (MESG_HOME);
		message (MESG_READY);
		message (MESG_HIGHLIGHT_ON);
		message (MESG_TEXT, "Type Return to Continue...  ");
		message (MESG_DONE);
		ii = dev.interact (INT_PAUSE, controltty, string);
		if (ii == DOVPLOT_EXIT)
		    return;
		message (MESG_HIGHLIGHT_OFF);
		if (!allowecho)
		{
		    message (MESG_READY);
		    message (MESG_TEXT, CRLF);
		}
		message (MESG_DONE);
		message (MESG_OFF);
		message (MESG_ERASE);
	    }
	    else
	    {
		if (epause > 0)
		    sleep ((unsigned) epause);
	    }
	    /*
	     * Inquire point back from device
	     */
	    if (interact[0] != '\0')
	    {
		getapoint ();
	    }

	    if (erase & DO_LITERALS)
		if ((c == VP_ERASE) || brake)
		{
		    dev.erase (ERASE_MIDDLE);
		    nplots++;
		}
		else
		{
		    dev.erase (ERASE_BREAK);
		    nplots++;
		}

	    new_style = default_style;
	    setstyle (new_style);
	    reset ();

	    /*
	     * Start a new group level 0 to contain the next frame (separated
	     * from the previous one by either an erase or break)
	     */
	    ii = ftell (pltin);
	    sprintf (group_name, "%s.%d", pltname, nplots);
	    dev.attributes (BEGIN_GROUP, group_number, ii, 0, 0);
	    group_number++;

	    if (framewindows)
		outline_window ();

/*    fprintf(stderr,"in dovplot z\n");*/
	    break;
	case VP_WINDOW:	/* window */
/*    fprintf(stderr,"in dovplot aa\n");*/
	    if (window)
	    {
		xwmin = geth (pltin);
		ywmin = geth (pltin);
		xwmax = geth (pltin);
		ywmax = geth (pltin);

		if (xwmin > xwmax || ywmin > ywmax)
		{
/*
 * vptodevw will always return a window that is the "right way around",
 * even if the original window was inverted. So produce an inside-out
 * window which will clip everything away to nothing.
 */

		    xwmin = xWmax + 1;
		    xwmax = xWmin - 1;
		    ywmin = yWmax + 1;
		    ywmax = yWmin - 1;
		}
		else
		{

		    vptodevw (xwmin, ywmin, xwmax, ywmax, &xwmin, &ywmin, &xwmax, &ywmax);

		    wlimit (xWmin, xWmax, &xwmin, &xwmax);
		    wlimit (yWmin, yWmax, &ywmin, &ywmax);
		}
	    }
	    else
	    {
		geth (pltin);
		geth (pltin);
		geth (pltin);
		geth (pltin);
	    }
	    reset_windows ();
	    if (framewindows)
		outline_window ();
/*    fprintf(stderr,"in dovplot aa\n");*/
	    break;
	case VP_NOOP:		/* no op */
/*    fprintf(stderr,"in dovplot ab\n");*/
	    break;
	case VP_TXALIGN:	/* set text alignment */
/*    fprintf(stderr,"in dovplot ac\n");*/
	    /*
	     * These are made available to the device text routine
	     */
	    txalign.hor = geth (pltin);
	    txalign.ver = geth (pltin);
	    dev.attributes (NEW_ALIGN, txalign.hor, txalign.ver, 0, 0);
/*    fprintf(stderr,"in dovplot ad\n");*/
	    break;
	case VP_TXFONTPREC:	/* set text font */
/*    fprintf(stderr,"in dovplot ae\n");*/
	    /*
	     * These are made available to the device text routine
	     */
	    ii = geth (pltin);
	    if (ii >= 0)
		txfont = ii;
	    else
		ii = -1;

	    jj = geth (pltin);
	    if (jj >= 0)
		txprec = jj;
	    else
		jj = -1;

	    kk = geth (pltin);
	    if (kk >= 0)
		txovly = kk;
	    else
		kk = -1;

	    /*
	     * Another way for the device to keep track of changes.
	     */
	    dev.attributes (NEW_FONT, ii, jj, kk, 0);
/*    fprintf(stderr,"in dovplot af\n");*/
	    break;
	case VP_OVERLAY:	/* change overlay mode */
/*    fprintf(stderr,"in dovplot ag\n");*/
	    /*
	     * This is made available to the device dependent subroutines
	     */
	    overlay = geth (pltin);
	    dev.attributes (NEW_OVERLAY, overlay, 0, 0, 0);
/*    fprintf(stderr,"in dovplot ah\n");*/
	    break;
	case VP_PATLOAD:	/* load a pattern */
/*    fprintf(stderr,"in dovplot ai\n");*/
	    nmul = geth (pltin);
	    ny = geth (pltin);

	    /* See whether Raster or Hatch */
	    if (ny >= 0)
	    {
		/* Raster */
		/* nmul gives pixels_per_inch pattern is designed for */
		ny_mult = ny * (pixels_per_inch / nmul) * patternmult / aspect_ratio;
		if (ny_mult == 0 && ny > 0)
		    ny_mult = 1;
		nx = geth (pltin);
		nx_mult = nx * (pixels_per_inch / nmul) * patternmult;
		if (nx_mult == 0 && nx > 0)
		    nx_mult = 1;

		ipat = geth (pltin) + 1;
		if (ipat > NPAT || ipat < 1)
		    ERR (FATAL, name, "bad pattern number %d (max %d, min 1)", ipat, NPAT);
		/*
		 * Free up pattern that may already be there
		 */
		if (pat[ipat].patbits != NULL)
		{
		    free ((char *) pat[ipat].patbits);
		}

		if (nx_mult * ny_mult > 0)
		{
		    if ((ptr = (int *) malloc ((unsigned) (nx_mult * ny_mult * sizeof (int)))) == NULL)
			ERR (FATAL, name, "cannot alloc memory to load pattern");
		    pat[ipat].patbits = ptr;
		}
		else
		{
		    if ((ptr = (int *) malloc ((unsigned) (1 * sizeof (int)))) == NULL)
			ERR (FATAL, name, "cannot alloc memory to load dummy pattern");
		    pat[ipat].patbits = ptr;
		    ptr[0] = 0;
		}

		if (nx * ny > 0)
		{
		    if ((tempbuf = (int *) malloc ((unsigned) (nx * ny * sizeof (int)))) == NULL)
			ERR (FATAL, name,
			     "cannot alloc memory to load pattern's temporary buffer");
		}
		else
		    tempbuf = NULL;

		/*
		 * read in pattern
		 */
		ptemp = tempbuf;
		for (j = 0; j < nx * ny; j++)
		{
		    k = geth (pltin);
		    if (k > MAX_COL || k < 0)
		    {
			ERR (WARN, name, "bad color number in pattern %d (max %d, min 0)",
			     k, MAX_COL);
			k = DEFAULT_COLOR;
		    }
		    *ptemp++ = COLOR_MAP (k);
		}

		/*
		 * copy into patbits, with "stretching"
		 */
		for (i = 0; i < ny_mult; i++)
		{
		    ny_temp = i * ny / ny_mult;
		    for (j = 0; j < nx_mult; j++)
		    {
			nx_temp = j * nx / nx_mult;
			ptr[j + nx_mult * i] = tempbuf[nx_temp + nx * ny_temp];
		    }
		}
		/*
		 * set dimensions of pattern
		 */
		pat[ipat].xdim = nx_mult;
		pat[ipat].ydim = ny_mult;
		pat[ipat].xdim_orig = nx_mult;
		pat[ipat].ydim_orig = ny_mult;

		if (tempbuf != NULL)
		    free ((char *) tempbuf);

		dev.attributes (NEW_PAT, ipat, 0, 0, 0);
	    }
	    else
	    {
		/* Hatch Pattern */
		/* nmul gives angle, ny is merely a flag */
		nx = geth (pltin);
		if (nx <= 0 || nx * 2 > NHATCH)
		    ERR (FATAL, name, "bad numhatch %d (max %d/2, min 1)", nx, NHATCH);
		ipat = geth (pltin) + 1;
		if (ipat > NPAT || ipat < 1)
		    ERR (FATAL, name, "bad pattern number %d (max %d, min 1)", ipat, NPAT);
		/*
		 * Free up pattern that may already be there
		 */
		if (pat[ipat].patbits != NULL)
		{
		    free ((char *) pat[ipat].patbits);
		}
		if ((ptr = (int *) malloc ((unsigned) (nx * 4 * 2 * sizeof (int)))) == NULL)
		    ERR (FATAL, name, "cannot alloc memory to load pattern");
		pat[ipat].patbits = ptr;

		for (i = 0; i < nx * 2; i++)
		{
		    hafat[i] = geth (pltin);
		    if (hafat[i] >= 0)
		    {
			hafat[i] = fatmult * (fatbase + hafat[i]);
		    }
		    k = geth (pltin);
		    if (k > MAX_COL || k < 0)
		    {
			ERR (WARN, name, "bad color number in hatch %d (max %d, min 0)",
			     k, MAX_COL);
			k = DEFAULT_COLOR;
		    }
		    hacol[i] = COLOR_MAP (k);
		    haoff[i] = geth (pltin) * patternmult * pixels_per_inch / RPERIN;
		    hasiz[i] = geth (pltin);
		}
		/*
		 * Find the smallest hatch interval. 1/2 that, and then force
		 * everything to be a multiple of that. If this were not
		 * done, then hatch intervals that originally differed by
		 * some simple fraction might end up slightly off, causing
		 * very unsightly beats.
		 */
		/*
		 * Upper quantization limit of 1/10 inch
		 */
		k = (RPERIN / 10) * 2;
		for (i = 0; i < nx * 2; i++)
		{
		    if (hasiz[i] > 0 && hasiz[i] < k)
			k = hasiz[i];
		}
		j = k * (patternmult * pixels_per_inch / RPERIN) / 2.;
		if (j < 1)
		    j = 1;
		for (i = 0; i < nx * 2; i++)
		{
		    hasiz[i] = ((int) ((hasiz[i] * 2) / k)) * j;
		}
		/*
		 * The above algorithm also means that you can't have a hatch
		 * pattern come out on any device with a repetition rate of
		 * faster than once per two pixels.
		 */

		for (i = 0; i < nx * 2; i++)
		{
/*
 * Make sure haoff < hasiz
 */
		    if (haoff[i] >= hasiz[i])
			haoff[i] = 0;
		    *ptr++ = hafat[i];
		    *ptr++ = hacol[i];
		    *ptr++ = haoff[i];
		    *ptr++ = hasiz[i];
		}
		/*
		 * set numhatch and angle... dimensions are 2 * numhatch by 4
		 * . pat[ipat].xdim negative is a flag that this is a hatch
		 * pattern, not raster, and so must be treated differently.
		 */
		/*
		 * numhatch, with neg as flag for hatch numhatch = 0 is OK,
		 * because that means don't fill, same as nx = 0.
		 */
		pat[ipat].xdim = -nx;
		pat[ipat].ydim = nmul;	/* angle */
	    }
	    break;
	case VP_SHORT_RASTER:	/* bit raster data */
	case VP_BIT_RASTER:	/* bit raster data */
	case VP_BYTE_RASTER:	/* byte raster data */
/*   fprintf(stderr,"in byte cxase \n");*/
      if (c == VP_SHORT_RASTER) byte2=1;
      else byte2=0;
	    ras_orient = geth (pltin);
	    if (rotate % 90 != 0)
	    {
		if (wantras)
		{
		    ERR (WARN, name, "Raster only possible in 4 principal orientations");
		    wantras = NO;
		}
	    }
	    else
	    {
		ras_orient += rotate / 90;
		if (ras_orient >= 0)
		    ras_orient = ras_orient % 4;
		else
		    ras_orient = ((ras_orient % 4) + 4) % 4;
	    }

	    /*
	     * They're on their honor to not go out of bounds. This check is
	     * just for things that HAVE to go out of bounds.
	     */
	    ras_offset = geth (pltin);

	    if (ras_offset + 0 > MAX_COL || ras_offset + 255 < 0)
	    {
		ERR (FATAL, name, "Absurd raster offset %d", ras_offset);
	    }

	    xvr_min = geth (pltin);
	    yvr_min = geth (pltin);
	    xvr_max = geth (pltin);
	    yvr_max = geth (pltin);
/*      fprintf(stderr,"what the problem %d %d %d %d \n",xvr_min,yvr_min,xvr_max,yvr_max);*/
	    vptodevw (xvr_min, yvr_min, xvr_max, yvr_max,
		      &xvr_min, &yvr_min, &xvr_max, &yvr_max);
	    xvru_min = xvr_min;
	    yvru_min = yvr_min;
	    xvru_max = xvr_max;
	    yvru_max = yvr_max;

	    xpix = geth (pltin);
	    ypix = geth (pltin);
/*     fprintf(stderr,"retrieve x and y %d %d \n",xpix,ypix);*/

	    switch (ras_orient)
	    {
	    case 0:
		xrasmult = (float) xpix / (float) (xvr_max - xvr_min);
		yrasmult = (float) ypix / (float) (yvr_max - yvr_min);
		xvr_max--;
		yvr_max--;
		break;
	    case 1:
		yrasmult = (float) ypix / (float) (xvr_max - xvr_min);
		xrasmult = (float) xpix / (float) (yvr_max - yvr_min);
		xvr_max--;
		yvr_min++;
		break;
	    case 2:
		xrasmult = (float) xpix / (float) (xvr_max - xvr_min);
		yrasmult = (float) ypix / (float) (yvr_max - yvr_min);
		xvr_min++;
		yvr_min++;
		break;
	    case 3:
		yrasmult = (float) ypix / (float) (xvr_max - xvr_min);
		xrasmult = (float) xpix / (float) (yvr_max - yvr_min);
		xvr_min++;
		yvr_max--;
		break;
	    }
/*    fprintf(stderr,"ORRIENT=%d xpix=%d ypix=%d xvr_min=%d xvr_max=%d yvr_min=%d yvr_max=%d \n",orient,xpix,ypix,xvr_min,xvr_max,yvr_min,yvr_max);*/

	    if (wantras && smart_raster)
	    {
		rasterline = (unsigned short *) malloc ((unsigned)((xpix + 7 + 2) * sizeof (unsigned short)));
		rasterline2 = (unsigned short *) malloc ((unsigned) ((xpix + 7 + 2) * sizeof (unsigned short)));
		if (rasterline == NULL || rasterline2 == NULL)
		    ERR (FATAL, name, "cannot alloc memory to load raster line");

		outraster = (unsigned short *) malloc ((xpix * ypix) * sizeof (unsigned short));
		if (outraster == NULL)
		    ERR (FATAL, name, "cannot alloc memory to load raster image");

/*
 * See whether we want to dither or not.
 * Dither if monochrome device and dithering has been asked for.
 * It's up to the device to decide whether to actually do this.
 */
		dither_it = dither && mono;

/*
 * If the device is color, but the raster happens to all be shades
 * of grey, let the device know so it can take advantage of the
 * fact. This is available to the routine via the external variable
 * "ras_allgrey".	joe & david 10/03/94
 */

              ras_allgrey = YES;
              if (!mono)
              {
                  if (c == VP_BIT_RASTER)
                  {
                      for (ii=0; ii<2; ii++)
                      {
                          if (
                          (color_set [ras_offset * ii][_RED] !=
                          color_set [ras_offset * ii][_GREEN]) ||
                          (color_set [ras_offset * ii][_GREEN] !=
                          color_set [ras_offset * ii][_BLUE])
                          )
                          {
                              ras_allgrey = NO;
                              break;
                          }
                      }
                  }
                  else
                  {
                      for (ii=0; ii<256; ii++)
                      {
                          if (
                          (color_set [ras_offset + ii][_RED] !=
                          color_set [ras_offset + ii][_GREEN]) ||
                          (color_set [ras_offset + ii][_GREEN] !=
                          color_set [ras_offset + ii][_BLUE])
                          )
                          {
                              ras_allgrey = NO;
                              break;
                          }
                      }
                  }
              }

		/*
		 * Read in the Raster data for "Smart" devices, ie, those
		 * which can stretch (and dither) their own raster.
		 */
		num_rep = 0;
		for (yrast = 0; yrast < ypix; yrast++)
		{
		    /*
		     * Read in the next raster line, if we have a new one
		     */
		    if (num_rep <= 0)
		    {
			num_rep = geth (pltin);
			if (num_rep <= 0){
          fprintf(stderr,"In multiplier error 1 \n");
			    ERR (FATAL, name, "Bad Raster line multiplier(1)");
      }
			pos = 0;
		new_pat:num_pat = geth (pltin);
			num_byte = geth (pltin);
			if (num_pat <= 0 || num_byte <= 0 ||
			    pos + num_pat * num_byte > xpix){
/*          fprintf(stderr,"problem 1 %d %d %d %d \n",num_pat,num_byte,pos,xpix);*/
			    ERR (FATAL, name, "Raster line not length promised");
        }

			if (num_pat > 1)
			{
			    if (dither_it)
			    {
/*fprintf(stderr,"in this read 1 \n");*/
				READ_RASTER (
					     rasterline2[j] = GREY_MAP (ras_offset + (unsigned short) fgetc (pltin)),
					     rasterline2[j] = GREY_MAP (ras_offset + (unsigned short) fgetc (pltin)*256+fgetc(pltin)),
					     rasterline2[j + jj] = GREY_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
				 );
			    }
			    else
			    {
/*fprintf(stderr,"in this read 2 \n");*/
				READ_RASTER (
					     rasterline2[j] = COLOR_MAP (ras_offset + (unsigned short) fgetc (pltin)),
					     rasterline2[j] = COLOR_MAP (ras_offset + (unsigned short) fgetc (pltin)*256+fgetc(pltin)),
					     rasterline2[j + jj] = COLOR_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
				 );
			    }
			    for (j = 0; j < num_pat; j++)
			    {
				for (k = 0; k < num_byte; k++)
				{
				    rasterline[pos] = rasterline2[k];
				    pos++;
				}
			    }
			}
			else
			{
			    if (dither_it)
			    {
/*fprintf(stderr,"in this read 3 \n");*/
				READ_RASTER (
					     rasterline[pos + j] = GREY_MAP (ras_offset + (unsigned short) fgetc (pltin)),
					     rasterline[pos + j] = GREY_MAP (ras_offset + (unsigned short) fgetc (pltin)*256+fgetc(pltin)),
					     rasterline[pos + j + jj] = GREY_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
				 );
			    }
			    else
			    {
/*fprintf(stderr,"in this read 4 \n");*/
				READ_RASTER (
					     rasterline[pos + j] = COLOR_MAP (ras_offset + (unsigned short) fgetc (pltin)),
					     rasterline[pos + j] = COLOR_MAP (ras_offset + (unsigned short) fgetc (pltin)*256+fgetc(pltin)),
					     rasterline[pos + j + jj] = COLOR_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
				 );
			    }
			    pos += num_byte;
			}

			if (pos < xpix)
			    goto new_pat;
			if (pos != xpix){
/*          fprintf(stderr,"problem 2 %d %d \n",pos,xpix);*/
			    ERR (FATAL, name, "Raster line not length promised");
      }
		    }
		    num_rep--;
/*       fprintf(stderr,"where do i die 1a \n");*/
		    for (ii = 0; ii < xpix; ii++)
		    {
			outraster[ii + yrast * xpix] = rasterline[ii];
/*            fprintf(stderr,"out raster check 2 %d %d %d \n",*/
/*              ii,(int) ((ii)),outraster2[ii+yrast*xpix]);*/
		    }
		}

		/* Smart form */
		dev.raster (xpix, ypix, xvr_min, yvr_min, xvr_max, yvr_max,
			    outraster, ras_orient, dither_it,byte2);

		free ((char *) rasterline);
		free ((char *) rasterline2);
		free ((char *) outraster);
	    }
	    else
	    {
		wlimit (xwmin, xwmax, &xvr_min, &xvr_max);
		wlimit (ywmin, ywmax, &yvr_min, &yvr_max);

		switch (ras_orient)
		{
		case 0:
		    xvr_max++;
		    yvr_max++;
		    xr_max = xvr_max - xvru_min;
		    xr_min = xvr_min - xvru_min;
		    yr_max = yvr_max - yvru_min;
		    yr_min = yvr_min - yvru_min;
		    yru_max = yvru_max - yvru_min;
		    break;
		case 1:
		    xvr_max++;
		    yvr_min--;
		    xr_max = yvru_max - yvr_min;
		    xr_min = yvru_max - yvr_max;
		    yr_max = xvr_max - xvru_min;
		    yr_min = xvr_min - xvru_min;
		    yru_max = xvru_max - xvru_min;
		    break;
		case 2:
		    xvr_min--;
		    yvr_min--;
		    xr_max = xvru_max - xvr_min;
		    xr_min = xvru_max - xvr_max;
		    yr_max = yvru_max - yvr_min;
		    yr_min = yvru_max - yvr_max;
		    yru_max = yvru_max - yvru_min;
		    break;
		case 3:
		    xvr_min--;
		    yvr_max++;
		    xr_max = yvr_max - yvru_min;
		    xr_min = yvr_min - yvru_min;
		    yr_max = xvru_max - xvr_min;
		    yr_min = xvru_max - xvr_max;
		    yru_max = xvru_max - xvru_min;
		    break;
		}
		xru_min = 0;

		if (yr_max < yr_min || xr_max < xr_min || !wantras)
		{
		    /*
		     * We need to read through all the raster stuff, even if
		     * we never use it.
		     */
		    yr_max = yr_min;
		    xr_max = xr_min;
		}

		rasterline = (unsigned short *) malloc (((xpix + 7 + 2) * sizeof (unsigned short)));
		rasterline2 = (unsigned short *) malloc (((xpix + 7 + 2) * sizeof (unsigned short)));
		if (rasterline == NULL || rasterline2 == NULL)
		    ERR (FATAL, name, "cannot alloc memory to load raster line");

/*
 * See whether we need to dither or not.
 * Dither if monochrome device and dithering has been asked for.
 */
		dither_it = dither && mono;

		if (xr_max > xr_min)
		{
		    outraster2 = (unsigned short *) malloc ( ((xr_max - xr_min) * sizeof (unsigned short)));
		    if (dither_it)
		    {
			outraster = (unsigned short *) malloc ( ((xr_max - xr_min) *  sizeof (unsigned short)));
		    }
		    else
		    {
			outraster = outraster2;
		    }
		    if (outraster2 == NULL || outraster == NULL)
			ERR (FATAL, name, "cannot alloc memory to load raster line");
		}
		else
		{
		    outraster2 = NULL;
		    outraster = NULL;
		}

		/*
		 * Read in the Raster data
		 */
		lastrast = -1;
		num_rep = 0;
		for (i = yr_max - 1; i >= yr_min - 1; i--)
		{
		    yrast = (yru_max - 1 - i) * yrasmult;
		    if (i == yr_min - 1)
		    {
			/*
			 * Assure that the last bit of unused raster, if any,
			 * is read. This last time through the loop is a
			 * "dummy".
			 */
			yrast = ypix - 1;
		    }

		    for (ii = 0; ii < (yrast - lastrast); ii++)
		    {
			/*
			 * Read in the next raster line, if we have a new one
			 */
			if (num_rep <= 0)
			{
			    num_rep = geth (pltin);
			    if (num_rep <= 0){
/*           fprintf(stderr,"in  bad raster2 %d %d %d %d \n",num_rep,yr_min,i,yr_max); */
				ERR (FATAL, name, "Bad Raster line multiplier(2)");
        }
			    pos = 0;
		    new_pat2:num_pat = geth (pltin);
			    num_byte = geth (pltin);
/*          fprintf(stderr,"problem check %d %d %d %d %d \n",num_pat,num_byte,pos,xpix,num_rep);*/
			    if (num_pat <= 0 || num_byte <= 0 ||
				pos + num_pat * num_byte > xpix){
/*          fprintf(stderr,"problem 3 %d %d %d %d \n",num_pat,num_byte,pos,xpix);*/
				ERR (FATAL, name, "Raster line not length promised");
       }

/*      fprintf(stderr,"IF CHECK %d+%d >= %d -%d && %d >%d && %d >= %d \n",*/
/*          ii,num_rep,yrast-lastrast,xr_max,xr_min,i,yr_min);*/
/*
 * Only bother with it if we're actually going to use it
 */
/*			    if ((ii + num_rep >= yrast - lastrast*/
/*				&& xr_max > xr_min && i >= yr_min) || c==VP_SHORT_RASTER)*/
			    if ((ii + num_rep >= yrast - lastrast
				&& xr_max > xr_min && i >= yr_min) )
			    {
				if (num_pat > 1)
				{
				    if (dither_it)
				    {
/*fprintf(stderr,"in this read 5 \n");*/
					READ_RASTER (
						     rasterline2[j] = GREY_MAP (ras_offset + (unsigned short) fgetc (pltin)),
						     rasterline2[j] = GREY_MAP (ras_offset + (unsigned short) fgetc (pltin)*256+fgetc(pltin)),
						     rasterline2[j + jj] = GREY_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
					 );
				    }
				    else
				    {
/*fprintf(stderr,"in this read 6 \n");*/
					READ_RASTER (
						     rasterline2[j] = COLOR_MAP (ras_offset + (unsigned short) fgetc (pltin)*256+fgetc(pltin)),
						     rasterline2[j] = COLOR_MAP (ras_offset + (unsigned short) fgetc (pltin)),
						     rasterline2[j + jj] = COLOR_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
					 );
				    }
				    for (j = 0; j < num_pat; j++)
				    {
					for (k = 0; k < num_byte; k++)
					{
					    rasterline[pos] = rasterline2[k];
					    pos++;
					}
				    }
				}
				else
				{
				    if (dither_it)
				    {
/*fprintf(stderr,"in this read 7 \n");*/
					READ_RASTER (
						     rasterline[pos + j] = GREY_MAP (ras_offset + (unsigned short) fgetc (pltin)),
						     rasterline[pos + j] = GREY_MAP (ras_offset + (unsigned short)fgetc (pltin)*256+fgetc(pltin)),
						     rasterline[pos + j + jj] = GREY_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
					 );
				    }
				    else
				    {
/*						     rasterline[pos + j] = fgetc (pltin),*/
/*						     rasterline[pos + j] =( (unsigned short) fgetc (pltin)*256+fgetc(pltin)),*/
/*fprintf(stderr,"in this read 8 \n");*/
					READ_RASTER (
						     rasterline[pos + j] = COLOR_MAP (ras_offset + (unsigned short) fgetc (pltin)),
						     rasterline[pos + j] = COLOR_MAP (ras_offset + (unsigned short) fgetc (pltin)*256+fgetc(pltin)),
						     rasterline[pos + j + jj] = COLOR_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
					 );
				    }
				    pos += num_byte;
				}
			    }
			    else
			    {
/*
 * We're just going to turn right around and read another one,
 * So throw this away!
 */
				if (c == VP_BYTE_RASTER || c==VP_SHORT_RASTER)
				{
				    for (j = 0; j < num_byte*(byte2+1); j++)
				    {
					fgetc (pltin);
				    }
				}
				else
				{
				    for (j = 0; j < num_byte; j += 8)
				    {
					fgetc (pltin);
				    }
				}
				pos += num_pat * num_byte;
			    }
			    if (pos < xpix)
				goto new_pat2;
			    if (pos != xpix){
/*          fprintf(stderr,"problem 3 %d %d \n",pos,xpix);*/
				ERR (FATAL, name, "Raster line not length promised");
         }

			    if (wantras && xr_max > xr_min && i >= yr_min)
			    {
				for (j = xr_min; j < xr_max; j++)
				{
				    outraster2[j - xr_min] =
				     rasterline[(int) ((j - xru_min) * xrasmult)];
/*            fprintf(stderr,"out raster check 1 %d %d %d -%f \n",*/
/*              j,(int) ((j - xru_min) * xrasmult),outraster2[j - xr_min],xrasmult);*/
				}
			    }
			}
			num_rep--;
		    }
		    lastrast = yrast;

		    if (xr_max > xr_min && i >= yr_min)
		    {
			if (dither_it)
			{
			    dithline (outraster2, outraster, xr_max - xr_min, yr_max - 1 - i, dither);
			}
			/* Dumb forms */
      
			switch (ras_orient)
			{
			case 0:
			    dev.raster (yr_max - 1 - i, yr_max - yr_min,
					xvr_min, i + yvru_min,
				     xr_max - xr_min, ras_orient, outraster,
					0, 0,byte2);
			    break;
			case 1:
			    dev.raster (yr_max - 1 - i, yr_max - yr_min,
					xvru_min + i, yvr_max,
				     xr_max - xr_min, ras_orient, outraster,
					0, 0,byte2);
			    break;
			case 2:
			    dev.raster (yr_max - 1 - i, yr_max - yr_min,
					xvr_max, yvru_max - i,
				     xr_max - xr_min, ras_orient, outraster,
					0, 0,byte2);
			    break;
			case 3:
			    dev.raster (yr_max - 1 - i, yr_max - yr_min,
					xvru_max - i, yvr_min,
				     xr_max - xr_min, ras_orient, outraster,
					0, 0,byte2);
			    break;
			}
		    }
		}
		free ((char *) rasterline);
		free ((char *) rasterline2);
		if (outraster2 != NULL)
		{
		    free ((char *) outraster2);
		    if (dither_it)
			free ((char *) outraster);
		}
		if (!wantras)
		{
		    xxx[0] = xvru_min;
		    yyy[0] = yvru_min;
		    xxx[1] = xvru_min;
		    yyy[1] = yvru_max;
		    xxx[2] = xvru_max;
		    yyy[2] = yvru_max;
		    xxx[3] = xvru_max;
		    yyy[3] = yvru_min;
		    drawpolygon (4, xxx, yyy);
		}
	    }
/*   fprintf(stderr,"out byte cxase \n");*/
	    break;
	case VP_MESSAGE:	/* Text message */
	    getvpstring ();
	    message (MESG_READY);
	    message (MESG_MESSAGE);
	    message (MESG_TEXT, txbuffer);
	    message (MESG_TEXT, CRLF);
	    message (MESG_DONE);
	    break;
	case VP_SETDASH:
	    npts = geth (pltin);
	    if (npts > MAXDASH)
	    {
		ERR (FATAL, name, "Too complicated a dash line pattern.");
	    }
	    dashon = npts;
	    k = 0;
	    dashsum = 0.;
	    for (ii = 0; ii < npts * 2; ii++)
	    {
		dashes[ii] = dashscale * (float) geth (pltin) / RPERIN;
		if (dashes[ii] < 0.)
		    ERR (FATAL, name, "Negative dash distance.");

		if (dashes[ii] != 0.)
		    k = 1;

		dashsum += dashes[ii];
	    }
	    if (!k)
		dashon = NO;
	    dev.attributes (NEW_DASH, dashon, 0, 0, 0);
	    break;
	default:		/* error */
	    ERR (FATAL, name,
		 "invalid VPLOT command decimal %d character %c",
		 (int) c, (char) c);
	    break;
	}
    }
End_of_file:
 end_of_file();
 return;
}
void end_of_file(){
int ii;
/*
 * End of main while loop. Either fall out here or jump here when you hit
 * the end of the file somewhere.
 */

    dev.close (CLOSE_FLUSH);	/* force last vector out */

/*
 * End the group for this frame
 */
    group_number--;
    if (group_number != 0)
    {
	ERR (WARN, name,
	     "group left unclosed at end of file");
	group_number = 0;
    }
    else
    {
	if ((erase & FORCE_INITIAL) && (erase & DO_LITERALS))
	    ii = dev.attributes (END_GROUP, group_number, EOF_ERASE_EGROUP, 0, 0);
	else
	    ii = dev.attributes (END_GROUP, group_number, EOF_IGNORE_EGROUP, 0, 0);

	if (ii == DOVPLOT_EXIT)
	    return;
	if (ii != DOVPLOT_CONT)
	    ERR (WARN, name, "dev.attributes(END_GROUP,...) returned junk.");
    }

/*
 * Exit of dovplot
 */
    return;


}

/*
 * reset variables that can be affected by vplot commands when
 * processing multiple plots, and don't stay set across pages
 */
reset ()
{
int             ii, jj, kk;

    xwmin = xWmin;		/* plot window parameters defaulted */
    xwmax = xWmax;		/* to maximum size	 */
    ywmin = yWmin;
    ywmax = yWmax;

    reset_windows ();

    fat = fatmult * (float) (fatbase);
    dev.attributes (NEW_FAT, fat, 0, 0, 0);

    if (cur_color != COLOR_MAP (DEFAULT_COLOR));
    {
	need_devcolor = YES;
	cur_color = COLOR_MAP (DEFAULT_COLOR);
    }

    txalign.hor = TH_NORMAL;
    txalign.ver = TV_NORMAL;
    dev.attributes (NEW_ALIGN, txalign.hor, txalign.ver, 0, 0);

    ii = -1;
    jj = -1;
    kk = -1;
    if (txfont != default_txfont)
    {
	txfont = default_txfont;
	ii = txfont;
    }
    if (txprec != default_txprec)
    {
	txprec = default_txprec;
	jj = txprec;
    }
    if (txovly != default_txovly)
    {
	txovly = default_txovly;
	kk = txovly;
    }
    dev.attributes (NEW_FONT, ii, jj, kk, 0);

    dashon = NO;
    dev.attributes (NEW_DASH, dashon, 0, 0, 0);

    overlay = default_overlay;
    dev.attributes (NEW_OVERLAY, overlay, 0, 0, 0);
}

reset_windows ()
{
extern int      xwmax_last, ywmax_last, xwmin_last, ywmin_last;

    if (xwmax != xwmax_last || ywmax != ywmax_last
	|| xwmin != xwmin_last || ywmin != ywmin_last)
	dev.attributes (SET_WINDOW, xwmin, ywmin, xwmax, ywmax);

    xwmin_last = xwmin;
    ywmin_last = ywmin;
    xwmax_last = xwmax;
    ywmax_last = ywmax;
}

outline_window ()
{
    if (need_devcolor == YES || cur_color != COLOR_MAP (DEFAULT_COLOR))
    {
	dev.attributes (SET_COLOR, COLOR_MAP (DEFAULT_COLOR), 0, 0, 0);
	need_devcolor = NO;
    }
    dev.vector (xwmin, ywmin, xwmax, ywmin, 0, 0);
    dev.vector (xwmax, ywmin, xwmax, ywmax, 0, 0);
    dev.vector (xwmax, ywmax, xwmin, ywmax, 0, 0);
    dev.vector (xwmin, ywmax, xwmin, ywmin, 0, 0);
    if (cur_color != COLOR_MAP (DEFAULT_COLOR))
    {
	dev.attributes (SET_COLOR, cur_color, 0, 0, 0);
	need_devcolor = NO;
    }
}

getvpstring ()
{
char           *txptr;
int             ii;

    txptr = txbuffer;

    while (1)
    {
	ii = getc (pltin);

	if (ii == EOF)
	{
	    ERR (FATAL, name,
		 "Unexpected EOF encountered in Text string");
	}

	*txptr = ii;
	txptr++;

	if (ii == 0)
	{
	    break;
	}

	if ((txptr - txbuffer) == txbuflen)
	{
	    txbuffer = (char *) realloc (txbuffer, (unsigned) (txbuflen + TXBUFLEN));
	    txptr = txbuffer + txbuflen;
	    txbuflen += TXBUFLEN;
	}
    }
    return;
}

drawpolygon (npts, x, y)
    int             npts, *x, *y;
{
int             i, j;
struct vertex  *vertex;
static int      point;

    j = 0;
    if (npts > (vxbuflen - 1))
    {
	free ((char *) vxbuffer);
	vxbuffer =
	 (struct vertex *) malloc ((unsigned) ((npts + 1) * sizeof (struct vertex)));
    }
    vertex = vxbuffer;
    xnew = x[j];
    ynew = y[j];
    j++;
    vertex->x = xnew;
    vertex->y = ynew;
    vertex->next = vertex + 1;
    vertex++;
    i = npts - 1;
    while (i--)
    {
	vertex->x = x[j];
	vertex->y = y[j];
	j++;
	if ((vertex->x != xnew) || (vertex->y != ynew))
	{
	    xnew = vertex->x;
	    ynew = vertex->y;
	    vertex->next = vertex + 1;
	    vertex->last = vertex - 1;
	    vertex++;
	    continue;
	}
	npts--;
    }
    vertex--;
    vxbuffer->last = vertex;
    vertex->next = vxbuffer;

    ipat = 0;
    point = cur_color;
    pat[ipat].patbits = &point;
    pat[ipat].xdim = 1;
    pat[ipat].ydim = 1;
    pat[ipat].xdim_orig = 1;
    pat[ipat].ydim_orig = 1;
    update_color ();
    if (npts > 2)
    {
	if (shade)
	    dev.area (npts, vxbuffer);
	else
	    vecoutline (vxbuffer);
    }
}

getpolygon (npts)
    int             npts;
{
int             i;
struct vertex  *vertex;

    if (npts > (vxbuflen - 1))
    {
	free ((char *) vxbuffer);
	vxbuffer =
	 (struct vertex *) malloc ((unsigned) ((npts + 1) * sizeof (struct vertex)));
    }
    vertex = vxbuffer;
    GETXY (xnew, ynew);
    vertex->x = xnew;
    vertex->y = ynew;
    vertex->next = vertex + 1;
    vertex++;
    i = npts - 1;
    while (i--)
    {
	GETXY (vertex->x, vertex->y);
	if ((vertex->x != xnew) || (vertex->y != ynew))
	{
	    xnew = vertex->x;
	    ynew = vertex->y;
	    vertex->next = vertex + 1;
	    vertex->last = vertex - 1;
	    vertex++;
	    continue;
	}
	npts--;
    }
    vertex--;
    vxbuffer->last = vertex;
    vertex->next = vxbuffer;
    return (npts);
}

update_color ()
{
    if (need_devcolor)
    {
	dev.attributes (SET_COLOR, cur_color, 0, 0, 0);
	need_devcolor = NO;
    }
}

getapoint ()
{
    while (1)
    {
	xret = dev_xmax + 1;
	yret = dev_ymax + 1;
	if ((int) dev.getpoint (controltty, &xret, &yret) ||
	    (xret > dev_xmax - 5 && yret > dev_ymax - 5))
	    break;
	devtovpxy (xret, yret, &xret, &yret);
	add_a_cor (interact, xret, yret);
    }
}


void set_table(void){
int             col_tab_no, red, green, blue, grey, dist, min_dist, best_col;
float           red_f, green_f, blue_f, red_fm, green_fm, blue_fm;
int             red_mask, green_mask, blue_mask;
int             red_0, green_0, blue_0;
int             red_7, green_7, blue_7;
int             c,ii,k,i;
/*    fprintf(stderr,"in dovplot s\n");*/
	    /*
	     * The logic here is a bit tricky. Basically, it goes like this:
	     * If the device actually has a settable color of that number,
	     * then reset the device's color as asked. Otherwise, take the
	     * closest color that we have and map to that. Color 0 is the
	     * background color, and is special. Merely "close" colors are
	     * not mapped to Color 0... it has to be exact. Even devices with
	     * NO settable colors are still sent colors 0 through 7, and
	     * these colors are then always left set to their original
	     * defaults. In this way you can handle terminals like the GIGI
	     * which have colors but not SETTABLE colors. Also remember that
	     * the device itself will THEN have to permute colors 0 through 7
	     * on top of all this to correspond with Vplot's definitions!
	     */
	    col_tab_no = geth (pltin);
	    if (col_tab_no > MAX_COL || col_tab_no < 0)
	    {
		ERR (FATAL, name, "Bad color table value %d (%d max)",
		     col_tab_no, MAX_COL);
	    }
	    red = geth (pltin);
	    green = geth (pltin);
	    blue = geth (pltin);
	    if (red > MAX_GUN || green > MAX_GUN || blue > MAX_GUN
		|| red < 0 || green < 0 || blue < 0)
	    {
		ERR (FATAL, name,
		     "Bad color level in color %d (%d,%d,%d), %d max",
		     col_tab_no, red, green, blue, MAX_GUN);
	    }

/*
 * Fun and games with color values
 */
	    red_f = (float) red / (float) MAX_GUN;
	    green_f = (float) green / (float) MAX_GUN;
	    blue_f = (float) blue / (float) MAX_GUN;

	    red_fm =
	     redmap[3] + redmap[0] * red_f + redmap[1] * green_f + redmap[2] * blue_f;
	    if (red_fm < 0.)
		red_fm = 0.;
	    if (red_fm > 1.)
		red_fm = 1.;
	    red_fm = pow (red_fm, redpow);
	    red = ROUND (red_fm * MAX_GUN);

	    green_fm =
	     greenmap[3] + greenmap[0] * red_f + greenmap[1] * green_f + greenmap[2] * blue_f;
	    if (green_fm < 0.)
		green_fm = 0.;
	    if (green_fm > 1.)
		green_fm = 1.;
	    green_fm = pow (green_fm, greenpow);
	    green = ROUND (green_fm * MAX_GUN);

	    blue_fm =
	     bluemap[3] + bluemap[0] * red_f + bluemap[1] * green_f + bluemap[2] * blue_f;
	    if (blue_fm < 0.)
		blue_fm = 0.;
	    if (blue_fm > 1.)
		blue_fm = 1.;
	    blue_fm = pow (blue_fm, bluepow);
	    blue = ROUND (blue_fm * MAX_GUN);

	    /*
	     * Do color mapping (colormask option) Leave color 0 alone.
	     */
	    if (col_tab_no != 0)
	    {
		red_mask = red;
		green_mask = green;
		blue_mask = blue;

/* Background color */
		red_0 = color_set[0][_RED];
		green_0 = color_set[0][_GREEN];
		blue_0 = color_set[0][_BLUE];

/* Opposite of background color */
		red_7 = MAX_GUN - color_set[0][_RED];
		green_7 = MAX_GUN - color_set[0][_GREEN];
		blue_7 = MAX_GUN - color_set[0][_BLUE];

		if (red == red_7 && green == green_7 && blue == blue_7)
		{
		    if (colormask[3] == NO)
		    {
			red_mask = red_0;
			green_mask = green_0;
			blue_mask = blue_0;
		    }
		}
		else
		{
		    if (colormask[4] == YES)
		    {
			if (colormask[0] == NO)
			    red_mask = red_0;
			if (colormask[1] == NO)
			    green_mask = green_0;
			if (colormask[2] == NO)
			    blue_mask = blue_0;
		    }
		    else
		    {
			if (colormask[0] == YES)
			{
			    green_mask = red_mask;
			    blue_mask = red_mask;
			}
			if (colormask[1] == YES)
			{
			    red_mask = green_mask;
			    blue_mask = green_mask;
			}
			if (colormask[2] == YES)
			{
			    red_mask = blue_mask;
			    green_mask = blue_mask;
			}
		    }
		}

		red = red_mask;
		green = green_mask;
		blue = blue_mask;
	    }

	    /*
	     * Process the new color table value
	     */
/*       fprintf(stderr,"IN DO a \n");*/
	    if (col_tab_no < num_col)
	    {
		/*
		 * A settable color.
		 */
		dev.attributes (SET_COLOR_TABLE, col_tab_no,
				red, green, blue);
		color_set[col_tab_no][STATUS] = SET;
	    }
	    else if (col_tab_no > 7)
	    {
		color_set[col_tab_no][STATUS] = MAPPED;
	    }
	    else
	    {
		if (col_tab_no > 0)
		{
		    color_set[col_tab_no][STATUS] = MAP_SET;
		    /*
		     * Means that we need to map these but they are still set
		     * to the original colors and can't be changed. Color 0
		     * is always background, and any other color exactly the
		     * same color as color 0 "wants to be", even if color 0
		     * actually can't be and hasn't been reset, is mapped to
		     * color zero if it can't be set.
		     */
		}
/*       fprintf(stderr,"IN DO b \n");*/
	    }
	    /*
	     * Save the color that this color table number wants to be
	     */
	    color_set[col_tab_no][_RED] = red;
	    color_set[col_tab_no][_GREEN] = green;
	    color_set[col_tab_no][_BLUE] = blue;
	    /*
	     * grey level is "Black and White TV style" mapped from color,
	     * with corrections added by the subroutine greycorr.
	     */
	    grey = (blue * 1 + red * 2 + green * 4 + 6) / 7;
	    color_set[col_tab_no][_GREY] = greycorr (grey);
	    /*
	     * If the next command is also a set color table command, we can
	     * postpone doing this and kill 2 (or more) birds with one stone.
	     */
	    c = getc (pltin);
	    if (c == EOF) end_of_file();
	    ungetc ((char) c, pltin);
	    if (c != VP_SET_COLOR_TABLE)
	    {
		if (mono)
		{
		    for (ii = 1; ii <= MAX_COL; ii++)
		    {
			if (color_set[ii][_RED] == color_set[0][_RED] &&
			    color_set[ii][_GREEN] == color_set[0][_GREEN] &&
			    color_set[ii][_BLUE] == color_set[0][_BLUE])
			{
			    color_set[ii][MAP] = 0;
			}
			else
			{
			    color_set[ii][MAP] = 7;
			}
		    }
		}
		else
		{
		    /*
		     * For all color table entries that aren't set, (because
		     * we ran out of colors, for example) find the best color
		     * that IS set and use that instead.
		     */
/*       fprintf(stderr,"IN DO 1 \n");*/
		    for (ii = num_col; ii <= MAX_COL; ii++)
		    {
			if (color_set[ii][STATUS] & MAPPED)
			{
			    min_dist = MAX_GUN * MAX_GUN * 8;
			    for (i = num_col_8 - 1; i >= 0; i--)
			    {
				/*
				 * Colors 1 through 7 are guaranteed SET, So
				 * we always get an answer. Color zero is
				 * background and special. To map to it you
				 * have to hit its color exactly, and no
				 * other color matched exactly first.
				 */
				if (color_set[i][STATUS] & SET)
				{
				    if (color_set[i][STATUS] == SET)
				    {
					k = color_set[i][_RED] - color_set[ii][_RED];
					dist = 2 * k * k;
					k = color_set[i][_GREEN] - color_set[ii][_GREEN];
					dist += 4 * k * k;
					k = color_set[i][_BLUE] - color_set[ii][_BLUE];
					dist += k * k;
				    }
				    else
				    {
					k = MAX_GUN * ((i & 2) / 2) - color_set[ii][_RED];
					dist = 2 * k * k;
					k = MAX_GUN * ((i & 4) / 4) - color_set[ii][_GREEN];
					dist += 4 * k * k;
					k = MAX_GUN * ((i & 1) / 1) - color_set[ii][_BLUE];
					dist += k * k;
				    }
				    if (dist < min_dist && (i != 0 || dist == 0))
				    {
					min_dist = dist;
					best_col = i;
					if (dist == 0)
					{
					    /*
					     * Might as well look no further
					     */
					    break;
					}
				    }
				}
			    }
			    color_set[ii][MAP] = best_col;
/*       fprintf(stderr,"IN MAP %d -%d \n",ii,best_col);*/
			}
		    }
		}
	    }
}
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
 *  source file:   ./filters/init_vplot.c
 *
 * Joe Dellinger (SEP), Feb 18 1988
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger, Feb 20 1988
 *	#define GETPAR to be fetch if SEPlib version.
 * Joe Dellinger Feb 24 1988
 * 	Moved gen_do_dovplot to genlib, where it belongs.
 *	txfont, txprec, txovly, overlay defaults remembered.
 * Joe Dellinger Feb 28 1988
 *	Changed fatness to scale with size of screen and scale,
 *	like any regular geometric attribute.
 * Joe Dellinger Mar 4 1988
 *	Fixed slight bug with fatness scaling and style=old.
 * Joe Dellinger Mar 17 1988
 *	Added "colormask" option.
 * Joe Dellinger April 15 1988
 *	Check for illegal colormask values.
 * Joe Dellinger Oct 28 1988
 *	Put in "no_stretch_text" flag.
 * Steve Cole Feb 2 1989
 *      Added default x and y scale parameters for use by vppen.
 * Joe Dellinger May 7 1989
 *	Added break=ignore option.
 *	Made txsquare=yes the default.
 * Joe Dellinger May 20 1989
 *	"/dev/tty" may NOT be openable. Behave gracefully in this case.
 * Joe Dellinger August 2 1989
 *	Put in redpow, greenpow, bluepow, redmap, greenmap, bluemap
 * W. Bauske IBM 03-26-91 
 *	Apply the SysV fixes to the code for the RS/6000
 *	Removed re-declare of system alloc routines for RS/6000
 * Dave Nichols 9-17-92
 *      Include stdlib.h if it is available instead of defining malloc.  
 * Stew Levin  2/27/95 
 *	Solaris mods.
 * Bob Clapp 1/24/96
 *      Changed a whole bunch of fetches to getches
 * Bob Clapp 10/98 
 *      Changed signals from bsd to posix1 for linux
 */
#include<sitedef.h>
#include	<stdio.h>
#include	<math.h>
#if defined(HAVE_TERMIO_H)
#include	<termio.h>
#endif
#if defined(SYS_IOCTL_H)
#include	<sys/ioctl.h>
#endif
#if defined(HAVE_SGTTY_H)
#include	<sgtty.h>
#endif
#if  defined(HAVE_SYS_TYPES_H)
#include	<sys/types.h>
#endif

#if  defined(HAVE_SYS_STAT_H)
#include	<sys/stat.h>
#endif
#include	<ctype.h>
#include	<string.h>

#include	"./include/vplot.h"

#include	"./include/params.h"	/* for machine dependencies */
#include	"./include/enum.h"
#include	"./include/err.h"
#include	"./include/attrcom.h"
#include	"./include/intcom.h"
#include	"./include/mesgcom.h"
#include	"./include/erasecom.h"
#include	"./include/closestat.h"
#include	"./include/pat.h"
#include	"./include/vertex.h"
#include	"./include/round.h"
#include	"./include/extern.h"

#define		OPEN_ERROR	-1

#ifdef SEP
#define		GETPAR	fetch
#else
#define		GETPAR	getpar
#endif /* SEP */

#if defined(HAVE_TERMIO_H)
struct termio	tty_clean_state;
struct termio	tty_plot_state;
#else
struct sgttyb   tty_clean_state;
struct sgttyb   tty_plot_state;
int             tty_clean_local_mode;
int             tty_plot_local_mode;
#endif /* USG */

/* 
 * The following variables must ALWAYS
 * be set in the device open or device reset file.
 * All their defaults are set to absurd values so that
 * we'll crash quickly it they're not set as required.
 */
/* Screen dimensions */
int             dev_xmax = 0, dev_ymax = 0, dev_xmin = 0, dev_ymin = 0;
/* Number of pixels per inch in the horizontal (X) direction on the device */
float           pixels_per_inch = 0.;
/* vertical height (on screen) / horizontal width for a single pixel */
float           aspect_ratio = 0.;
/*
 * Number of SETTABLE COLORS that the device has.
 * Non settable colors don't count.
 */
int             num_col = -1;

/*
 * Other things that may need to be reset in dev.open
 * (Can't reset them some of them in dev.reset, because the user may
 * override from the command line.)
 */
/* Does device need erase at end? */
int             need_end_erase = NO;
/* should the output be buffered? */
int             buffer_output = YES;
/* should the input be buffered? */
int             buffer_input = YES;
/* should pipes be allowed for input? */
int             allow_pipe = YES;
/* Can the device do its own clipping? (of vectors and polygons.) */
int             smart_clip = NO;
/* Can the device stretch AND clip its own raster? */
int             smart_raster = NO;

/*
 * These may be reset to device-dependent defaults.
 */
float           fatmult = 1.;
float           patternmult = 1.;

/*
 * Look in extern.h for a better-organized list of the preceding.
 * Here I have left things which are also filter options with the
 * other filter options, even though they also belong in the above
 * categories.
 */

/*
 * flags and variables
 */
char            wstype[25];
char            callname[25];
int             xcenterflag, ycenterflag;
int             ever_called = NO;
int             first_time = YES;
int             device_open = NO;	/* used by ERR */
int             out_isatty = YES;	/* If NO, output is a file or pipe */
int             nplots = 0;	/* number of plots made */
int             no_stretch_text = YES;	/* Don't allow stretched text? */

/*
 * coordinate variables
 */
int             xwmax, xwmin, ywmax, ywmin;	/* window */
int             xnew, ynew;	/* new pen location */
int             xold, yold;	/* old pen location */
int             xorigin = 0, yorigin = 0;	/* global "origin" */
char           *txbuffer;
int             txbuflen, vxbuflen;
struct vertex  *vxbuffer;
int             xret, yret;
float           mxx, mxy, myx, myy;
/*
 * This is so the device can throw in a coordinate transformation
 * for its convenience ON TOP OF whatever transformation the user
 * wants.
 */
int             default_rotate = 0, default_hshift = 0, default_vshift = 0;
float           default_scale = 1., default_xscale = 1., default_yscale = 1.;

/*
 * attribute variables
 */
int             linestyle;
int             cur_color = DEFAULT_COLOR;
int             pat_color = DEFAULT_COLOR + 1;
int             next_color;
struct txalign  txalign;

int             txfont = DEFAULT_FONT;
int             txprec = DEFAULT_PREC;
int             txovly = OVLY_NORMAL;
int             default_txfont, default_txprec, default_txovly;
int             fat = 0;
int             afat = 0;
int             ipat = 0;	/* currently loaded pattern */
struct pat      pat[NPAT + 1];
float           dashsum = 0.;
int             dashon = NO;	/* Dashed lines? */
float           dashes[MAXDASH * 2];
float           dashpos = 0.;	/* Position in current dashing pattern */
int             color_set[MAX_COL + 1][_NUM_PRIM];
extern int      greycorr ();
int             num_col_8;


/*
 * filter options - flags
 */
/* Monochrome device? */
int             mono = NO;
/*
 * Invras determines the polarity of dithered raster.
 * invras=n means raster works the same way as vectors; what is normally
 * WHITE on most devices is BLACK on a hardcopy black and white device.
 * invras=y is the proper default to make dithered images not come out as negatives
 */
int             invras = YES;
int             window = YES;
int             shade = YES;
int             brake = BREAK_BREAK;
int             framewindows = NO;
int             endpause = NO;
int             cachepipe = NO;
/* Setting cachepipe = YES will copy any piped files to a temporary file,
 * this may get done in dev.open.
 * This is useful if the program may want to reverse seek. e.g. Xvpen, Sunpen
 */
int             allowecho = NEVER_SET;
/*
 * setting allowecho NO (which may get done
 * in dev.open) means that the output device
 * is the user's terminal and echoing chars
 * back at the user would insert nonsense
 * into the plot stream.  In this case we
 * explicitly shut down echoing and output
 * translation until the end of the job.
 * Arguably, output translation should be
 * set/unset in dev.open or dev.reset as we may
 * be plotting on a terminal not our own or
 * the device filter may have built in
 * workarounds for known output translation
 * effects. I use it here as a safety measure.
 */
int             wantras = YES;
int             colormask[5] = {1, 1, 1, 1, 1};

/*
 * filter options - enumerated
 */
int             style = NO_STYLE_YET;
int             default_style = STYLE;
int             rotate;
int             size = RELATIVE;
int             erase = FORCE_INITIAL | DO_LITERALS;

/*
 * filter options - valued
 */
int             xcenter, ycenter;	/* Vplot of point to force as center */
int             fatbase = 0;
int             epause = 0;	/* time to pause before erasing screen */
int             overlay = 0;	/* 1=overlay 0=replace */
int             default_overlay;

/*
 * 0 = none
 * 1 = random
 * 2 = Ordered
 * 3 = Floyd-Steinberg
 * 4 = Halftone
 */
int             dither = 1;	/* Dithering type */

int             xWmax, xWmin, yWmax, yWmin;
float           hshift, vshift;	/* Allow global translation of plot */
float           scale;	/* global scale */
float           xscale;	/* global x-scale */
float           yscale;	/* global y-scale */
float           hdevscale;	/* Vplot units to device units for x */
float           vdevscale;	/* Vplot units to device units for y */
float           txscale;/* global text scale */
float           mkscale;/* global marker scale */
float           dashscale;	/* global dashed line scale */
char            interact[MAXFLEN + 1] = "";	/* Where to store coordinate
						 * file */
float           greyc = 1.;	/* Nonlinear correction */
float           pixc = 1.;	/* Pixel overlap correction */
float           redpow = 1., redmap[4] = {1., 0., 0., 0.};
float           greenpow = 1., greenmap[4] = {0., 1., 0., 0.};
float           bluepow = 1., bluemap[4] = {0., 0., 1., 0.};

/* filter options - resettable between plots */
int             user_rotate;
float           user_txscale;
float           user_mkscale;
float           user_dashscale;
float           user_scale;
float           user_xscale;
float           user_yscale;
int             user_size;
float           user_hshift;
float           user_vshift;
int             user_xwmax_flag;
float           user_xwmax;
int             user_ywmax_flag;
float           user_ywmax;
int             user_xwmin_flag;
float           user_xwmin;
int             user_ywmin_flag;
float           user_ywmin;

int             ifat = 0;
float           fatmult_orig;

/*
 * file and terminal control variables
 */
int             pltoutfd, stderrfd, controlfd;
extern int      genmessage ();
int             (*message) () = genmessage;
struct stat     stderrstat;
struct stat     pltoutstat;
FILE           *pltout, *pltin;
FILE           *fopen ();
FILE           *fdopen ();
FILE           *controltty;
char            outbuf[BUFSIZ];
char           *getenv ();
#if defined(HAVE_STDLIB_H)
#include <stdlib.h>
#else
char           *malloc ();
char           *realloc ();
#endif
char            group_name[MAXFLEN + 1];
int             group_number = 0;
FILE           *pltinarray[MAXIN];
char            pltinname[MAXIN][MAXFLEN + 1];
char            pltname[MAXFLEN + 1] = "";
int             infileno = 0;
extern int      gen_do_dovplot ();
int             (*genreader) () = gen_do_dovplot;

/*
 * Initialize and declare global variables.
 * Use getpar and getenv to search for command-line options.
 */

init_vplot ()
{
char           *stringptr;
int             ii;
char            string[MAXFLEN + 1];
float           ftemp;


    txbuffer = (char *) malloc (TXBUFLEN);
    txbuflen = TXBUFLEN;
/*
 * Sun III lint complains about "pointer alignment problem" here,
 * although the documentation makes it sound like that shouldn't
 * be possible.
 */
    vxbuffer = (struct vertex *) malloc (sizeof (struct vertex) * VXBUFLEN);
    vxbuflen = VXBUFLEN;

    /*
     * If the device can't do color, it can set mono to YES If device is
     * mono, this overrides the user's insistence that it isn't. So GETPAR
     * before we call dev.open. 
     */
    GETPAR ("mono", "1", &mono);

    strcpy (wstype, "default");
    if ((stringptr = getenv ("WSTYPE")) != NULL)
	strcpy (wstype, stringptr);
    /*
     * Don't want to accidentally get this from a history file. Must come for
     * the command line, so getpar, not GETPAR 
     */
    getpar ("wstype", "s", wstype);

    /*
     * Initialize all patterns to be undefined. (device can create a
     * device-dependent default set if it wishes) 
     */
    for (ii = 0; ii <= NPAT; ii++)
    {
	pat[ii].patbits = NULL;
    }

/*
 * The device-independent code MAY NOT actually read from
 * the terminal itself. However, these variables are 
 * declared so that the device-dependent routines can use
 * them if they wish as a courtesy.
 */
    controlfd = open ("/dev/tty", 0);
    if (controlfd == OPEN_ERROR)
    {
	controltty = NULL;
    }
    else
    {
	controltty = fdopen (controlfd, "r");
    }

    /*
     * Call device open before doing anything else. this finalizes pltout 
     */
    dev.open ();
    device_open = YES;

/*
 * Beware trying to print out error messages until
 * we've figured out where to connect the "message" routine!
 */

    if (buffer_output)
	setbuf (pltout, outbuf);
    else
	setbuf (pltout, (char *) NULL);

    /*
     * If graphics output going to control terminal, disable echoing and
     * arrange for use of device's message routine 
     */
    pltoutfd = fileno (pltout);
    stderrfd = fileno (stderr);
    out_isatty = isatty (pltoutfd);

    if (allowecho == NEVER_SET)
    {
	allowecho = YES;
	if (out_isatty)
	{
	    fstat (pltoutfd, &pltoutstat);
	    fstat (stderrfd, &stderrstat);
#ifdef SEP
	    if ((pltoutstat.st_dev == stderrstat.st_dev))
		/*
		 * Something in seplib makes the next 2 tests not work, I
		 * don't know why. Unfortunately, this means the SEP versions
		 * aren't as good about detecting the case of sending the
		 * output to a different terminal than the one you're logged
		 * in at. We should probably get this fixed sometime. - Joe
		 * D. 
		 */
#else
	    if ((pltoutstat.st_dev == stderrstat.st_dev) &&
		(pltoutstat.st_ino == stderrstat.st_ino) &&
		(pltoutstat.st_rdev == stderrstat.st_rdev))
#endif /* SEP */
	    {
		allowecho = NO;
		message = dev.message;
	    }
	}
    }

    /*
     * process YES or NO option arguments. GETPAR means it could have come
     * from a seplib history file, getpar means it had to come from the
     * command line. Use "getpar" when having the right value of the
     * parameter is critical because setting the value overrides frontend's
     * judgement. 
     */
    getpar ("echo", "1", &allowecho);
    if (!allowecho)
    {
#if defined(HAVE_TERMIO_H)
	if (ioctl (pltoutfd, TCGETA, &tty_clean_state) == -1)
	{
		ERR (FATAL, name, "Bad ioctl call!");
	}
	tty_plot_state = tty_clean_state;
	tty_plot_state.c_lflag &= ~ECHO;
	if (ioctl (pltoutfd, TCSETAW, &tty_plot_state) == -1)
	{
		ERR (FATAL, name, "Bad ioctl call! (2)");
	}
#else /* USG */
	ioctl (pltoutfd, TIOCGETP, (char *) (&tty_clean_state));
	bcopy ((char *) (&tty_clean_state), (char *) (&tty_plot_state),
	       sizeof (struct sgttyb));
	ioctl (pltoutfd, TIOCLGET, (char *) (&tty_clean_local_mode));
	tty_plot_local_mode = tty_clean_local_mode | LLITOUT;
	ioctl (pltoutfd, TIOCLSET, (char *) (&tty_plot_local_mode));
	tty_plot_state.sg_flags &= (~ECHO);
	tty_plot_state.sg_flags &= (~LCASE);
	tty_plot_state.sg_flags |= (CBREAK);
	ioctl (pltoutfd, TIOCSETN, (char *) (&tty_plot_state));
#endif /* USG */
    }
    getpar ("endpause", "1", &endpause);
    GETPAR ("cachepipe", "1", &cachepipe);
    GETPAR ("shade", "1", &shade);
    GETPAR ("wantras", "1", &wantras);
    getch ("window", "1", &window);
    getch ("frame", "1", &framewindows);
    GETPAR ("overlay", "1", &overlay);
    default_overlay = overlay;
    GETPAR ("invras", "1", &invras);

    GETPAR ("txsquare", "1", &no_stretch_text);

    GETPAR ("colormask", "1", colormask);
    if (colormask[4] == NO)
    {
	if (colormask[0] + colormask[1] + colormask[2] + colormask[3] != 1)
	    ERR (FATAL, name, "Invalid colormask option.");
    }
    if (colormask[4] == NO && colormask[3] == YES)
    {
	colormask[4] = YES;
    }

    GETPAR ("redpow", "f", &redpow);
    GETPAR ("greenpow", "f", &greenpow);
    GETPAR ("bluepow", "f", &bluepow);

    if (redpow <= 0. || greenpow <= 0. || bluepow <= 0.)
	ERR (FATAL, name, "Invalid color pow option.");

    GETPAR ("red", "f", redmap);
    GETPAR ("green", "f", greenmap);
    GETPAR ("blue", "f", bluemap);

/*
 * Valued arguments
 */

    getpar ("dither", "d", &dither);
    getpar ("greyc", "f", &greyc);
    getpar ("pixc", "f", &pixc);

    GETPAR ("txfont", "d", &txfont);
    GETPAR ("txprec", "d", &txprec);
    GETPAR ("txovly", "d", &txovly);
    default_txfont = txfont;
    default_txprec = txprec;
    default_txovly = txovly;

    if (getch ("erase", "s", string))
    {
	if ((string[0] == 'n') || (string[0] == 'N'))
	    erase = NO;
	else
	if ((string[0] == 'o') || (string[0] == 'O'))
	    erase = FORCE_INITIAL;
	else
	if ((string[0] == 'l') || (string[0] == 'L'))
	    erase = DO_LITERALS;
	else
	    erase = FORCE_INITIAL | DO_LITERALS;
    }

    if (getch ("break", "s", string))
    {
	if (string[0] == 'b')
	    brake = BREAK_BREAK;
	else
	if (string[0] == 'e')
	    brake = BREAK_ERASE;
	else
	    brake = BREAK_IGNORE;
    }

    xcenter = 0;
    ycenter = 0;
    xcenterflag = NO;
    ycenterflag = NO;
    if (GETPAR ("xcenter", "f", &ftemp))
    {
	xcenterflag = YES;
	xcenter = ROUND (ftemp * RPERIN);
    }
    if (GETPAR ("ycenter", "f", &ftemp))
    {
	ycenterflag = YES;
	ycenter = ROUND (ftemp * RPERIN);
    }

    if ((stringptr = getenv ("PATTERNMULT")) != NULL)
    {
	sscanf (stringptr, "%f", &patternmult);
    }
    GETPAR ("patternmult", "f", &patternmult);

    GETPAR ("interact", "s", interact);
    GETPAR ("pause", "d", &epause);

    if (interact[0] != '\0')
    {
	epause = 0;		/* interact makes it own sort of pausing */
	endpause = NO;
    }

    /*
     * Find the default style 
     */
    stringptr = NULL;
    if (GETPAR ("style", "s", string))
	stringptr = string;
    else
	stringptr = getenv ("PLOTSTYLE");
    if (stringptr != NULL)
    {
	if ((stringptr[0] == 'r') || (stringptr[0] == 'R') ||
	    (stringptr[0] == 'm') || (stringptr[0] == 'M'))
	    default_style = ROTATED;
	else
	if ((stringptr[0] == 'o') || (stringptr[0] == 'O'))
	    default_style = OLD;
	else
	if ((stringptr[0] == 'a') || (stringptr[0] == 'A'))
	    default_style = ABSOLUTE;
	else
	    default_style = STANDARD;
    }

/*
 * Options changeable between calls to dovplot
 * (By calling reset_parameters or reset())
 * Dovplot calls reset every time it starts. It only calls
 * reset_parameters the first time.
 *
 * Things in this category:
 * user_*
 * fatmult_orig, fatbase
 */

    if ((stringptr = getenv ("FATMULT")) != NULL)
    {
	sscanf (stringptr, "%f", &fatmult);
    }
    GETPAR ("fatmult", "f", &fatmult);
    fatmult_orig = fatmult;

    user_rotate = 0;
    getch ("rotate", "d", &user_rotate);

    user_txscale = 1.0;
    user_mkscale = 1.0;
    user_dashscale = 1.0;
    GETPAR ("txscale", "f", &user_txscale);
    GETPAR ("mkscale", "f", &user_mkscale);
    GETPAR ("dashscale", "f", &user_dashscale);

    user_scale = 1.0;
    user_xscale = 1.0;
    user_yscale = 1.0;
    getch ("scale", "f", &user_scale);
    getch ("xscale", "f", &user_xscale);
    getch ("yscale", "f", &user_yscale);

    user_size = size;
    if (getch ("size", "s", string))
    {
	if ((string[0] == 'a') || (string[0] == 'A'))
	    user_size = ABSOLUTE;
	else
	    user_size = RELATIVE;
    }

    user_hshift = 0.;
    getch ("hshift xshift", "f", &user_hshift);
    user_vshift = 0.;
    getch ("vshift yshift", "f", &user_vshift);

    user_xwmax_flag = GETPAR ("xwmax", "f", &user_xwmax);
    user_ywmax_flag = GETPAR ("ywmax", "f", &user_ywmax);
    user_xwmin_flag = GETPAR ("xwmin", "f", &user_xwmin);
    user_ywmin_flag = GETPAR ("ywmin", "f", &user_ywmin);

    GETPAR ("fat", "d", &fatbase);



/*
 *  These parameters can simply be changed at any time with no
 *  need for re-initialization of anything else:
 *
 *  shade, wantras
 *
 */

}


setstyle (new_style)
    int             new_style;
{
    /*
     * Check to see if the style has changed 
     */
    if (new_style == style)
	return;

    style = new_style;
    reset_parameters ();
}

reset_parameters ()
{
float           inches;	/* scaling base for y axis */
float           screenheight, screenwidth;	/* true size of the screen */
int             ix, iy;

    xorigin = 0;
    yorigin = 0;

    rotate = default_rotate;
    rotate += user_rotate;

    txscale = user_txscale;
    mkscale = user_mkscale;
    dashscale = user_dashscale;

    scale = default_scale * user_scale;
    xscale = default_xscale * user_xscale;
    yscale = default_yscale * user_yscale;

    fatmult = fatmult_orig;

    size = user_size;

    switch (style)
    {
/*
 * The old standard on the machine erebus.
 * Pretty much dead now. Still useful for some old programs nobody's
 * wanted to update.
 */
    case OLD:
	txscale /= 3.;
	fatmult /= 3.;
	scale *= 3.;
	inches = STANDARD_HEIGHT;
	break;
/*
 * The old standard on the machine mazama. A useful coordinate system
 * for geological sorts of plots. 
 * The Y axis goes the long way, across the screen, which is the device's
 * horizontal axis.
 */
    case ROTATED:
	rotate += 90;
	inches = ROTATED_HEIGHT;
	break;
    case ABSOLUTE:
	size = ABSOLUTE;
    case STANDARD:
    default:
	inches = STANDARD_HEIGHT;
	break;
    }

    if (rotate >= 0)
	rotate = rotate % 360;
    else
	rotate = ((rotate % 360) + 360) % 360;

    mxx = cos (2. * 3.14159 * rotate / 360.);
    myy = cos (2. * 3.14159 * rotate / 360.);
    mxy = sin (2. * 3.14159 * rotate / 360.);
    myx = -sin (2. * 3.14159 * rotate / 360.);

    if (size == ABSOLUTE)
    {
	vdevscale = pixels_per_inch / (float) (RPERIN * aspect_ratio);
	hdevscale = pixels_per_inch / (float) RPERIN;
    }
    else
    {
	/*
	 * Fit the inches x inches unit square into a displayable box with
	 * aspect ratio SCREEN_RATIO 
	 */
	screenwidth = (dev_xmax - dev_xmin) * pixels_per_inch;
	screenheight =
	 (dev_ymax - dev_ymin) * pixels_per_inch * aspect_ratio;
	if ((screenheight / screenwidth) > SCREEN_RATIO)
	{
	    vdevscale = (SCREEN_RATIO * ((dev_xmax - dev_xmin) / aspect_ratio)) /
	     (inches * RPERIN);
	    hdevscale = vdevscale * aspect_ratio;
	}
	else
	{
	    vdevscale = (dev_ymax - dev_ymin) / (inches * RPERIN);
	    hdevscale = vdevscale * aspect_ratio;
	}
    }

    hshift = default_hshift;
    vshift = default_vshift;

    if (style == ROTATED)
    {
	vshift += dev_ymax - dev_ymin;
    }

    yscale *= scale;
    xscale *= scale;
    mkscale *= scale;

/*
 * Set up fatness multiplication factor
 */
    fatmult *= scale * hdevscale * RPERIN / FATPERIN;

    /*
     * The point (xcenter,ycenter) in vplot coordinates is to be centered in
     * the screen. 
     */
    if (xcenterflag || ycenterflag)
    {
	vptodevxy (xcenter, ycenter, &ix, &iy);
	if (xcenterflag)
	{
	    hshift += (dev_xmax + dev_xmin) / 2 - ix;
	}
	if (ycenterflag)
	{
	    vshift += (dev_ymax + dev_ymin) / 2 - iy;
	}
    }

    hshift += user_hshift * pixels_per_inch;
    vshift += user_vshift * pixels_per_inch / aspect_ratio;

/* plot window parameters defaulted */
/* to maximum size */

    devtovpw (dev_xmin, dev_ymin, dev_xmax, dev_ymax,
	      &xWmin, &yWmin, &xWmax, &yWmax);

    if (user_xwmax_flag)
	xWmax = ROUND (user_xwmax * RPERIN);
    if (user_ywmax_flag)
	yWmax = ROUND (user_ywmax * RPERIN);
    if (user_xwmin_flag)
	xWmin = ROUND (user_xwmin * RPERIN);
    if (user_ywmin_flag)
	yWmin = ROUND (user_ywmin * RPERIN);

    vptodevw (xWmin, yWmin, xWmax, yWmax, &xWmin, &yWmin, &xWmax, &yWmax);

    wlimit (dev_xmin, dev_xmax, &xWmin, &xWmax);
    wlimit (dev_ymin, dev_ymax, &yWmin, &yWmax);

    xwmax = xWmax;		/* plot window parameters defaulted */
    xwmin = xWmin;		/* to maximum size 		 */
    ywmax = yWmax;
    ywmin = yWmin;
    reset_windows ();
}

add_a_cor (filename, xcor, ycor)
    char           *filename;
    int             xcor, ycor;
{
static int      first_time = YES;
static FILE    *outfp;

    if (first_time == YES)
    {
	outfp = fopen (filename, "w");
	if (outfp == NULL)
	{
	    ERR (FATAL, name, "Can't open interact output file %s!", filename);
	}
	first_time = NO;
    }
    fprintf (outfp, "%f\t%f\n", (float) xcor / RPERIN, (float) ycor / RPERIN);
}

wlimit (min, max, xmin, xmax)
    int             min, max;
    int            *xmin, *xmax;
{
    if (*xmin < min)
	*xmin = min;

    if (*xmax > max)
	*xmax = max;
}
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
 *  source file:   ./filters/init_vplot.c
 *
 * Joe Dellinger (SEP), Feb 18 1988
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger, Feb 20 1988
 *	#define GETPAR to be fetch if SEPlib version.
 * Joe Dellinger Feb 24 1988
 * 	Moved gen_do_dovplot to genlib, where it belongs.
 *	txfont, txprec, txovly, overlay defaults remembered.
 * Joe Dellinger Feb 28 1988
 *	Changed fatness to scale with size of screen and scale,
 *	like any regular geometric attribute.
 * Joe Dellinger Mar 4 1988
 *	Fixed slight bug with fatness scaling and style=old.
 * Joe Dellinger Mar 17 1988
 *	Added "colormask" option.
 * Joe Dellinger April 15 1988
 *	Check for illegal colormask values.
 * Joe Dellinger Oct 28 1988
 *	Put in "no_stretch_text" flag.
 * Steve Cole Feb 2 1989
 *      Added default x and y scale parameters for use by vppen.
 * Joe Dellinger May 7 1989
 *	Added break=ignore option.
 *	Made txsquare=yes the default.
 * Joe Dellinger May 20 1989
 *	"/dev/tty" may NOT be openable. Behave gracefully in this case.
 * Joe Dellinger August 2 1989
 *	Put in redpow, greenpow, bluepow, redmap, greenmap, bluemap
 * W. Bauske IBM 03-26-91 
 *	Apply the SysV fixes to the code for the RS/6000
 *	Removed re-declare of system alloc routines for RS/6000
 * Dave Nichols 9-17-92
 *      Include stdlib.h if it is available instead of defining malloc.  
 * Stew Levin  2/27/95 
 *	Solaris mods.
 * Bob Clapp 1/24/96
 *      Changed a whole bunch of fetches to getches
 * Bob Clapp 10/98 
 *      Changed signals from bsd to posix1 for linux
 */
#include<sitedef.h>
#include	<stdio.h>
#include	<math.h>
#if defined(HAVE_TERMIO_H)
#include	<termio.h>
#endif
#if defined(SYS_IOCTL_H)
#include	<sys/ioctl.h>
#endif
#if defined(HAVE_SGTTY_H)
#include	<sgtty.h>
#endif
#if  defined(HAVE_SYS_TYPES_H)
#include	<sys/types.h>
#endif

#if  defined(HAVE_SYS_STAT_H)
#include	<sys/stat.h>
#endif
#include	<ctype.h>
#include	<string.h>

#include	"./include/vplot.h"

#include	"./include/params.h"	/* for machine dependencies */
#include	"./include/enum.h"
#include	"./include/err.h"
#include	"./include/attrcom.h"
#include	"./include/intcom.h"
#include	"./include/mesgcom.h"
#include	"./include/erasecom.h"
#include	"./include/closestat.h"
#include	"./include/pat.h"
#include	"./include/vertex.h"
#include	"./include/round.h"
#include	"./include/extern.h"

#define		OPEN_ERROR	-1

#ifdef SEP
#define		GETPAR	fetch
#else
#define		GETPAR	getpar
#endif /* SEP */

#if defined(HAVE_TERMIO_H)
struct termio	tty_clean_state;
struct termio	tty_plot_state;
#else
struct sgttyb   tty_clean_state;
struct sgttyb   tty_plot_state;
int             tty_clean_local_mode;
int             tty_plot_local_mode;
#endif /* USG */

/* 
 * The following variables must ALWAYS
 * be set in the device open or device reset file.
 * All their defaults are set to absurd values so that
 * we'll crash quickly it they're not set as required.
 */
/* Screen dimensions */
int             dev_xmax = 0, dev_ymax = 0, dev_xmin = 0, dev_ymin = 0;
/* Number of pixels per inch in the horizontal (X) direction on the device */
float           pixels_per_inch = 0.;
/* vertical height (on screen) / horizontal width for a single pixel */
float           aspect_ratio = 0.;
/*
 * Number of SETTABLE COLORS that the device has.
 * Non settable colors don't count.
 */
int             num_col = -1;

/*
 * Other things that may need to be reset in dev.open
 * (Can't reset them some of them in dev.reset, because the user may
 * override from the command line.)
 */
/* Does device need erase at end? */
int             need_end_erase = NO;
/* should the output be buffered? */
int             buffer_output = YES;
/* should the input be buffered? */
int             buffer_input = YES;
/* should pipes be allowed for input? */
int             allow_pipe = YES;
/* Can the device do its own clipping? (of vectors and polygons.) */
int             smart_clip = NO;
/* Can the device stretch AND clip its own raster? */
int             smart_raster = NO;

/*
 * These may be reset to device-dependent defaults.
 */
float           fatmult = 1.;
float           patternmult = 1.;

/*
 * Look in extern.h for a better-organized list of the preceding.
 * Here I have left things which are also filter options with the
 * other filter options, even though they also belong in the above
 * categories.
 */

/*
 * flags and variables
 */
char            wstype[25];
char            callname[25];
int             xcenterflag, ycenterflag;
int             ever_called = NO;
int             first_time = YES;
int             device_open = NO;	/* used by ERR */
int             out_isatty = YES;	/* If NO, output is a file or pipe */
int             nplots = 0;	/* number of plots made */
int             no_stretch_text = YES;	/* Don't allow stretched text? */

/*
 * coordinate variables
 */
int             xwmax, xwmin, ywmax, ywmin;	/* window */
int             xnew, ynew;	/* new pen location */
int             xold, yold;	/* old pen location */
int             xorigin = 0, yorigin = 0;	/* global "origin" */
char           *txbuffer;
int             txbuflen, vxbuflen;
struct vertex  *vxbuffer;
int             xret, yret;
float           mxx, mxy, myx, myy;
/*
 * This is so the device can throw in a coordinate transformation
 * for its convenience ON TOP OF whatever transformation the user
 * wants.
 */
int             default_rotate = 0, default_hshift = 0, default_vshift = 0;
float           default_scale = 1., default_xscale = 1., default_yscale = 1.;

/*
 * attribute variables
 */
int             linestyle;
int             cur_color = DEFAULT_COLOR;
int             pat_color = DEFAULT_COLOR + 1;
int             next_color;
struct txalign  txalign;

int             txfont = DEFAULT_FONT;
int             txprec = DEFAULT_PREC;
int             txovly = OVLY_NORMAL;
int             default_txfont, default_txprec, default_txovly;
int             fat = 0;
int             afat = 0;
int             ipat = 0;	/* currently loaded pattern */
struct pat      pat[NPAT + 1];
float           dashsum = 0.;
int             dashon = NO;	/* Dashed lines? */
float           dashes[MAXDASH * 2];
float           dashpos = 0.;	/* Position in current dashing pattern */
int             color_set[MAX_COL + 1][_NUM_PRIM];
extern int      greycorr ();
int             num_col_8;


/*
 * filter options - flags
 */
/* Monochrome device? */
int             mono = NO;
/*
 * Invras determines the polarity of dithered raster.
 * invras=n means raster works the same way as vectors; what is normally
 * WHITE on most devices is BLACK on a hardcopy black and white device.
 * invras=y is the proper default to make dithered images not come out as negatives
 */
int             invras = YES;
int             window = YES;
int             shade = YES;
int             brake = BREAK_BREAK;
int             framewindows = NO;
int             endpause = NO;
int             cachepipe = NO;
/* Setting cachepipe = YES will copy any piped files to a temporary file,
 * this may get done in dev.open.
 * This is useful if the program may want to reverse seek. e.g. Xvpen, Sunpen
 */
int             allowecho = NEVER_SET;
/*
 * setting allowecho NO (which may get done
 * in dev.open) means that the output device
 * is the user's terminal and echoing chars
 * back at the user would insert nonsense
 * into the plot stream.  In this case we
 * explicitly shut down echoing and output
 * translation until the end of the job.
 * Arguably, output translation should be
 * set/unset in dev.open or dev.reset as we may
 * be plotting on a terminal not our own or
 * the device filter may have built in
 * workarounds for known output translation
 * effects. I use it here as a safety measure.
 */
int             wantras = YES;
int             colormask[5] = {1, 1, 1, 1, 1};

/*
 * filter options - enumerated
 */
int             style = NO_STYLE_YET;
int             default_style = STYLE;
int             rotate;
int             size = RELATIVE;
int             erase = FORCE_INITIAL | DO_LITERALS;

/*
 * filter options - valued
 */
int             xcenter, ycenter;	/* Vplot of point to force as center */
int             fatbase = 0;
int             epause = 0;	/* time to pause before erasing screen */
int             overlay = 0;	/* 1=overlay 0=replace */
int             default_overlay;

/*
 * 0 = none
 * 1 = random
 * 2 = Ordered
 * 3 = Floyd-Steinberg
 * 4 = Halftone
 */
int             dither = 1;	/* Dithering type */

int             xWmax, xWmin, yWmax, yWmin;
float           hshift, vshift;	/* Allow global translation of plot */
float           scale;	/* global scale */
float           xscale;	/* global x-scale */
float           yscale;	/* global y-scale */
float           hdevscale;	/* Vplot units to device units for x */
float           vdevscale;	/* Vplot units to device units for y */
float           txscale;/* global text scale */
float           mkscale;/* global marker scale */
float           dashscale;	/* global dashed line scale */
char            interact[MAXFLEN + 1] = "";	/* Where to store coordinate
						 * file */
float           greyc = 1.;	/* Nonlinear correction */
float           pixc = 1.;	/* Pixel overlap correction */
float           redpow = 1., redmap[4] = {1., 0., 0., 0.};
float           greenpow = 1., greenmap[4] = {0., 1., 0., 0.};
float           bluepow = 1., bluemap[4] = {0., 0., 1., 0.};

/* filter options - resettable between plots */
int             user_rotate;
float           user_txscale;
float           user_mkscale;
float           user_dashscale;
float           user_scale;
float           user_xscale;
float           user_yscale;
int             user_size;
float           user_hshift;
float           user_vshift;
int             user_xwmax_flag;
float           user_xwmax;
int             user_ywmax_flag;
float           user_ywmax;
int             user_xwmin_flag;
float           user_xwmin;
int             user_ywmin_flag;
float           user_ywmin;

int             ifat = 0;
float           fatmult_orig;

/*
 * file and terminal control variables
 */
int             pltoutfd, stderrfd, controlfd;
extern int      genmessage ();
int             (*message) () = genmessage;
struct stat     stderrstat;
struct stat     pltoutstat;
FILE           *pltout, *pltin;
FILE           *fopen ();
FILE           *fdopen ();
FILE           *controltty;
char            outbuf[BUFSIZ];
char           *getenv ();
#if defined(HAVE_STDLIB_H)
#include <stdlib.h>
#else
char           *malloc ();
char           *realloc ();
#endif
char            group_name[MAXFLEN + 1];
int             group_number = 0;
FILE           *pltinarray[MAXIN];
char            pltinname[MAXIN][MAXFLEN + 1];
char            pltname[MAXFLEN + 1] = "";
int             infileno = 0;
extern int      gen_do_dovplot ();
int             (*genreader) () = gen_do_dovplot;

/*
 * Initialize and declare global variables.
 * Use getpar and getenv to search for command-line options.
 */

init_vplot ()
{
char           *stringptr;
int             ii;
char            string[MAXFLEN + 1];
float           ftemp;


    txbuffer = (char *) malloc (TXBUFLEN);
    txbuflen = TXBUFLEN;
/*
 * Sun III lint complains about "pointer alignment problem" here,
 * although the documentation makes it sound like that shouldn't
 * be possible.
 */
    vxbuffer = (struct vertex *) malloc (sizeof (struct vertex) * VXBUFLEN);
    vxbuflen = VXBUFLEN;

    /*
     * If the device can't do color, it can set mono to YES If device is
     * mono, this overrides the user's insistence that it isn't. So GETPAR
     * before we call dev.open. 
     */
    GETPAR ("mono", "1", &mono);

    strcpy (wstype, "default");
    if ((stringptr = getenv ("WSTYPE")) != NULL)
	strcpy (wstype, stringptr);
    /*
     * Don't want to accidentally get this from a history file. Must come for
     * the command line, so getpar, not GETPAR 
     */
    getpar ("wstype", "s", wstype);

    /*
     * Initialize all patterns to be undefined. (device can create a
     * device-dependent default set if it wishes) 
     */
    for (ii = 0; ii <= NPAT; ii++)
    {
	pat[ii].patbits = NULL;
    }

/*
 * The device-independent code MAY NOT actually read from
 * the terminal itself. However, these variables are 
 * declared so that the device-dependent routines can use
 * them if they wish as a courtesy.
 */
    controlfd = open ("/dev/tty", 0);
    if (controlfd == OPEN_ERROR)
    {
	controltty = NULL;
    }
    else
    {
	controltty = fdopen (controlfd, "r");
    }

    /*
     * Call device open before doing anything else. this finalizes pltout 
     */
    dev.open ();
    device_open = YES;

/*
 * Beware trying to print out error messages until
 * we've figured out where to connect the "message" routine!
 */

    if (buffer_output)
	setbuf (pltout, outbuf);
    else
	setbuf (pltout, (char *) NULL);

    /*
     * If graphics output going to control terminal, disable echoing and
     * arrange for use of device's message routine 
     */
    pltoutfd = fileno (pltout);
    stderrfd = fileno (stderr);
    out_isatty = isatty (pltoutfd);

    if (allowecho == NEVER_SET)
    {
	allowecho = YES;
	if (out_isatty)
	{
	    fstat (pltoutfd, &pltoutstat);
	    fstat (stderrfd, &stderrstat);
#ifdef SEP
	    if ((pltoutstat.st_dev == stderrstat.st_dev))
		/*
		 * Something in seplib makes the next 2 tests not work, I
		 * don't know why. Unfortunately, this means the SEP versions
		 * aren't as good about detecting the case of sending the
		 * output to a different terminal than the one you're logged
		 * in at. We should probably get this fixed sometime. - Joe
		 * D. 
		 */
#else
	    if ((pltoutstat.st_dev == stderrstat.st_dev) &&
		(pltoutstat.st_ino == stderrstat.st_ino) &&
		(pltoutstat.st_rdev == stderrstat.st_rdev))
#endif /* SEP */
	    {
		allowecho = NO;
		message = dev.message;
	    }
	}
    }

    /*
     * process YES or NO option arguments. GETPAR means it could have come
     * from a seplib history file, getpar means it had to come from the
     * command line. Use "getpar" when having the right value of the
     * parameter is critical because setting the value overrides frontend's
     * judgement. 
     */
    getpar ("echo", "1", &allowecho);
    if (!allowecho)
    {
#if defined(HAVE_TERMIO_H)
	if (ioctl (pltoutfd, TCGETA, &tty_clean_state) == -1)
	{
		ERR (FATAL, name, "Bad ioctl call!");
	}
	tty_plot_state = tty_clean_state;
	tty_plot_state.c_lflag &= ~ECHO;
	if (ioctl (pltoutfd, TCSETAW, &tty_plot_state) == -1)
	{
		ERR (FATAL, name, "Bad ioctl call! (2)");
	}
#else /* USG */
	ioctl (pltoutfd, TIOCGETP, (char *) (&tty_clean_state));
	bcopy ((char *) (&tty_clean_state), (char *) (&tty_plot_state),
	       sizeof (struct sgttyb));
	ioctl (pltoutfd, TIOCLGET, (char *) (&tty_clean_local_mode));
	tty_plot_local_mode = tty_clean_local_mode | LLITOUT;
	ioctl (pltoutfd, TIOCLSET, (char *) (&tty_plot_local_mode));
	tty_plot_state.sg_flags &= (~ECHO);
	tty_plot_state.sg_flags &= (~LCASE);
	tty_plot_state.sg_flags |= (CBREAK);
	ioctl (pltoutfd, TIOCSETN, (char *) (&tty_plot_state));
#endif /* USG */
    }
    getpar ("endpause", "1", &endpause);
    GETPAR ("cachepipe", "1", &cachepipe);
    GETPAR ("shade", "1", &shade);
    GETPAR ("wantras", "1", &wantras);
    getch ("window", "1", &window);
    getch ("frame", "1", &framewindows);
    GETPAR ("overlay", "1", &overlay);
    default_overlay = overlay;
    GETPAR ("invras", "1", &invras);

    GETPAR ("txsquare", "1", &no_stretch_text);

    GETPAR ("colormask", "1", colormask);
    if (colormask[4] == NO)
    {
	if (colormask[0] + colormask[1] + colormask[2] + colormask[3] != 1)
	    ERR (FATAL, name, "Invalid colormask option.");
    }
    if (colormask[4] == NO && colormask[3] == YES)
    {
	colormask[4] = YES;
    }

    GETPAR ("redpow", "f", &redpow);
    GETPAR ("greenpow", "f", &greenpow);
    GETPAR ("bluepow", "f", &bluepow);

    if (redpow <= 0. || greenpow <= 0. || bluepow <= 0.)
	ERR (FATAL, name, "Invalid color pow option.");

    GETPAR ("red", "f", redmap);
    GETPAR ("green", "f", greenmap);
    GETPAR ("blue", "f", bluemap);

/*
 * Valued arguments
 */

    getpar ("dither", "d", &dither);
    getpar ("greyc", "f", &greyc);
    getpar ("pixc", "f", &pixc);

    GETPAR ("txfont", "d", &txfont);
    GETPAR ("txprec", "d", &txprec);
    GETPAR ("txovly", "d", &txovly);
    default_txfont = txfont;
    default_txprec = txprec;
    default_txovly = txovly;

    if (getch ("erase", "s", string))
    {
	if ((string[0] == 'n') || (string[0] == 'N'))
	    erase = NO;
	else
	if ((string[0] == 'o') || (string[0] == 'O'))
	    erase = FORCE_INITIAL;
	else
	if ((string[0] == 'l') || (string[0] == 'L'))
	    erase = DO_LITERALS;
	else
	    erase = FORCE_INITIAL | DO_LITERALS;
    }

    if (getch ("break", "s", string))
    {
	if (string[0] == 'b')
	    brake = BREAK_BREAK;
	else
	if (string[0] == 'e')
	    brake = BREAK_ERASE;
	else
	    brake = BREAK_IGNORE;
    }

    xcenter = 0;
    ycenter = 0;
    xcenterflag = NO;
    ycenterflag = NO;
    if (GETPAR ("xcenter", "f", &ftemp))
    {
	xcenterflag = YES;
	xcenter = ROUND (ftemp * RPERIN);
    }
    if (GETPAR ("ycenter", "f", &ftemp))
    {
	ycenterflag = YES;
	ycenter = ROUND (ftemp * RPERIN);
    }

    if ((stringptr = getenv ("PATTERNMULT")) != NULL)
    {
	sscanf (stringptr, "%f", &patternmult);
    }
    GETPAR ("patternmult", "f", &patternmult);

    GETPAR ("interact", "s", interact);
    GETPAR ("pause", "d", &epause);

    if (interact[0] != '\0')
    {
	epause = 0;		/* interact makes it own sort of pausing */
	endpause = NO;
    }

    /*
     * Find the default style 
     */
    stringptr = NULL;
    if (GETPAR ("style", "s", string))
	stringptr = string;
    else
	stringptr = getenv ("PLOTSTYLE");
    if (stringptr != NULL)
    {
	if ((stringptr[0] == 'r') || (stringptr[0] == 'R') ||
	    (stringptr[0] == 'm') || (stringptr[0] == 'M'))
	    default_style = ROTATED;
	else
	if ((stringptr[0] == 'o') || (stringptr[0] == 'O'))
	    default_style = OLD;
	else
	if ((stringptr[0] == 'a') || (stringptr[0] == 'A'))
	    default_style = ABSOLUTE;
	else
	    default_style = STANDARD;
    }

/*
 * Options changeable between calls to dovplot
 * (By calling reset_parameters or reset())
 * Dovplot calls reset every time it starts. It only calls
 * reset_parameters the first time.
 *
 * Things in this category:
 * user_*
 * fatmult_orig, fatbase
 */

    if ((stringptr = getenv ("FATMULT")) != NULL)
    {
	sscanf (stringptr, "%f", &fatmult);
    }
    GETPAR ("fatmult", "f", &fatmult);
    fatmult_orig = fatmult;

    user_rotate = 0;
    getch ("rotate", "d", &user_rotate);

    user_txscale = 1.0;
    user_mkscale = 1.0;
    user_dashscale = 1.0;
    GETPAR ("txscale", "f", &user_txscale);
    GETPAR ("mkscale", "f", &user_mkscale);
    GETPAR ("dashscale", "f", &user_dashscale);

    user_scale = 1.0;
    user_xscale = 1.0;
    user_yscale = 1.0;
    getch ("scale", "f", &user_scale);
    getch ("xscale", "f", &user_xscale);
    getch ("yscale", "f", &user_yscale);

    user_size = size;
    if (getch ("size", "s", string))
    {
	if ((string[0] == 'a') || (string[0] == 'A'))
	    user_size = ABSOLUTE;
	else
	    user_size = RELATIVE;
    }

    user_hshift = 0.;
    getch ("hshift xshift", "f", &user_hshift);
    user_vshift = 0.;
    getch ("vshift yshift", "f", &user_vshift);

    user_xwmax_flag = GETPAR ("xwmax", "f", &user_xwmax);
    user_ywmax_flag = GETPAR ("ywmax", "f", &user_ywmax);
    user_xwmin_flag = GETPAR ("xwmin", "f", &user_xwmin);
    user_ywmin_flag = GETPAR ("ywmin", "f", &user_ywmin);

    GETPAR ("fat", "d", &fatbase);



/*
 *  These parameters can simply be changed at any time with no
 *  need for re-initialization of anything else:
 *
 *  shade, wantras
 *
 */

}

init_colors ()
{
int             ii;

    /*
     * Set up the default color table 
     */

/*
 * First 7 colors are assumed to map to the standard 7 pen colors
 * on ALL devices, even those with no settable colors!
 */
    for (ii = 0; ii < 8; ii++)
    {
/*     if(color_set[ii][STATUS]!=SET){*/
	    color_set[ii][STATUS] = SET;
	    color_set[ii][_RED] = MAX_GUN * ((ii & 2) / 2);
	    color_set[ii][_GREEN] = MAX_GUN * ((ii & 4) / 4);
	    color_set[ii][_BLUE] = MAX_GUN * ((ii & 1) / 1);
	    color_set[ii][_GREY] = greycorr ((int) ((MAX_GUN * ii) / 7));
/*   }*/
    }

    if (mono)
    {
	/* Monochrome devices are assumed to have no settable colors */
	num_col = 0;
    }

    num_col_8 = (num_col > 8) ? num_col : 8;

    if (mono)
    {
	color_set[0][MAP] = 0;
	for (ii = 1; ii <= MAX_COL; ii++)
	{
/*     if(color_set[ii][STATUS]!=SET)*/
	    color_set[ii][MAP] = 7;
	}
	for (ii = num_col_8; ii < 256; ii++)
	{
/*     if(color_set[ii][STATUS]!=SET)*/
	    color_set[ii][_GREY] = color_set[((ii - 8) % 7) + 1][_GREY];
	}
	/* Put a grey scale in the upper half of the color table */
	for (ii = 256; ii <= MAX_COL; ii++)
	{
/*     if(color_set[ii][STATUS]!=SET)*/
	    color_set[ii][_GREY] = greycorr (ii - 256);
	}
    }
    else
    {
	/*
	 * Unmapped colors shouldn't be mapped; ie, they map to themselves 
	 */

	for (ii = 0; ii < num_col_8; ii++)
	{
/*     if(color_set[ii][STATUS]!=SET)*/
	    color_set[ii][MAP] = ii;
/*  fprintf(stderr,"MAPPING %d %d - %d \n",ii,MAP,color_set[ii][MAP]);*/
	 }
	/*
	 * Colors outside the range of this terminal map cyclically back into
	 * colors 1 through 7 
	 */
	for (ii = num_col_8; ii <= MAX_COL; ii++)
	{
/*     if(color_set[ii][STATUS]!=SET)*/
	    color_set[ii][MAP] = ((ii - 8) % 7) + 1;
	}
    }
}

setstyle (new_style)
    int             new_style;
{
    /*
     * Check to see if the style has changed 
     */
    if (new_style == style)
	return;

    style = new_style;
    reset_parameters ();
}

reset_parameters ()
{
float           inches;	/* scaling base for y axis */
float           screenheight, screenwidth;	/* true size of the screen */
int             ix, iy;

    xorigin = 0;
    yorigin = 0;

    rotate = default_rotate;
    rotate += user_rotate;

    txscale = user_txscale;
    mkscale = user_mkscale;
    dashscale = user_dashscale;

    scale = default_scale * user_scale;
    xscale = default_xscale * user_xscale;
    yscale = default_yscale * user_yscale;

    fatmult = fatmult_orig;

    size = user_size;

    switch (style)
    {
/*
 * The old standard on the machine erebus.
 * Pretty much dead now. Still useful for some old programs nobody's
 * wanted to update.
 */
    case OLD:
	txscale /= 3.;
	fatmult /= 3.;
	scale *= 3.;
	inches = STANDARD_HEIGHT;
	break;
/*
 * The old standard on the machine mazama. A useful coordinate system
 * for geological sorts of plots. 
 * The Y axis goes the long way, across the screen, which is the device's
 * horizontal axis.
 */
    case ROTATED:
	rotate += 90;
	inches = ROTATED_HEIGHT;
	break;
    case ABSOLUTE:
	size = ABSOLUTE;
    case STANDARD:
    default:
	inches = STANDARD_HEIGHT;
	break;
    }

    if (rotate >= 0)
	rotate = rotate % 360;
    else
	rotate = ((rotate % 360) + 360) % 360;

    mxx = cos (2. * 3.14159 * rotate / 360.);
    myy = cos (2. * 3.14159 * rotate / 360.);
    mxy = sin (2. * 3.14159 * rotate / 360.);
    myx = -sin (2. * 3.14159 * rotate / 360.);

    if (size == ABSOLUTE)
    {
	vdevscale = pixels_per_inch / (float) (RPERIN * aspect_ratio);
	hdevscale = pixels_per_inch / (float) RPERIN;
    }
    else
    {
	/*
	 * Fit the inches x inches unit square into a displayable box with
	 * aspect ratio SCREEN_RATIO 
	 */
	screenwidth = (dev_xmax - dev_xmin) * pixels_per_inch;
	screenheight =
	 (dev_ymax - dev_ymin) * pixels_per_inch * aspect_ratio;
	if ((screenheight / screenwidth) > SCREEN_RATIO)
	{
	    vdevscale = (SCREEN_RATIO * ((dev_xmax - dev_xmin) / aspect_ratio)) /
	     (inches * RPERIN);
	    hdevscale = vdevscale * aspect_ratio;
	}
	else
	{
	    vdevscale = (dev_ymax - dev_ymin) / (inches * RPERIN);
	    hdevscale = vdevscale * aspect_ratio;
	}
    }

    hshift = default_hshift;
    vshift = default_vshift;

    if (style == ROTATED)
    {
	vshift += dev_ymax - dev_ymin;
    }

    yscale *= scale;
    xscale *= scale;
    mkscale *= scale;

/*
 * Set up fatness multiplication factor
 */
    fatmult *= scale * hdevscale * RPERIN / FATPERIN;

    /*
     * The point (xcenter,ycenter) in vplot coordinates is to be centered in
     * the screen. 
     */
    if (xcenterflag || ycenterflag)
    {
	vptodevxy (xcenter, ycenter, &ix, &iy);
	if (xcenterflag)
	{
	    hshift += (dev_xmax + dev_xmin) / 2 - ix;
	}
	if (ycenterflag)
	{
	    vshift += (dev_ymax + dev_ymin) / 2 - iy;
	}
    }

    hshift += user_hshift * pixels_per_inch;
    vshift += user_vshift * pixels_per_inch / aspect_ratio;

/* plot window parameters defaulted */
/* to maximum size */

    devtovpw (dev_xmin, dev_ymin, dev_xmax, dev_ymax,
	      &xWmin, &yWmin, &xWmax, &yWmax);

    if (user_xwmax_flag)
	xWmax = ROUND (user_xwmax * RPERIN);
    if (user_ywmax_flag)
	yWmax = ROUND (user_ywmax * RPERIN);
    if (user_xwmin_flag)
	xWmin = ROUND (user_xwmin * RPERIN);
    if (user_ywmin_flag)
	yWmin = ROUND (user_ywmin * RPERIN);

    vptodevw (xWmin, yWmin, xWmax, yWmax, &xWmin, &yWmin, &xWmax, &yWmax);

    wlimit (dev_xmin, dev_xmax, &xWmin, &xWmax);
    wlimit (dev_ymin, dev_ymax, &yWmin, &yWmax);

    xwmax = xWmax;		/* plot window parameters defaulted */
    xwmin = xWmin;		/* to maximum size 		 */
    ywmax = yWmax;
    ywmin = yWmin;
    reset_windows ();
}

add_a_cor (filename, xcor, ycor)
    char           *filename;
    int             xcor, ycor;
{
static int      first_time = YES;
static FILE    *outfp;

    if (first_time == YES)
    {
	outfp = fopen (filename, "w");
	if (outfp == NULL)
	{
	    ERR (FATAL, name, "Can't open interact output file %s!", filename);
	}
	first_time = NO;
    }
    fprintf (outfp, "%f\t%f\n", (float) xcor / RPERIN, (float) ycor / RPERIN);
}

wlimit (min, max, xmin, xmax)
    int             min, max;
    int            *xmin, *xmax;
{
    if (*xmin < min)
	*xmin = min;

    if (*xmax > max)
	*xmax = max;
}
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
 *  source file:   ./filters/main_vplot.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Stew Levin (SEP), December 7, 1993
 *      Dummy routine MAIN_ for shared library linkage
 * Stew Levin (Mobil) May 8, 1996
 *      Updated signal handling #ifdef's for LINUX
 *
 */

/*
 * generic pen -  VPLOT filter for whatever
 * Keyword: graphics vplot pen
 */

/*
 *  Edit History
 *
 *  "Pen" written 1979 for a PDP-11, Rob Clayton
 *  Eventually modified into "tekpen" by Glenn Kroeger early 1984
 *  Made device independent and extensively modified into "screenpen"
 *  by Joe Dellinger 1984 through 1985
 *  Reworked to be more GKS-like by Glenn and Michel Debiche, 1985
 *  Cleaned up, Joe 1986
 *  Added 'wstype' variable so one program can support multiple
 *		terminal types		-- Chuck Karish, Nov 1986
 *  Raster capability added -- Joe and Steve Cole 1987
 *  Cleaned up and documented for public release, Joe Dellinger 1987
 *  Whew!
 *
 *  Changed system("stty ...") calls to ioctls.  The system() calls
 *  are hanging for jon under 6.0.  Added sigblock() to cleanup
 *  error handler to avoid exit() problems.  Also shut down output
 *  translations to tty graphics terminals.  Stewart A. Levin  6-23-87
 *
 *  Shut down output translations with LLITOUT now that I've gotten
 *  the magic incantation that makes it stick on buggy BSD systems.
 *  Stewart A. Levin  7-5-87
 *
 *  Made "scale" a global so dovplot could use it.
 *  Joe Dellinger Oct 18 1987
 *
 *  Gather up all the incoming plot files, and then have "genreader"
 *  process them all. This allows someone to make a vplot-editor that
 *  can re-order plot files, repeat them, etc, etc. Added "buffer_input"
 *  and "allow_pipe" as part of this.
 *  Joe Dellinger Dec 16 1987
 *
 *  Added "hclose_done" to SEPlib version.
 *  Joe Dellinger Dec 19 1987
 *
 *  Added "reset_parameters" to frontend and dovplot
 *  Joe Dellinger Jan 8 1988
 *
 *  Added group_name, group_number, pltname
 *  Joe Dellinger Jan 20 1988
 *
 *  Inverted window should mean that everything gets clipped.
 *  Make it possible to override standard defaults for SEP.
 *  Joe Dellinger Feb 12 1988
 *
 *  Just use getpar, not getpar_.
 *  Joe Dellinger Feb 16 1988
 *
 *  Split "frontend.c" into 3 files,
 *  main_vplot, init_vplot, and proc_vplot.
 *  This so that other programs can use vplot and still keep their own
 *  mains.
 *  Joe Dellinger Feb 18 1988
 *
 *  Made some SEPlib warning messages a bit clearer.
 *  Joe Dellinger June 30 1988
 *
 * Dave Nichols (SEP), March 25 1989
 * Added cachepipe option to write piped input to a temporary file.
 * This allows them to be reread and replotted under X-window or Sunview.
 *
 * Dave Nichols (SEP), March 27 1989
 * Make sure the temporary file is cleaned up and that seplib input really
 * is coming from a pipe before it is copied.  
 * W. Bauske IBM, 03-27-91
 *	Applied SysV mods for RS/6000 
 *	Removed external re-declare of malloc for RS/6000 
 *
 * Dave Nichols (SEP), March 30 1992
 *  	Use sepxargv and sepxargc for traversing argument list to be consistent.
 *	This can be a problem is Xt initialisation deletes arguments.
 *
 * Dave Nichols (SEP), May 2 1992
 *	Check for environment variable VPLOTSPOOLDIR when creating temporary
 *	files.
 * Dave Nichols (SEP), Dec 11 1992
 *	Added check for zero length files on stdin (e.g. /dev/null )
 * Stew Levin 2-27-95
 *	Solaris mods
 *  Bob Clapp 10-98  Switched to POSIX (ala Sloaris) for LINUX signals
 */

#include<sitedef.h>
#include <stdlib.h>

#ifdef SEP
extern void     sepwhere ();
extern char     sepoutwhere[];
extern char     sepheadwhere[];

#define		OUT	sepoutwhere
#define		HEAD	sepheadwhere
#if !defined(RS6000) && !defined(SOURCE)
#define		SOURCE  "\014Joe Dellinger, Stanford Exploration Project\014"
#endif
#include	<string.h>
#include	<sep.main>
#define		GETPAR	fetch

#else /* SEP */
#include	<stdio.h>
#include	<math.h>
#include	<string.h>
#define		GETPAR	getpar
#endif /* SEP */


#if defined(HAVE_TERMIO_H)
#include	<termio.h>
#else
#if defined (HAVE_SGTTY_H)
#include	<sgtty.h>
#else
#include	<sys/ioctl.h>
#include	<sgtty.h>
#endif
#endif /* USG */
#include	<sys/types.h>
#include	<sys/stat.h>
#include	<ctype.h>
#include	<signal.h>

#include	<vplot.h>

#include	"./include/params.h"	/* for machine dependencies */
#include	"./include/enum.h"
#include	"./include/err.h"
#include	"./include/attrcom.h"
#include	"./include/intcom.h"
#include	"./include/mesgcom.h"
#include	"./include/erasecom.h"
#include	"./include/closestat.h"
#include	"./include/pat.h"
#include	"./include/vertex.h"
#include	"./include/round.h"
#include	"./include/extern.h"


#if defined(HAVE_TERMIO_H)
#else /* USG */
/*
 * signal catching
 */
#ifdef SIGFNC_RTN_VOID
void            cleanup ();
#else
int             cleanup ();
#endif
int             signum[] =
{
#ifdef LINUX
 SIGHUP, SIGINT, SIGQUIT, SIGIOT, SIGBUS, SIGPIPE, SIGTERM, SIGXCPU, SIGXFSZ
#else
 SIGHUP, SIGINT, SIGQUIT, SIGIOT, SIGEMT, SIGPIPE, SIGTERM, SIGXCPU, SIGXFSZ
#endif
};
#define NOSIG (sizeof (signum)/sizeof (int))	/* number of signals caught */
#if defined(SOLARIS ) || defined(LINUX)
struct sigaction   errhandler =
{
#ifdef LINUX
 cleanup, 0, 0
#else
 0, cleanup, 0
#endif
};
struct sigaction   ignored =
{
#ifdef LINUX
 SIG_IGN, 0, 0
#else
 0, SIG_IGN, 0
#endif
};
struct sigaction   oldvec;
#else /*SOLARIS*/
int             sigvec ();
struct sigvec   errhandler =
{
 cleanup, 0, 0
};
struct sigvec   ignored =
{
 SIG_IGN, 0, 0
};
struct sigvec   oldvec;
#endif /*SOLARIS*/
#endif /* USG */

#if defined(HAVE_TERMIO_H)
extern struct termio tty_clean_state;
#else /* USG */
extern struct sgttyb tty_clean_state;	/* external for utilities */
extern int      tty_clean_local_mode;
#endif /* USG */
extern int      allow_pipe;
extern char     callname[];
extern int      nplots;
extern int      allowecho;
extern int      cachepipe;

/*
 * file and terminal control variables
 */
extern int      genmessage ();
extern int      (*message) ();
extern FILE    *pltout;

FILE           *fopen ();
FILE           *fdopen ();
FILE	       *tempcopy();
extern int	unlink();

extern FILE    *pltinarray[MAXIN];
extern char     pltinname[MAXIN][MAXFLEN + 1];
extern int      infileno;
extern int      pltoutfd;

/*
 * This routine is responsible for finding the input files,
 * setting up the input and output, and calling init_vplot
 * and proc_vplot. You can link to vplot without using this
 * routine, and so have your own main. See vplothacker.doc to
 * learn how to do this. (Especially you, Jon!)
 */

#ifdef SEP
int             xsepxargc;
char          **xsepxargv;
char            scrap[MAXFLEN + 1];
int             hclose_done = NO;
int             fake_header = NO;
MAIN ()

#else /* SEP */
int             sepxargc;
char          **sepxargv;	/* for getpar */
MAIN_(){} /* dummy for shared linkage */
main (argc, argv)
    int             argc;
    char           *argv[];
#endif /* SEP */
{
#ifndef SEP
int             in_isatty, num_vplot, docflag;
char            instring[MAXFLEN + 1];
#endif /* SEP */

char           *cptr;
char           *stringptr;
int             ii;
FILE           *temp;
char            string[MAXFLEN + 1];
int		tempfileindex = -1;

    nulldev ();			/* Just to make sure it gets loaded */

#ifndef SEP
    if (stringptr = strrchr(argv[0], '/'))
	strncpy (callname, ++stringptr, 24);
    else
	strncpy (callname, argv[0], 24);
#else /* SEP */
    if (stringptr = strrchr(sepxargv[0], '/'))
	strncpy (callname, ++stringptr, 24);
    else
	strncpy (callname, sepxargv[0], 24);
#endif /* SEP */

#ifdef SEP
    pltout = outstream;
    if (redout ())
    {
	getch ("head", "s", scrap);
	if (strcmp (scrap, "/dev/null") == 0 &&
	    strcmp (sepheadwhere, "/dev/null") == 0)
	{
	    fake_header = YES;
	    getch ("out", "s", scrap);
	    if (strcmp (scrap, "stdout") != 0)
	    {
		/*
		 * They probably want the header output into the redirected
		 * output. (SEP only) 
		 */
		headstream = stdout;
		headfd = fileno (headstream);
		Puthead ("Fake header for special device-dependent data only.\n");
		Puthead ("To get input history passed along, over-ride default with head = stdout\n");
		/*
		 * If the output is going into a file, then put this
		 * information into the header file. 
		 */
		if (strcmp (scrap, "/dev/tty") != 0)
		{
		    fullnm (scrap, MAXFLEN + 1);
		    Puthead ("\tin=%s\n", scrap);
		}
		else
		{
		    Puthead ("\tin= nowhere\n");
		    Puthead ("\t(Sorry, the data went to the terminal; it's gone now!)\n");
		}

		/*
		 * Have to use my own puthead routine. Standard ones such as
		 * putch, etc, don't work with this. They remember where the
		 * header used to be. Puthead does this, checking for
		 * headstream not null. This is so this will work without
		 * sep, as well. 
		 */
	    }
	}
    }

    if (!fake_header)
    {
	Puthead ("\tn3 = unknown\n\t(Sorry, SEPlib requires the header be closed now.)\n");
	Puthead ("\t(Normally you won't care that n3 is unknown unless you're using Raspen!\n)");
	hclose ();
	hclose_done = YES;
    }
#else /* SEP */

    /*
     * If no arguments, and not in a pipeline, self document "wstype="
     * doesn't count as an argument for our purposes 
     */
    in_isatty = isatty ((int) (fileno (stdin)));
    sepxargc = argc;
    sepxargv = argv;
    docflag = 0;
    if (argc == 1)
	docflag = 1;
    if ((argc == 2) && !strncmp ("wstype=", argv[1], 7))
	docflag = 1;
    getpar ("selfdoc", "1", &docflag);
    if (in_isatty && docflag)
    {
	for (ii = 0; ii < doclength; ii++)
	    printf ("%s\n", documentation[ii]);
	exit (0);
    }

    pltout = stdout;
#endif /* SEP */

#if defined(HAVE_TERMIO_H)

#else /* USG */
    /*
     * This getpar for signal is only included for debugging purposes. By
     * using a signal option, one can stop any signals from being caught. 
     */
    if (getpar ("signal", "s", string) == 0)
    {
/*#ifdef SOLARIS*/
#if defined(SOLARIS) || defined(LINUX)
        sigfillset(&(errhandler.sa_mask));
#endif
	for (ii = 0; ii < NOSIG; ++ii)
	{
#if defined(SOLARIS) || defined(LINUX)
	    if (-1 == sigaction (signum[ii], &ignored, &oldvec))
	    {
		ERR (FATAL, name, "Bad sigvec call!");
	    }
	    if (oldvec.sa_handler == ignored.sa_handler)
		(void) sigaction (signum[ii], &oldvec, (struct sigaction *) NULL);
	    else
		(void) sigaction (signum[ii], &errhandler, (struct sigaction *) NULL);
#else
	    if (-1 == sigvec (signum[ii], &ignored, &oldvec))
	    {
		ERR (FATAL, name, "Bad sigvec call!");
	    }
	    if (oldvec.sv_handler == ignored.sv_handler)
		(void) sigvec (signum[ii], &oldvec, (struct sigvec *) NULL);
	    else
		(void) sigvec (signum[ii], &errhandler, (struct sigvec *) NULL);
#endif
	}
    }
#endif /* USG */

/*
 ****************************************************************************
 * Set all global variables, open the device.
 ****************************************************************************
 */

    init_vplot ();

/*
 ****************************************************************************
 * Start processing input files
 ****************************************************************************
 */


#ifdef SEP
    if (instream != NULL)
    {
	if (infileno >= MAXIN)
	{
	    ERR (FATAL, name, "too many input files (%d max)", MAXIN);
	}
	if( cachepipe && isapipe(fileno( instream)) )
	{
	    if( (pltinarray[infileno] = tempcopy(instream,string) ) == NULL )
	    {
	    	ERR( FATAL, name, "copy of piped input failed");
	    }
	 
            /* check for zero length input (e.g. /dev/null ) */
	    if( pltinarray[infileno] != (FILE*) -1 ) {

	    strcpy( pltinname[infileno], string );
	    /*remember what number this file is so we can delete it later*/
	    tempfileindex = infileno;
	    infileno++;
	    }
	}
	else
	{
	    if (!allow_pipe && isapipe (fileno (instream)))
	    {
	    	ERR (WARN, name, "cannot use pipes with this device, try cachepipe=y ");
	    }
	    else
	    {
	        strcpy (pltinname[infileno], "Pipe");
	    	pltinarray[infileno] = instream;
	    	infileno++;
	    }
	}
    }
    else
	ERR (WARN, name, "cannot read input pipe");

    xsepxargc = sepxargc;
    xsepxargv = sepxargv;

    for (xsepxargc--, xsepxargv++; xsepxargc; xsepxargc--, xsepxargv++)
    {
	cptr = *xsepxargv;
	while (*cptr)
	{
	    if (*cptr == '=')
		break;
	    cptr++;
	}
	if (*cptr)
	    continue;
	/* Ignore dummy arguments */
	if (strcmp (*xsepxargv, "dummy") == 0)
	    continue;
	if ((temp = fopen (*xsepxargv, "r")) == NULL)
	{
	    ERR (WARN, name, "cannot open header file %s", *xsepxargv);
	    continue;
	}
	fclose (temp);
	if (getch2 ("in", "s", string, *xsepxargv))
	{
	    if ((temp = fopen (string, "r")) != NULL)
	    {
		Puthead ("   +  %s --> in = %s\n", *xsepxargv, string);
		if (infileno >= MAXIN)
		{
		    ERR (FATAL, name, "too many input files (%d max)", MAXIN);
		}
		strcpy (pltinname[infileno], string);
		pltinarray[infileno] = temp;
		infileno++;
	    }
	    else
	    {
		ERR (WARN, name, "cannot open input file %s", string);
	    }
	}
    }

#else /* SEP */
    /*
     * first process pipe input 
     */
    if (!in_isatty)
    {
	if (infileno >= MAXIN)
	{
	    ERR (FATAL, name, "too many input files (%d max)", MAXIN);
	}

  	if( cachepipe )
        {
            if( (pltinarray[infileno] = tempcopy( stdin,string) ) == NULL )
            {
                ERR( FATAL, name, "copy of piped input failed");
	    }

            /* check for zero length input (e.g. /dev/null ) */
	    if( pltinarray[infileno] != (FILE*)-1 ) {

	    strcpy( pltinname[infileno], string );
	    /* remember what number this file is so we can delete it later*/
	    tempfileindex = infileno;
            infileno++;
	    }
        }
        else
        {
	    if (!allow_pipe)
	    {
	    	ERR (WARN, name, "cannot use pipes with this device, try cachepipe=y ");
	    }
	    else
	    {
	    	strcpy (pltinname[infileno], "stdin");
	    	pltinarray[infileno] = stdin;
	    	infileno++;
	    }
	}
    }

    /*
     * next process in= inputfiles If they set num_vplot, also look for in1=
     * in2= etc 
     */

    num_vplot = 0;
    getpar ("numvplot", "d", &num_vplot);

    for (ii = 0; ii <= num_vplot; ii++)
    {
	if (ii == 0)
	    strcpy (instring, "in");
	else
	    sprintf (instring, "in%d", ii);

	if (getpar (instring, "s", string))
	{
	    if ((temp = fopen (string, "r")) != NULL)
	    {
		if (infileno >= MAXIN)
		{
		    ERR (FATAL, name, "too many input files (%d max)", MAXIN);
		}
		strcpy (pltinname[infileno], string);
		pltinarray[infileno] = temp;
		infileno++;
	    }
	    else
	    {
		ERR (WARN, name, "cannot open %s", string);
	    }
	}
    }

    /*
     * finally process input line for non-getpar arguments and assume they
     * are also input files 
     */
    for (sepxargc--, sepxargv++; sepxargc; sepxargc--, sepxargv++)
    {
	cptr = *sepxargv;
	while (*cptr)
	{
	    if (*cptr == '=')
		break;
	    cptr++;
	}
	if (*cptr)
	    continue;
	cptr = *sepxargv;
	if ((temp = fopen (cptr, "r")) != NULL)
	{
	    if (infileno >= MAXIN)
	    {
		ERR (FATAL, name, "too many input files (%d max)", MAXIN);
	    }
	    strcpy (pltinname[infileno], cptr);
	    pltinarray[infileno] = temp;
	    infileno++;
	}
	else
	{
	    ERR (WARN, name, "cannot open %s", cptr);
	}
    }
#endif /* SEP */

/*
 ****************************************************************************
 * Go do the plots
 ****************************************************************************
 */

    proc_vplot ();

#ifdef SEP
    if (!hclose_done)
    {
	Puthead ("\tn3=%d\n", nplots);
	hclose ();
	hclose_done = YES;
    }
#endif /* SEP */

    /*  delete the temporary copy of piped input if there is one*/
    removtemp();

    /* All done */
    exit (0);
    return 0;
}

#if defined(HAVE_TERMIO_H)
#else /* USG */
#ifdef SIGFNC_RTN_VOID
void cleanup ()
#else
cleanup ()
#endif
{
#ifndef SOLARIS
#ifndef LINUX
    sigblock (~(SIGKILL | SIGSTOP | SIGCONT));
#endif
#endif
    dev.close (CLOSE_INTERRUPT);
    message (MESG_ON);
    ERR (COMMENT, name, "Interrupted out.");
    dev.close (CLOSE_DONE);
    /*  delete the temporary copy of piped input if there is one*/
    removtemp();
    /*
     * Let them see what they are doing again 
     */
    if (!allowecho)
    {
/*#ifdef SOLARIS*/
#if defined(SOLARIS) || defined(LINUX)
	ioctl (pltoutfd, TCSETAW, (char *) (&tty_clean_state));
#else
	ioctl (pltoutfd, TIOCLSET, (char *) (&tty_clean_local_mode));
	ioctl (pltoutfd, TIOCSETN, (char *) (&tty_clean_state));
#endif
    }
    exit (0);
}
#endif /* USG */

/* routine to copy a file to a temporary file, used to copy stdin so that
 * it can be reread as required
 */

#define XFER_SIZE 1024

/* the name of the temporary file */
static char *tspoolnm=(char*)NULL;

removtemp(){
    if( tspoolnm != (char*)NULL )
        if( unlink(tspoolnm) )
               ERR (WARN, name, "unable to delete temporary file");
    tspoolnm = (char*)NULL;
}

FILE* tempcopy( infile, filename )
FILE* infile;
char* filename;
{
    FILE *temp;
    int len,total;
    char xfer_buf[ XFER_SIZE];
    char* spooldirnm;


    /* make up a temporary file name 
     * I know there are beter routines to make up the name 
     * but this one is more portable.
     * PEN_SPOOL is a directory we can put temporary files in
     * it is defined in params.h. It can be overridden by the
     * environment variable VPLOTSPOOLDIR.
     */
    if( (spooldirnm = getenv("VPLOTSPOOLDIR")) != NULL ){
	tspoolnm = (char *) malloc( strlen( spooldirnm ) +13 );
        strcpy( tspoolnm, spooldirnm );
    }else{
        tspoolnm = (char *) malloc( strlen( PEN_SPOOL ) + 13 );
        strcpy( tspoolnm, PEN_SPOOL );
    }
    tspoolnm = strcat( tspoolnm, "/vplotXXXXXX" );
    tspoolnm = mktemp( tspoolnm );

    if ( (temp = fopen(tspoolnm, "w")) == NULL) {
	ERR (WARN, name, "unable to create temporary file");
	return NULL;
    }

    /* now copy everything from the input stream to our temporary file */
    total=0;
    while( (len=fread( xfer_buf, sizeof(char), XFER_SIZE, infile ) ) != 0 ){
	if( ferror(infile) ) {
	    ERR (WARN, name, "read error from pipe ");
	    return NULL;
	}
	    
        if( ( fwrite(xfer_buf,sizeof(char),len,temp) != len ) || ferror(temp) ) 	{
	    ERR (WARN, name, "write error to temp ");
	    return NULL;
	}
	total += len;
    }

    fclose(temp);

    /* We could unlink the file after we reopen it to make it really
     * temporary but we don't know for sure that a filter isn't going
     * to close the file and try to reopen it.
     * Hence the buisness with tempfileindex in the main routine.
     */

    /* check for zero length input */
    if( total == 0 ) return (FILE*)-1;

    strcpy( filename, tspoolnm );
    return fopen( tspoolnm, "r");
}

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
 *  source file:   ./filters/main_vplot.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Stew Levin (SEP), December 7, 1993
 *      Dummy routine MAIN_ for shared library linkage
 * Stew Levin (Mobil) May 8, 1996
 *      Updated signal handling #ifdef's for LINUX
 *
 */

/*
 * generic pen -  VPLOT filter for whatever
 * Keyword: graphics vplot pen
 */

/*
 *  Edit History
 *
 *  "Pen" written 1979 for a PDP-11, Rob Clayton
 *  Eventually modified into "tekpen" by Glenn Kroeger early 1984
 *  Made device independent and extensively modified into "screenpen"
 *  by Joe Dellinger 1984 through 1985
 *  Reworked to be more GKS-like by Glenn and Michel Debiche, 1985
 *  Cleaned up, Joe 1986
 *  Added 'wstype' variable so one program can support multiple
 *		terminal types		-- Chuck Karish, Nov 1986
 *  Raster capability added -- Joe and Steve Cole 1987
 *  Cleaned up and documented for public release, Joe Dellinger 1987
 *  Whew!
 *
 *  Changed system("stty ...") calls to ioctls.  The system() calls
 *  are hanging for jon under 6.0.  Added sigblock() to cleanup
 *  error handler to avoid exit() problems.  Also shut down output
 *  translations to tty graphics terminals.  Stewart A. Levin  6-23-87
 *
 *  Shut down output translations with LLITOUT now that I've gotten
 *  the magic incantation that makes it stick on buggy BSD systems.
 *  Stewart A. Levin  7-5-87
 *
 *  Made "scale" a global so dovplot could use it.
 *  Joe Dellinger Oct 18 1987
 *
 *  Gather up all the incoming plot files, and then have "genreader"
 *  process them all. This allows someone to make a vplot-editor that
 *  can re-order plot files, repeat them, etc, etc. Added "buffer_input"
 *  and "allow_pipe" as part of this.
 *  Joe Dellinger Dec 16 1987
 *
 *  Added "hclose_done" to SEPlib version.
 *  Joe Dellinger Dec 19 1987
 *
 *  Added "reset_parameters" to frontend and dovplot
 *  Joe Dellinger Jan 8 1988
 *
 *  Added group_name, group_number, pltname
 *  Joe Dellinger Jan 20 1988
 *
 *  Inverted window should mean that everything gets clipped.
 *  Make it possible to override standard defaults for SEP.
 *  Joe Dellinger Feb 12 1988
 *
 *  Just use getpar, not getpar_.
 *  Joe Dellinger Feb 16 1988
 *
 *  Split "frontend.c" into 3 files,
 *  main_vplot, init_vplot, and proc_vplot.
 *  This so that other programs can use vplot and still keep their own
 *  mains.
 *  Joe Dellinger Feb 18 1988
 *
 *  Made some SEPlib warning messages a bit clearer.
 *  Joe Dellinger June 30 1988
 *
 * Dave Nichols (SEP), March 25 1989
 * Added cachepipe option to write piped input to a temporary file.
 * This allows them to be reread and replotted under X-window or Sunview.
 *
 * Dave Nichols (SEP), March 27 1989
 * Make sure the temporary file is cleaned up and that seplib input really
 * is coming from a pipe before it is copied.  
 * W. Bauske IBM, 03-27-91
 *	Applied SysV mods for RS/6000 
 *	Removed external re-declare of malloc for RS/6000 
 *
 * Dave Nichols (SEP), March 30 1992
 *  	Use sepxargv and sepxargc for traversing argument list to be consistent.
 *	This can be a problem is Xt initialisation deletes arguments.
 *
 * Dave Nichols (SEP), May 2 1992
 *	Check for environment variable VPLOTSPOOLDIR when creating temporary
 *	files.
 * Dave Nichols (SEP), Dec 11 1992
 *	Added check for zero length files on stdin (e.g. /dev/null )
 * Stew Levin 2-27-95
 *	Solaris mods
 *  Bob Clapp 10-98  Switched to POSIX (ala Sloaris) for LINUX signals
 */

#include<sitedef.h>
#include <stdlib.h>

#ifdef SEP
extern void     sepwhere ();
extern char     sepoutwhere[];
extern char     sepheadwhere[];

#define		OUT	sepoutwhere
#define		HEAD	sepheadwhere
#if !defined(RS6000) && !defined(SOURCE)
#define		SOURCE  "\014Joe Dellinger, Stanford Exploration Project\014"
#endif
#include	<string.h>
#include	<sep.main>
#define		GETPAR	fetch

#else /* SEP */
#include	<stdio.h>
#include	<math.h>
#include	<string.h>
#define		GETPAR	getpar
#endif /* SEP */


#if defined(HAVE_TERMIO_H)
#include	<termio.h>
#else
#if defined (HAVE_SGTTY_H)
#include	<sgtty.h>
#else
#include	<sys/ioctl.h>
#include	<sgtty.h>
#endif
#endif /* USG */
#include	<sys/types.h>
#include	<sys/stat.h>
#include	<ctype.h>
#include	<signal.h>

#include	<vplot.h>

#include	"./include/params.h"	/* for machine dependencies */
#include	"./include/enum.h"
#include	"./include/err.h"
#include	"./include/attrcom.h"
#include	"./include/intcom.h"
#include	"./include/mesgcom.h"
#include	"./include/erasecom.h"
#include	"./include/closestat.h"
#include	"./include/pat.h"
#include	"./include/vertex.h"
#include	"./include/round.h"
#include	"./include/extern.h"


#if defined(HAVE_TERMIO_H)
#else /* USG */
/*
 * signal catching
 */
#ifdef SIGFNC_RTN_VOID
void            cleanup ();
#else
int             cleanup ();
#endif
int             signum[] =
{
#ifdef LINUX
 SIGHUP, SIGINT, SIGQUIT, SIGIOT, SIGBUS, SIGPIPE, SIGTERM, SIGXCPU, SIGXFSZ
#else
 SIGHUP, SIGINT, SIGQUIT, SIGIOT, SIGEMT, SIGPIPE, SIGTERM, SIGXCPU, SIGXFSZ
#endif
};
#define NOSIG (sizeof (signum)/sizeof (int))	/* number of signals caught */
#if defined(SOLARIS ) || defined(LINUX)
struct sigaction   errhandler =
{
#ifdef LINUX
 cleanup, 0, 0
#else
 0, cleanup, 0
#endif
};
struct sigaction   ignored =
{
#ifdef LINUX
 SIG_IGN, 0, 0
#else
 0, SIG_IGN, 0
#endif
};
struct sigaction   oldvec;
#else /*SOLARIS*/
int             sigvec ();
struct sigvec   errhandler =
{
 cleanup, 0, 0
};
struct sigvec   ignored =
{
 SIG_IGN, 0, 0
};
struct sigvec   oldvec;
#endif /*SOLARIS*/
#endif /* USG */

#if defined(HAVE_TERMIO_H)
extern struct termio tty_clean_state;
#else /* USG */
extern struct sgttyb tty_clean_state;	/* external for utilities */
extern int      tty_clean_local_mode;
#endif /* USG */
extern int      allow_pipe;
extern char     callname[];
extern int      nplots;
extern int      allowecho;
extern int      cachepipe;

/*
 * file and terminal control variables
 */
extern int      genmessage ();
extern int      (*message) ();
extern FILE    *pltout;

FILE           *fopen ();
FILE           *fdopen ();
FILE	       *tempcopy();
extern int	unlink();

extern FILE    *pltinarray[MAXIN];
extern char     pltinname[MAXIN][MAXFLEN + 1];
extern int      infileno;
extern int      pltoutfd;

/*
 * This routine is responsible for finding the input files,
 * setting up the input and output, and calling init_vplot
 * and proc_vplot. You can link to vplot without using this
 * routine, and so have your own main. See vplothacker.doc to
 * learn how to do this. (Especially you, Jon!)
 */

#ifdef SEP
int             xsepxargc;
char          **xsepxargv;
char            scrap[MAXFLEN + 1];
int             hclose_done = NO;
int             fake_header = NO;
MAIN ()

#else /* SEP */
int             sepxargc;
char          **sepxargv;	/* for getpar */
MAIN_(){} /* dummy for shared linkage */
main (argc, argv)
    int             argc;
    char           *argv[];
#endif /* SEP */
{
#ifndef SEP
int             in_isatty, num_vplot, docflag;
char            instring[MAXFLEN + 1];
#endif /* SEP */

char           *cptr;
char           *stringptr;
int             ii;
FILE           *temp;
char            string[MAXFLEN + 1];
int		tempfileindex = -1;

    nulldev ();			/* Just to make sure it gets loaded */

#ifndef SEP
    if (stringptr = strrchr(argv[0], '/'))
	strncpy (callname, ++stringptr, 24);
    else
	strncpy (callname, argv[0], 24);
#else /* SEP */
    if (stringptr = strrchr(sepxargv[0], '/'))
	strncpy (callname, ++stringptr, 24);
    else
	strncpy (callname, sepxargv[0], 24);
#endif /* SEP */

#ifdef SEP
    pltout = outstream;
    if (redout ())
    {
	getch ("head", "s", scrap);
	if (strcmp (scrap, "/dev/null") == 0 &&
	    strcmp (sepheadwhere, "/dev/null") == 0)
	{
	    fake_header = YES;
	    getch ("out", "s", scrap);
	    if (strcmp (scrap, "stdout") != 0)
	    {
		/*
		 * They probably want the header output into the redirected
		 * output. (SEP only) 
		 */
		headstream = stdout;
		headfd = fileno (headstream);
		Puthead ("Fake header for special device-dependent data only.\n");
		Puthead ("To get input history passed along, over-ride default with head = stdout\n");
		/*
		 * If the output is going into a file, then put this
		 * information into the header file. 
		 */
		if (strcmp (scrap, "/dev/tty") != 0)
		{
		    fullnm (scrap, MAXFLEN + 1);
		    Puthead ("\tin=%s\n", scrap);
		}
		else
		{
		    Puthead ("\tin= nowhere\n");
		    Puthead ("\t(Sorry, the data went to the terminal; it's gone now!)\n");
		}

		/*
		 * Have to use my own puthead routine. Standard ones such as
		 * putch, etc, don't work with this. They remember where the
		 * header used to be. Puthead does this, checking for
		 * headstream not null. This is so this will work without
		 * sep, as well. 
		 */
	    }
	}
    }

    if (!fake_header)
    {
	Puthead ("\tn3 = unknown\n\t(Sorry, SEPlib requires the header be closed now.)\n");
	Puthead ("\t(Normally you won't care that n3 is unknown unless you're using Raspen!\n)");
	hclose ();
	hclose_done = YES;
    }
#else /* SEP */

    /*
     * If no arguments, and not in a pipeline, self document "wstype="
     * doesn't count as an argument for our purposes 
     */
    in_isatty = isatty ((int) (fileno (stdin)));
    sepxargc = argc;
    sepxargv = argv;
    docflag = 0;
    if (argc == 1)
	docflag = 1;
    if ((argc == 2) && !strncmp ("wstype=", argv[1], 7))
	docflag = 1;
    getpar ("selfdoc", "1", &docflag);
    if (in_isatty && docflag)
    {
	for (ii = 0; ii < doclength; ii++)
	    printf ("%s\n", documentation[ii]);
	exit (0);
    }

    pltout = stdout;
#endif /* SEP */

#if defined(HAVE_TERMIO_H)

#else /* USG */
    /*
     * This getpar for signal is only included for debugging purposes. By
     * using a signal option, one can stop any signals from being caught. 
     */
    if (getpar ("signal", "s", string) == 0)
    {
/*#ifdef SOLARIS*/
#if defined(SOLARIS) || defined(LINUX)
        sigfillset(&(errhandler.sa_mask));
#endif
	for (ii = 0; ii < NOSIG; ++ii)
	{
#if defined(SOLARIS) || defined(LINUX)
	    if (-1 == sigaction (signum[ii], &ignored, &oldvec))
	    {
		ERR (FATAL, name, "Bad sigvec call!");
	    }
	    if (oldvec.sa_handler == ignored.sa_handler)
		(void) sigaction (signum[ii], &oldvec, (struct sigaction *) NULL);
	    else
		(void) sigaction (signum[ii], &errhandler, (struct sigaction *) NULL);
#else
	    if (-1 == sigvec (signum[ii], &ignored, &oldvec))
	    {
		ERR (FATAL, name, "Bad sigvec call!");
	    }
	    if (oldvec.sv_handler == ignored.sv_handler)
		(void) sigvec (signum[ii], &oldvec, (struct sigvec *) NULL);
	    else
		(void) sigvec (signum[ii], &errhandler, (struct sigvec *) NULL);
#endif
	}
    }
#endif /* USG */

/*
 ****************************************************************************
 * Set all global variables, open the device.
 ****************************************************************************
 */

    init_vplot ();

/*
 ****************************************************************************
 * Start processing input files
 ****************************************************************************
 */


#ifdef SEP
    if (instream != NULL)
    {
	if (infileno >= MAXIN)
	{
	    ERR (FATAL, name, "too many input files (%d max)", MAXIN);
	}
	if( cachepipe && isapipe(fileno( instream)) )
	{
	    if( (pltinarray[infileno] = tempcopy(instream,string) ) == NULL )
	    {
	    	ERR( FATAL, name, "copy of piped input failed");
	    }
	 
            /* check for zero length input (e.g. /dev/null ) */
	    if( pltinarray[infileno] != (FILE*) -1 ) {

	    strcpy( pltinname[infileno], string );
	    /*remember what number this file is so we can delete it later*/
	    tempfileindex = infileno;
	    infileno++;
	    }
	}
	else
	{
	    if (!allow_pipe && isapipe (fileno (instream)))
	    {
	    	ERR (WARN, name, "cannot use pipes with this device, try cachepipe=y ");
	    }
	    else
	    {
	        strcpy (pltinname[infileno], "Pipe");
	    	pltinarray[infileno] = instream;
	    	infileno++;
	    }
	}
    }
    else
	ERR (WARN, name, "cannot read input pipe");

    xsepxargc = sepxargc;
    xsepxargv = sepxargv;

    for (xsepxargc--, xsepxargv++; xsepxargc; xsepxargc--, xsepxargv++)
    {
	cptr = *xsepxargv;
	while (*cptr)
	{
	    if (*cptr == '=')
		break;
	    cptr++;
	}
	if (*cptr)
	    continue;
	/* Ignore dummy arguments */
	if (strcmp (*xsepxargv, "dummy") == 0)
	    continue;
	if ((temp = fopen (*xsepxargv, "r")) == NULL)
	{
	    ERR (WARN, name, "cannot open header file %s", *xsepxargv);
	    continue;
	}
	fclose (temp);
	if (getch2 ("in", "s", string, *xsepxargv))
	{
	    if ((temp = fopen (string, "r")) != NULL)
	    {
		Puthead ("   +  %s --> in = %s\n", *xsepxargv, string);
		if (infileno >= MAXIN)
		{
		    ERR (FATAL, name, "too many input files (%d max)", MAXIN);
		}
		strcpy (pltinname[infileno], string);
		pltinarray[infileno] = temp;
		infileno++;
	    }
	    else
	    {
		ERR (WARN, name, "cannot open input file %s", string);
	    }
	}
    }

#else /* SEP */
    /*
     * first process pipe input 
     */
    if (!in_isatty)
    {
	if (infileno >= MAXIN)
	{
	    ERR (FATAL, name, "too many input files (%d max)", MAXIN);
	}

  	if( cachepipe )
        {
            if( (pltinarray[infileno] = tempcopy( stdin,string) ) == NULL )
            {
                ERR( FATAL, name, "copy of piped input failed");
	    }

            /* check for zero length input (e.g. /dev/null ) */
	    if( pltinarray[infileno] != (FILE*)-1 ) {

	    strcpy( pltinname[infileno], string );
	    /* remember what number this file is so we can delete it later*/
	    tempfileindex = infileno;
            infileno++;
	    }
        }
        else
        {
	    if (!allow_pipe)
	    {
	    	ERR (WARN, name, "cannot use pipes with this device, try cachepipe=y ");
	    }
	    else
	    {
	    	strcpy (pltinname[infileno], "stdin");
	    	pltinarray[infileno] = stdin;
	    	infileno++;
	    }
	}
    }

    /*
     * next process in= inputfiles If they set num_vplot, also look for in1=
     * in2= etc 
     */

    num_vplot = 0;
    getpar ("numvplot", "d", &num_vplot);

    for (ii = 0; ii <= num_vplot; ii++)
    {
	if (ii == 0)
	    strcpy (instring, "in");
	else
	    sprintf (instring, "in%d", ii);

	if (getpar (instring, "s", string))
	{
	    if ((temp = fopen (string, "r")) != NULL)
	    {
		if (infileno >= MAXIN)
		{
		    ERR (FATAL, name, "too many input files (%d max)", MAXIN);
		}
		strcpy (pltinname[infileno], string);
		pltinarray[infileno] = temp;
		infileno++;
	    }
	    else
	    {
		ERR (WARN, name, "cannot open %s", string);
	    }
	}
    }

    /*
     * finally process input line for non-getpar arguments and assume they
     * are also input files 
     */
    for (sepxargc--, sepxargv++; sepxargc; sepxargc--, sepxargv++)
    {
	cptr = *sepxargv;
	while (*cptr)
	{
	    if (*cptr == '=')
		break;
	    cptr++;
	}
	if (*cptr)
	    continue;
	cptr = *sepxargv;
	if ((temp = fopen (cptr, "r")) != NULL)
	{
	    if (infileno >= MAXIN)
	    {
		ERR (FATAL, name, "too many input files (%d max)", MAXIN);
	    }
	    strcpy (pltinname[infileno], cptr);
	    pltinarray[infileno] = temp;
	    infileno++;
	}
	else
	{
	    ERR (WARN, name, "cannot open %s", cptr);
	}
    }
#endif /* SEP */

/*
 ****************************************************************************
 * Go do the plots
 ****************************************************************************
 */

    proc_vplot ();

#ifdef SEP
    if (!hclose_done)
    {
	Puthead ("\tn3=%d\n", nplots);
	hclose ();
	hclose_done = YES;
    }
#endif /* SEP */

    /*  delete the temporary copy of piped input if there is one*/
    removtemp();

    /* All done */
    exit (0);
    return 0;
}

#if defined(HAVE_TERMIO_H)
#else /* USG */
#ifdef SIGFNC_RTN_VOID
void cleanup ()
#else
cleanup ()
#endif
{
#ifndef SOLARIS
#ifndef LINUX
    sigblock (~(SIGKILL | SIGSTOP | SIGCONT));
#endif
#endif
    dev.close (CLOSE_INTERRUPT);
    message (MESG_ON);
    ERR (COMMENT, name, "Interrupted out.");
    dev.close (CLOSE_DONE);
    /*  delete the temporary copy of piped input if there is one*/
    removtemp();
    /*
     * Let them see what they are doing again 
     */
    if (!allowecho)
    {
/*#ifdef SOLARIS*/
#if defined(SOLARIS) || defined(LINUX)
	ioctl (pltoutfd, TCSETAW, (char *) (&tty_clean_state));
#else
	ioctl (pltoutfd, TIOCLSET, (char *) (&tty_clean_local_mode));
	ioctl (pltoutfd, TIOCSETN, (char *) (&tty_clean_state));
#endif
    }
    exit (0);
}
#endif /* USG */

/* routine to copy a file to a temporary file, used to copy stdin so that
 * it can be reread as required
 */

#define XFER_SIZE 1024

/* the name of the temporary file */
static char *tspoolnm=(char*)NULL;

removtemp(){
    if( tspoolnm != (char*)NULL )
        if( unlink(tspoolnm) )
               ERR (WARN, name, "unable to delete temporary file");
    tspoolnm = (char*)NULL;
}

FILE* tempcopy( infile, filename )
FILE* infile;
char* filename;
{
    FILE *temp;
    int len,total;
    char xfer_buf[ XFER_SIZE];
    char* spooldirnm;


    /* make up a temporary file name 
     * I know there are beter routines to make up the name 
     * but this one is more portable.
     * PEN_SPOOL is a directory we can put temporary files in
     * it is defined in params.h. It can be overridden by the
     * environment variable VPLOTSPOOLDIR.
     */
    if( (spooldirnm = getenv("VPLOTSPOOLDIR")) != NULL ){
	tspoolnm = (char *) malloc( strlen( spooldirnm ) +13 );
        strcpy( tspoolnm, spooldirnm );
    }else{
        tspoolnm = (char *) malloc( strlen( PEN_SPOOL ) + 13 );
        strcpy( tspoolnm, PEN_SPOOL );
    }
    tspoolnm = strcat( tspoolnm, "/vplotXXXXXX" );
    tspoolnm = mktemp( tspoolnm );

    if ( (temp = fopen(tspoolnm, "w")) == NULL) {
	ERR (WARN, name, "unable to create temporary file");
	return NULL;
    }

    /* now copy everything from the input stream to our temporary file */
    total=0;
    while( (len=fread( xfer_buf, sizeof(char), XFER_SIZE, infile ) ) != 0 ){
	if( ferror(infile) ) {
	    ERR (WARN, name, "read error from pipe ");
	    return NULL;
	}
	    
        if( ( fwrite(xfer_buf,sizeof(char),len,temp) != len ) || ferror(temp) ) 	{
	    ERR (WARN, name, "write error to temp ");
	    return NULL;
	}
	total += len;
    }

    fclose(temp);

    /* We could unlink the file after we reopen it to make it really
     * temporary but we don't know for sure that a filter isn't going
     * to close the file and try to reopen it.
     * Hence the buisness with tempfileindex in the main routine.
     */

    /* check for zero length input */
    if( total == 0 ) return (FILE*)-1;

    strcpy( filename, tspoolnm );
    return fopen( tspoolnm, "r");
}

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
 *  source file:   ./filters/proc_vplot.c
 *
 * Joe Dellinger (SEP), Feb 19 1988
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 *
 * Joe Dellinger Feb 22 1988
 *	Created INT_PAUSE to be separate from INT_GET_STRING.
 * Joe Dellinger Feb 27 1988
 *	Interact option turns of endpausing.
 * W. Bauske IBM 03-26-91 
 *	Apply SysV fixes for RS/6000
 * Stew Levin MOBIL 2-27-95
 *	Solaris mods
 & Bob Clapp
 *  Changed signals from BSD to POSIX1 for LINUX
 */
#include <sitedef.h>

#include	<stdio.h>
#include	<math.h>
#define		GETPAR	getpar


#if defined(HAVE_TERMIO_H)
#include	<termio.h>
#else /* USG */
#ifdef LINUX
#include	<bsd/sgtty.h>
#else
#include	<sys/ioctl.h>
#include	<sgtty.h>
#endif
#endif /* USG */
#include	<sys/types.h>
#include	<sys/stat.h>
#include	<ctype.h>
#include	<string.h>

#include	<vplot.h>

#include	"./include/params.h"	/* for machine dependencies */
#include	"./include/enum.h"
#include	"./include/err.h"
#include	"./include/attrcom.h"
#include	"./include/intcom.h"
#include	"./include/mesgcom.h"
#include	"./include/erasecom.h"
#include	"./include/closestat.h"
#include	"./include/pat.h"
#include	"./include/vertex.h"
#include	"./include/round.h"
#include	"./include/extern.h"

#if defined (HAVE_TERMIO_H)
extern struct termio tty_clean_state;
#else /* USG */
extern struct sgttyb tty_clean_state;
extern int      tty_clean_local_mode;
#endif /* USG */
extern int      need_end_erase;
extern int      buffer_input;
extern int      ever_called;
extern int      out_isatty;
extern int      nplots;
extern int      endpause;
extern int      allowecho;
extern int      epause;
extern char     interact[];
extern int      pltoutfd;
extern int      (*message) ();
extern FILE    *pltin;
extern FILE    *controltty;
extern FILE    *pltinarray[];
extern char     pltinname[][MAXFLEN + 1];
extern char     pltname[];
extern int      infileno;
extern int      (*genreader) ();

FILE           *fopen ();

/*
 * This routine is responsible for processing the input files,
 * and performing the necessary pausing, etc, that may be needed
 * at the end before exiting.
 */

proc_vplot ()
{
int             ii;
char            string[MAXFLEN + 1];

/*
 * Finally, shove all the plot files off to be done!
 */

/*	fprintf(stderr,"i am in process plot \n");*/
    if (!buffer_input)
    {
	for (ii = 0; ii < infileno; ii++)
	{
/*	fprintf(stderr,"2 am in process plot \n");*/
	    setbuf (pltinarray[ii], (char *) NULL);
/*	fprintf(stderr,"3 am in process plot \n");*/
	}
    }
/*	fprintf(stderr,"3b am in process plot \n");*/

    (*genreader) (infileno, pltinarray, pltinname);

/*	fprintf(stderr,"4 am in process plot \n");*/
/*
 * Normally, *genreader will be gen_do_dovplot, found in genlib
 */

    if (ever_called)
    {
	dev.close (CLOSE_FLUSH);
	if (epause > 0)
	{
	    sleep ((unsigned) epause);
	}
	if (need_end_erase)
	{
	    dev.erase (ERASE_END);
	}
	/*
	 * Inquire point back from device. Skip endpause stuff if we do
	 * interact, Since that's a pause in itself. 
	 */
	if (interact[0] != '\0')
	{
	    getapoint ();
	}
	else
	{

/*
 * Pause at the end if the user specifically asks us to on the command line,
 * even if we don't think we should because we think it's a file.
 */
	    if (epause <= 0 &&
		(out_isatty || getpar ("endpause", "s", string))
		&& endpause)
	    {
		dev.close (CLOSE_PAUSE);
		dev.interact (INT_F_PAUSE, controltty, string);
	    }
	}
	message (MESG_ON);
	dev.close (CLOSE_NORMAL);
	dev.close (CLOSE_DONE);
	nplots++;
    }
    else
    {
	dev.close (CLOSE_NOTHING);
	ERR (COMMENT, name, "No input?");
	dev.close (CLOSE_DONE);
    }

    /*
     * Done, let them see what they are doing again 
     */
    if (!allowecho)
    {
#if defined(HAVE_TERMIO_H)
	if (ioctl (pltoutfd, TCSETAW, &tty_clean_state) == -1)
	{
		ERR (FATAL, name, "Bad ioctl call!");
	}
#else /* USG */
	ioctl (pltoutfd, TIOCLSET, (char *) (&tty_clean_local_mode));
	ioctl (pltoutfd, TIOCSETN, (char *) (&tty_clean_state));
#endif /* USG */
    }
}

#endif
