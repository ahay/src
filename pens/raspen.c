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
 *  source file:   ./filters/raslib/rasattr.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Stewart A. Levin (MRDC/DRL), September 8, 1988
 *      Return DOVPLOT_CONT on exit.
 */

/*
 * control graphics attributes
 */
#include <sitedef.h>
#include <stdio.h>
#include "../include/attrcom.h"
#include "raspen.h"
extern FILE    *pltout;
extern int      color_mult;
extern int      esize;

int             color_table[NCOLOR][3], rascolor;

rasattr (command, value, v1, v2, v3)
    register int    command, value;
    int             v1, v2, v3;
{
    switch (command)
    {
    case SET_COLOR:
	rascolor = color_mult * value;
/*fprintf(ntderr,"setting color %d  \n",rascolor);*/
	break;
    case SET_COLOR_TABLE:
	color_table[value * color_mult][0] = v1;
	color_table[value * color_mult][1] = v2;
	color_table[value * color_mult][2] = v3;
  if(value==0) zap_change(color_table);
/*fprintf(stderr,"setting color %d %d %d \n",v1,v2,v3);*/
	break;
    default:
	break;
    }

    return (DOVPLOT_CONT);
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
 *  source file:   ./filters/raslib/rasclose.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger (SEP), May 7, 1989
 *	RGB option added.
 * Chuck Karish, August 5, 1989
 * 	For non-SEP users, don't try to open colfile unless it has been
 *	specified!
 */

#include<sitedef.h>
#include <stdio.h>
#include "../include/closestat.h"
#include "../include/err.h"
#include "../include/extern.h"
#include "../include/params.h"
#include "raspen.h"
extern int      color_table[NCOLOR][3];
extern char     colfile[];

extern int      grf_format;

rasclose (status)
    int             status;
{
int             value;
FILE           *colout;

    switch (status)
    {
    case CLOSE_NORMAL:
#ifdef SEP
	Puthead ("Color table:\n");
	for (value = 0; value < NCOLOR; value++)
	{
	    if (color_table[value][0] != -1)
	    {
		Puthead ("%d\t\t%d\t%d\t%d\n", value,
			 color_table[value][0], color_table[value][1], color_table[value][2]);
	    }
	}
	Puthead ("\n");
#endif

	if (!grf_format && esize == 1 && colfile[0])
	{
	    colout = fopen (colfile, "w");
	    if (colout == NULL)
		ERR (WARN, name, "can't open colfile %s\n", colfile);
	    else
	    {
		for (value = 0; value < NCOLOR; value++)
		{
		    if (color_table[value][0] != -1)
		    {
			fprintf (colout, "%d\t\t%f\t%f\t%f\n", value,
				 (float) color_table[value][0] / MAX_GUN,
				 (float) color_table[value][1] / MAX_GUN,
				 (float) color_table[value][2] / MAX_GUN);
		    }
		}
	    }
	}
	break;
    case CLOSE_FLUSH:
	if (!grf_format)
	    fflush (pltout);
	break;
    default:
	break;
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
 *  source file:   ./filters/raslib/rasconf.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

/*
 * Keyword: vplot raster movie pen
 */
#include<sitedef.h>
#include <stdio.h>
#include "../include/enum.h"
#include "../include/extern.h"
#include "raspen.h"

/*
 * mandatory declarations and initializations
 */
#ifdef SEP
char            name[] = "Raspen";
#else
char            name[] = "raspen";
#endif
#include "rasdoc.h"

/*
 * device routine table
 */
extern int
genmessage (), raserase (), rasopen (), genmarker ();
extern int
rasvector (), gentext ();
extern int      genpoint ();
extern int
genpatarea (), genraster ();
extern int
nulldev (), rasattr (), rasreset (), rasclose ();

struct device   dev =
{

 /* control routines */
 rasopen,		/* open */
 rasreset,		/* reset */
 genmessage,		/* message */
 raserase,		/* erase */
 rasclose,		/* close */

 /* high level output */
 rasvector,		/* vector */
 genmarker,		/* marker */
 gentext,		/* text */
 genpatarea,		/* area */
 genraster,		/* raster */
 genpoint,		/* point */
 rasattr,		/* attributes */

 /* input */
 nulldev,		/* getpoint */
 nulldev,		/* interact */

 /* low level output */
 nulldev,		/* plot */
 nulldev,		/* startpoly */
 nulldev,		/* midpoly */
 nulldev		/* endpoly */
};
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
 *  source file:   ./filters/raslib/raserase.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger (SEP), May 7, 1989
 *	RGB option added. Stew Levin's Color tektronix code modified and
 *	incorporated.
 * Joe Dellinger (SEP), August 5, 1989
 *	Fixed bug that caused only first page of GRF format output to
 *	actually print.
 * Dave Nichols (SEP) May 2 1992
 *      Allow VPLOTSPOOLDIR environment variable to override PEN_SPOOL
 * Dave Nichols (SEP) July 10 1992
 *      Added ppm output.
 */

/*
 * Erase the graphics area
 */
#include<sitedef.h>

#include <sitedef.h>
#include	<stdio.h>
#include	"../include/erasecom.h"
#include	"../include/enum.h"
#include	"../include/err.h"
#include	"../include/extern.h"
#include	"../include/params.h"
#include "raspen.h"
extern FILE    *pltout;
extern int      zap ();

#if defined(__stdc__) || defined(__STDC__)
#include <stdlib.h>
#else
extern char    *getenv ();
#endif

#ifdef HAVE_PPM
/*#include <ppm.h>*/
#include "ppm.h"
#endif

raserase (command)
    int             command;
{

    switch (command)
    {
    case ERASE_START:
    default:
	break;
    case ERASE_MIDDLE:
    case ERASE_END:
/*
 * Output raster file, and then CLEAR it out (erase)
 */
	ras_write ();
	zap ();
	break;
    case ERASE_BREAK:
/*
 * Output raster file, but don't CLEAR it!
 */
	ras_write ();
	break;
    }
}

ras_write ()
{
static int      file_created = NO;
static int      counter = 0;
char            scratch_file[120];
char            system_call[120];
int             value;
extern int      system ();
extern int      grf_format;
extern int      ppm_format;
extern int      white;
extern int      default_out;
extern int      color_table[NCOLOR][3];
extern float    o_pixels_per_inch;

    counter++;

    if (!grf_format && ! ppm_format )
    {
	if (write (fileno (pltout), image, dev_xmax * dev_ymax * esize) != dev_xmax * dev_ymax * esize)
	{
	    ERR (FATAL, name, "Can't write raster image %d\n", counter);
	}
    }
#ifdef HAVE_PPM
    else if ( ppm_format )
    {
	static int called=0;
	pixel * pixrow;
	unsigned char r,g,b;
	unsigned char *ptr;
	register pixel *pP;
	int row, col;

        if( called ) {
		return;
	}
        called=1;
	
	ppm_writeppminit( pltout, dev_xmax, dev_ymax, (pixval)255, 0);
	pixrow = ppm_allocrow( dev_xmax );
	for ( row = 0, ptr=image;  row < dev_ymax; ++row )
        {
        for ( col = 0, pP = pixrow; col < dev_xmax; ++col, ++pP ){
	    r = *(ptr++);
	    g = *(ptr++);
	    b = *(ptr++);
/*fprintf(stderr,"check this dork %d %d %d \n",(int)r,(int)g,(int)b);*/
	    PPM_ASSIGN( *pP, r, g, b );
	}
/*			fprintf(stderr,"going in %d \n",PPM_GETR(pixrow[0]));*/
/*        bbb_writeppmrow( pltout, pixrow, dev_xmax, (pixval)255, 0 );*/
/*        dpm_writeppmrow( pltout, pixrow, dev_xmax, (pixval)255, 0 );*/
        ppm_writeppmrow( pltout, pixrow, dev_xmax, (pixval)255, 0 );
        }
	pm_close( pltout );
    }
#endif
    else
    {
	if (default_out)
	{
  	    char *spooldirnm;

	    file_created = YES;
	    if( (spooldirnm = getenv("VPLOTSPOOLDIR")) != NULL ){
               sprintf (scratch_file, "%s%s%d%s", spooldirnm, "/Tcpr_", counter, "_XXXXXX");
            }else{
               sprintf (scratch_file, "%s%s%d%s", PEN_SPOOL, "/Tcpr_", counter, "_XXXXXX");
            }

	    mkstemp (scratch_file);
	    pltout = fopen (scratch_file, "w");
	    if (pltout == NULL)
	    {
		ERR (FATAL, name, "could not open scratch file %s!",
		     scratch_file);
	    }
	}
	else
	{
	    strcpy (scratch_file, "(pipe)");
	}

	/* build image attributes header block */
	fprintf (pltout, "source_platform\tvplot\n");
	fprintf (pltout, "source_application\traspen\n");
	fprintf (pltout, "image_type\tSCREEN_DATA\n");
	fprintf (pltout, "horizontal_resolution\t%d\n",
		 (int) pixels_per_inch);
	fprintf (pltout, "vertical_resolution\t%d\n",
		 (int) (pixels_per_inch / aspect_ratio));
	fprintf (pltout, "color_units\tRGB\n");
	fprintf (pltout, "number_of_scans\t%d\n", dev_ymax);
	fprintf (pltout, "scan_length\t%d\n", dev_xmax);
	fprintf (pltout, "scan_direction\tLEFT_TO_RIGHT\n");
	fprintf (pltout, "scan_order\tTOP_TO_BOTTOM\n");
	fprintf (pltout, "bits_per_value\t8\n");
	fprintf (pltout, "value_field_width\t8\n");
/*
 * The documentation seems to be inconsistent in describing whether
 * pixel_field_width should be 8 or 24 in the esize=3 case. 24 seems to
 * work.
 */
	if (esize == 1)
	    fprintf (pltout, "pixel_field_width\t8\n");
	else
	    fprintf (pltout, "pixel_field_width\t24\n");

	fprintf (pltout, "packing_level\tNO_PACKING\n");
	if (esize == 3)
	    fprintf (pltout, "map\t0\n");
	else
	{
	    fprintf (pltout, "map\t%d\n", NCOLOR);
	    for (value = 0; value < NCOLOR; value++)
	    {
		if (color_table[value][0] != -1)
		{
		    fprintf (pltout, "%d\t%d\t%d\n", color_table[value][0],
			     color_table[value][1], color_table[value][2]);
		}
		else
		{
		    fprintf (pltout, "%d\t%d\t%d\n", value, value, value);
		}
	    }
	}
	if (white)
	{
	    fprintf (pltout, "invert_bw\tINVERT\n");
	}
	else
	{
	    fprintf (pltout, "invert_bw\tPASS\n");
	}

	fprintf (pltout, "header_end\n");
	fflush (pltout);

	/* send image data */
	if (write (fileno (pltout), image, dev_xmax * dev_ymax * esize) != dev_xmax * dev_ymax * esize)
	{
	    ERR (FATAL, name, "Can't write raster image #%d to \"%s\".\n", counter, scratch_file);
	}
	fflush (pltout);

	if (file_created)
	{
	    fclose (pltout);
	    sprintf (system_call, "lpr -Ptcpr -r %s\n", scratch_file);
	    if (0 != system (system_call))
	    {
		/* else figure we'll need to delete scratch file ourselves */
		ERR (WARN, name, "Couldn't get lpr to cooperate\n");
		unlink (scratch_file);
	    }
	    file_created = NO;
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
 *  source file:   ./filters/raslib/rasopen.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger (SEP), Sept 25 1988
 *	"ppi" undocumented, and read as an integer and not a float!
 * Joe Dellinger (SEP), May 7 1989
 *	Cleaned up to make raslib work with Tektronix color plotter.
 *	Incorporated portions of Stew Levin's color plotter code.
 * Joe Dellinger (SEP), August 5, 1989
 *	Use the variable "default_out" to remember the ORIGINAL
 *	value of isatty(fileno(pltout)).
 * Chuck Karish, 5 August 1989
 *	Non-SEP, colfile defaults to empty string instead of "colfile".
 * Joe Dellinger, 17 Jan 1990
 *	ARG!!! I'm changing the default behavior back so it matches
 *	the documentation. I'm also fixing the bug introduced by the
 *	change so that "colfile" was not even settable from the command
 *	line.
 * Dave Nichols, 10 July 1992
 *	Added support for ppm output, (ppmpen, Ppmpen).
 */

#include<sitedef.h>

#include <stdio.h>
#include "../include/enum.h"
#include "../include/extern.h"
#include "../include/err.h"
#include "../include/params.h"
#if defined(__stdc__) || defined(__STDC__)
#include <string.h>
#else
#include <strings.h>
#endif
#ifdef SEP
#include <ctype.h>
#endif
#define DEFAULT_OUT	isatty(fileno(pltout))
#include "raspen.h"

/*#include <ppm.h>*/
#ifdef HAVE_PPM
#include "ppm.h"
#endif
extern int sepxargc;
extern char** xragv;

unsigned char  *image;
extern float    aspect_ratio;
extern float    pixels_per_inch;
int             color_mult, rasor = 0;
char            colfile[60];
int             esize = 1;

int             grf_format = NO;
int             ppm_format = NO;
int             white = NO;
float           o_pixels_per_inch, o_aspect_ratio;
int             default_out = YES;

extern int      num_col;
extern int      rasvector2 ();
extern int      rasvector3 ();
extern int      rasvector4 ();

rasopen ()
{
extern char   **sepxargv;
extern FILE    *pltout;
extern char    *alloc ();
#ifdef SEP
extern int      headfd;
char            headname[30], fname[80];
char            path[50];
char           *front, *tail, *ptr;
#endif
char            newpath[60];

    txfont = DEFAULT_HARDCOPY_FONT;
    txprec = DEFAULT_HARDCOPY_PREC;

    allowecho = YES;

    if (
	strcmp (callname, "tcprpen") == 0 ||
	strcmp (callname, "Tcprpen") == 0 ||
	strcmp (wstype, "tcpr") == 0
     )
    {
	grf_format = YES;
	getpar ("white", "1", &white);
    }

    if ( 
	strcmp (callname, "ppmpen" ) == 0 ||    
	strcmp (callname, "Ppmpen" ) == 0 ||    
	strcmp (wstype, "ppm") == 0
     )
    {
#ifndef HAVE_PPM
	ERR(FATAL,name," ppm format not available ");
#else
	ppm_format = YES;
	ppm_init( &sepxargc, sepxargv );
#endif
    }


    if (grf_format)
    {
/*
 * physical device parameters for Tektronix color printer 4693D
 */
	o_pixels_per_inch = pixels_per_inch = 300.;
	dev_xmax = 8.53 * pixels_per_inch;
	dev_ymax = 8.13 * pixels_per_inch;
	dev_xmin = 0;
	dev_ymin = 0;
	o_aspect_ratio = aspect_ratio = 1.;
/*
 * Best compromise between resolution and speed
 */
	pixels_per_inch = 150.;

    }else if (ppm_format)
    {
/*
 * physical device parameters for ppm format output
 */
	o_pixels_per_inch = pixels_per_inch = 100.;
	dev_xmax = 10 * pixels_per_inch;
	dev_ymax = 7.5 * pixels_per_inch;
	dev_xmin = 0;
	dev_ymin = 0;
	o_aspect_ratio = aspect_ratio = 1.;
    }
    else
    {
/*
 * physical device parameters for RasterTek
 */
	o_pixels_per_inch = pixels_per_inch = 100.;
	dev_xmax = 1024;
	dev_ymax = 768;
	dev_xmin = 0;
	dev_ymin = 0;
	o_aspect_ratio = aspect_ratio = 1.;
    }

/*
 * device capabilities
 */
    need_end_erase = YES;
    buffer_output = YES;
    smart_clip = NO;

    if (grf_format)
    {
	getpar ("esize", "d", &esize);
	if (esize == 3)
	{
	    dev.vector = rasvector4;
	}
	else
	if (esize != 1)
	    ERR (FATAL, name, "Esize must be 1 or 3!\n");

	color_mult = 1;
	num_col = NCOLOR;
    }
    else if ( ppm_format )
    {
	/* always use esize 3 */
	esize=3;
	color_mult =1;
	num_col = NCOLOR;
	dev.vector = rasvector3;

    }
    else
    {
	getpar ("esize", "d", &esize);
	if (esize == 1)
	{
	    color_mult = 2;
	    getpar ("colormult", "d", &color_mult);
	    num_col = NCOLOR / color_mult;

	    getpar ("or", "1", &rasor);
	    if (rasor)
	    {
		dev.vector = rasvector2;
	    }
	}
	else
	if (esize == 3)
	{
	    color_mult = 1;
	    num_col = NCOLOR;
	    dev.vector = rasvector3;
	}
	else
	    ERR (FATAL, name, "Esize must be 1 or 3!\n");
    }


    getpar ("aspect", "f", &aspect_ratio);
    getpar ("ppi", "f", &pixels_per_inch);
    dev_xmax *= pixels_per_inch / o_pixels_per_inch;
    dev_ymax *= (o_aspect_ratio / aspect_ratio) *
     (pixels_per_inch / o_pixels_per_inch);
    getpar ("n1", "d", &dev_xmax);
    getpar ("n2", "d", &dev_ymax);
#ifdef SEP
    if (esize == 1)
    {
	Puthead ("\n\n# Raspen: VPLOT graphics via Movie,\n");
	Puthead ("#\tor any other byte-deep raster device.\n\n");
    }
    else
    {
	Puthead ("\n\n# Raspen: VPLOT graphics via RGB byte triples.\n\n");
    }
    Puthead ("\taspect_ratio=%f\n", aspect_ratio);
    Puthead ("\tppi=%f\n", pixels_per_inch);
    Puthead ("\tesize=%d\n", esize);
    Puthead ("\tn1=%d\n", dev_xmax);
    Puthead ("\tn2=%d\n", dev_ymax);
#endif


    /*
     * Allocate space for image 
     */
    if ((image = (unsigned char *) malloc (dev_xmax * dev_ymax * esize)) == NULL)
    {
	ERR (FATAL, name, "Can't allocate space for raster image\n");
    }

    default_out = DEFAULT_OUT;

    if (!grf_format && default_out)
    {
#ifdef SEP
	datapath (path);
/* Code stolen from output.c to get a reasonable raster file name. */
	if (0 < findnm (headfd, headname, sizeof (headname)))
	{
	    /* modify slightly */
	    strcpy (fname, "");
	    front = strrchr (headname, '/');
	    if (front == ((char *) NULL))
		front = headname;
	    else
		front++;
	    if ((*front) == 'H')
		strcat (fname, ++front);
	    else
	    {
		tail = strrchr (front, '.');
		if (tail == ((char *) NULL))
		    strcat (fname, front);
		else
		{
		    for (ptr = tail + 1; *ptr; ptr++)
			if (!isupper (*ptr))
			    break;
		    if (!(*ptr))/* drop suffix if all caps */
			*tail = '\0';
		    (void) strcat (fname, front);
		}
	    }
	    (void) strcat (fname, ".raster");
	}
	else
	{
	    strcpy (fname, "raster");
	}

	sprintf (newpath, "%s%s", path, fname);
	Puthead ("\tin=%s\n", newpath);
#else
	sprintf (newpath, "%s", "raster_file");
#endif
	pltout = fopen (newpath, "w");
	if (pltout == NULL)
	    ERR (FATAL, name, "can't open file %s\n", newpath);
    }

    if (!grf_format &&  !ppm_format  && esize == 1)
    {
	strcpy (colfile, "colfile");
	getpar ("colfile", "s", colfile);
#ifdef SEP
	Puthead ("\tcolfile=%s\n", colfile);
	Puthead ("\tcolor=T\n");
	if (color_mult == 1)
	    Puthead ("\tumask=255\n");
#endif
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
 *  source file:   ./filters/raslib/rasreset.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */
#include<stdio.h>
#include<sitedef.h>


#include "raspen.h"
#include "../include/params.h"
extern int      color_table[NCOLOR][3], rascolor;
extern int      color_mult;

rasreset ()
{
int             value;
    zap ();
    for (value = 0; value < NCOLOR; value++)
    {
/*why should i bother resetting this in 99% cases-RGC*/
/*	color_table[value][0] = -1;*/
    }
    for (value = 0; value < 8; value++)
    {
	color_table[value * color_mult][0] = MAX_GUN * ((value & 2) / 2);
	color_table[value * color_mult][1] = MAX_GUN * ((value & 4) / 4);
	color_table[value * color_mult][2] = MAX_GUN * ((value & 1) / 1);
    }
	color_table[0 * color_mult][0] =256;
	color_table[0 * color_mult][1] =0;
	color_table[0 * color_mult][2] =0;
}

/*
 * Zero image
 */
zap ()
{
unsigned char  *p;
    for (p = image; p < &image[dev_xmax * dev_ymax * esize]; p++)
	*p = 0;
}
zap_change(int *col)
{
unsigned char  *p;
  if(esize==3){
    for (p = image; p < &image[dev_xmax * dev_ymax*esize];){
	*p++ =(unsigned char) col[0];
	*p++ =(unsigned char)col[1];
	*p++ =(unsigned char)col[2];
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
 *  source file:   ./filters/raslib/rasvector2.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#define  OR
#include "rasvector.c"
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
 *  source file:   ./filters/raslib/rasvector3.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#define  RGB
#include "rasvector.c"
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
 *  source file:   ./filters/raslib/rasvector4.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#define  RGB
/*
 * The tektronix color plotter uses the order {Blue, Green, Red}, despite
 * what the documentation may say! Sigh.
 */
#define  BGR
#include "rasvector.c"
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
 *  source file:   ./filters/raslib/rasvector.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger (SEP), May 7 1989
 *	RGB version added.
 */
#include<sitedef.h>


#include <stdio.h>
#include <math.h>
#include "../include/extern.h"
#include "raspen.h"

/*
 * This code originally written by Jeff Thorson ("apenh"),
 * stolen by Joe Dellinger and somewhat modified for vplot standardization.
 */

#ifdef RGB
#ifndef BGR
extern int      color_table[NCOLOR][3], color_mult;
#define WRITEIT(A,B) \
Image3 (A, B, 0) = color_table[rascolor * color_mult][0]; \
Image3 (A, B, 1) = color_table[rascolor * color_mult][1]; \
Image3 (A, B, 2) = color_table[rascolor * color_mult][2]
#define RASVECTOR rasvector3
#else
extern int      color_table[NCOLOR][3], color_mult;
#define WRITEIT(A,B) \
Image3 (A, B, 0) = color_table[rascolor * color_mult][2]; \
Image3 (A, B, 1) = color_table[rascolor * color_mult][1]; \
Image3 (A, B, 2) = color_table[rascolor * color_mult][0]
#define RASVECTOR rasvector4
#endif
#else
#ifdef OR
#define WRITEIT(A,B) Image (A, B) |= rascolor
#define RASVECTOR rasvector2
#else
#define WRITEIT(A,B) Image (A, B) = rascolor
#define RASVECTOR rasvector
#endif
#endif

RASVECTOR (x1, y1, x2, y2, nfat, dashon)
    int             x1, y1, x2, y2;
    int             nfat, dashon;
{
int             test, tmp, x, y;
double          slope, fx, fx3, fy, fy3;

/*
 * Vector rasterizes the line defined by the endpoints (x1,y1) and (x2,y2).
 * If 'nfat' is nonzero then draw parallel lines to fatten the line, by
 * recursive calls to vector.
 */

    if (nfat < 0)
	return;

    if (dashon)
    {
	dashvec (x1, y1, x2, y2, nfat, dashon);
	return;
    }

    if (nfat)
    {
	if (clip (&x1, &y1, &x2, &y2))
	    return;

	fatvec (x1, y1, x2, y2, nfat, dashon);
	return;
    }

    if (clip (&x1, &y1, &x2, &y2))
	return;

/* Beware checks out of bounds, since the coordinate system may have rotated */

    test = (abs (x2 - x1) >= abs (y2 - y1));

    if (test)
    {
	if (x1 == x2)
	{
	    /* Just a point */
	    WRITEIT (x1, y1);
	    return;
	}
	else
	if (x1 > x2)
	{
	    tmp = x1;
	    x1 = x2;
	    x2 = tmp;
	    tmp = y1;
	    y1 = y2;
	    y2 = tmp;
	}
	slope = (double) (y2 - y1) / (double) (x2 - x1);
	fy3 = y1;

	for (x = x1, fy = fy3; x < x2; x++, fy += slope)
	{
	    y = fy + .5;	/* OK rounding, since always positive */
	    WRITEIT (x, y);
	}
	WRITEIT (x2, y2);
	return;
    }
    else
    {
	/* y1 can't equal y2 here */
	if (y1 > y2)
	{
	    tmp = x1;
	    x1 = x2;
	    x2 = tmp;
	    tmp = y1;
	    y1 = y2;
	    y2 = tmp;
	}
	slope = (double) (x2 - x1) / (double) (y2 - y1);
	fx3 = x1;

	for (y = y1, fx = fx3; y < y2; y++, fx += slope)
	{
	    x = fx + .5;
	    WRITEIT (x, y);
	}
	WRITEIT (x2, y2);
	return;
    }
}
