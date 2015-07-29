/* Vplot filter for X-Window */

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
 *  source file:   ./filters/x11penlib/x11penconf.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Stewart A. Levin (DRL), September 21, 1988
 *      Cribbed Xpen version
 * Dave Nichols (SEP), July 27, 1989
 *      Now uses new xraster routine, and corrected xarea .
 */

/*
 *
 *  source file:   ./filters/x11penlib/xarea.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Stewart A. Levin (DRL), September 21, 1988
 *	Cribbed xpen for X-11 Windows version.
 * Dave Nichols (SEP), July 27, 1989
 *	Corrected calling sequence for XFillPolygon
 */

/*
 *
 *  source file:   ./filters/x11penlib/xattr.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Stewart A. Levin (DRL), September 21, 1988
 *	Cribbed Xpen for XU01 Windows version.
 */

/*
 *
 *  source file:   ./filters/x11penlib/xclose.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Stewart A. Levin (DRL), September 21, 1988
 *      Cribbed xpen to support X-11 windows
 * Stewart A. Levin (DRL), December 30, 1988
 *      Exit on button press only, ignore keyboard.  This makes it easier
 *      to edit in another window without accidentally blowing plot away.
 */

/*
 *
 *  source file:   ./filters/x11penlib/xerase.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Stewart A. Levin (DRL), September 11, 1988
 *	Cribbed Xpen for X-11 Windows version.
 * Stewart A. Levin (Mobil), May 8, 1996
 *	Use R5 treatment of GC's.
 */

/*
 *
 *  source file:   ./filters/x11penlib/xgetpoint.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Stewart A. Levin (DRL), September 21, 1988
 *	Cribbed Xpen for X-11 Windows version.
 */

/*
 *
 *  source file:   ./filters/x11penlib/xopen.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Stewart A. Levin (DRL), September 21, 1988
 *      Cribbed xpen to support X-11 windows
 * Stewart A. Levin (DRL), January 16, 1989
 *	Added recognition of server subscreen notation
 * Steve Cole (SEP), July 6, 1989
 *      Added wait for exposure event to make it work with window manager.
 * Dave Nichols (SEP), July 27, 1989
 *      Changed colormap handling, if it can it gets its own colormap.
 */

/*
 *
 *  source file:   ./filters/x11penlib/xplot.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Stewart A. Levin (DRL), September 21, 1988
 *      Cribbed xpen to support X-11 Windows.
 */

/*
 *
 *  source file:   ./filters/x11penlib/xpoint.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Stewart A. Levin (DRL), September 21, 1988
 *	Cribbed Xpen for X-11 Windows.
 */

/*
 *
 *  source file:   ./filters/x11penlib/xraster.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Stewart A. Levin (DRL), September 21, 1988
 *	Cribbed Xpen for X-11 Windows version.
 * Dave Nichols (SEP), July 27, 1989
 *	Total rewrite, should support mono case and all orientations,
 *	and cases where we use a pulic colormap (e.g. DirectColor visuals).
 *      It is still slow(but more general), so go buy a faster machine.
 */

#include <string.h>
#include <stdio.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include <rsfplot.h>

#include "../include/extern.h"
#include "../include/attrcom.h"
#include "../include/closestat.h"
#include "../include/mesgcom.h"
#include "../include/erasecom.h"
#include "../include/enum.h"
#include "../include/params.h"

#include "dovplot.h"
#include "init_vplot.h"

#include "../utilities/util.h"
#include "../genlib/genpen.h"

#include "x11doc.h"
#include "x11pen.h"

#include "_x11.h"

char name[] = "x11pen";

#include "_device.h"

/* physical device size */
int             version = 1;

#define MAXVERT 1000
XPoint vlist[MAXVERT];

void xarea (int npts, struct vertex  *head)
/*< area >*/
{
    int             i;

    /* translate data structures */
    for (i = 0; i < npts && i < MAXVERT; i++)
    {
	vlist[i].x = XCORD (head->x);
	vlist[i].y = YCORD (head->y);
	head = head->next;
    }
    XDrawLines(pen_display, pen_win, pen_gc, vlist, npts, CoordModeOrigin);
    XFillPolygon(pen_display, pen_win, pen_gc, vlist, npts, Complex, CoordModeOrigin);
}


/* we cannot control pixel values for public color table entries in X,
   so mapping is used */
/* colors 0, 1 for windows */

int             color;	/* current color */
static int imap = 2; /* allocated color count */
#define U0 (~((unsigned long) 0))
unsigned long  map[256] = { U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0,
			    U0, U0, U0, U0, U0, U0, U0, U0 };

void xattributes (int command, int v, int v1, int v2, int v3)
/*< attributes >*/
{
    XColor           def;

    switch (command)
    {
	case SET_COLOR:
	    if (!dev.num_col)
	    {
		if (v)
		    color = 1;
		else
		    color = 0;
	    }
	    else
	    {
		if (v < dev.num_col)
		    color = map[v];
	    }
	    XSetForeground(pen_display, pen_gc, color);
	    break;
	case SET_COLOR_TABLE:
	    if(map[v] != U0) {/* dealloc that color first if using public colormap*/
		if( !own_colormap ){
		    XFreeColors(pen_display, pen_colormap, map+v, 1, AllPlanes); 
		}
		imap--;
	    }
	    
	    if (++imap >= dev.num_col)
		break;
	    def.red = v1 << 8;
	    def.green = v2 << 8;
	    def.blue = v3 << 8;
	    def.pixel = v;
	    def.flags = DoRed | DoGreen | DoBlue;
	    
	    /* put the color in the colormap, we can specify the pixel value if we
	       own the colormap, otherwise we must accept what pixel we are given */
	    if( own_colormap ){
		def.pixel = v;
		XStoreColor( pen_display, pen_colormap, &def );
		map[v] = def.pixel;
	    }else{
		
	 	if( !XAllocColor(pen_display, pen_colormap, &def)){
		    map[v]=1;
		}else{
		    map[v] = def.pixel;
		}
	    }
	    
	    
	    
	    if(v == 0) { /* change of background color */
		XSetBackground(pen_display, pen_gc, def.pixel);
		XSetWindowBackground(pen_display, pen_win, def.pixel);
		XClearWindow(pen_display, pen_win);
	    }
	    break;
    }
}


void xclose (int status)
/*< close >*/
{
    XEvent event;
    int             /* fd ,*/ pid;
    char           *option;
    bool persist;

    switch (status)
    {
	
	case CLOSE_FLUSH:
	    XFlush (pen_display);
	    break;
	    
	case CLOSE_DONE:
	    
	    if ((option = XGetDefault (pen_display, name, "Stay")) != NULL)
		if ((option[0] == 'y') || (option[0] == 'Y') || (option[0] == '1'))
		    persist = true;
	    if (!sf_getbool("stay",&persist)) persist=false;
	    /* open terminal to count keys */

	    /* fd = open ("/dev/tty", 0); */
	    if (!persist)
	    {
		dev.message (MESG_TEXT, "Click mouse in window to end.\n");
		XSync (pen_display, True);
		/* die on mouse click */
		XSelectInput (pen_display, pen_win, ButtonPressMask);
		/* wait for mouse click in display */
		XNextEvent(pen_display, &event);
	    }
	    else
	    {
		pid = fork ();
		if (pid < 0)
		{
		    fprintf (stderr, "x11pen: xclose: Couldn't fork a child\n");
		    XCloseDisplay (pen_display);
		    exit (-1);
		}
		else
		    if (pid > 0)
		    {
			_exit (0);
		    }
		    else
			if (pid == 0)
			{
			    dev.message (MESG_TEXT, "Click mouse in vplot window to end.\n");
			    /*	XSync (pen_display, True); */
			    /* die on key or mouse click */
			    XSelectInput (pen_display, pen_win, ButtonPressMask);
			    XNextEvent(pen_display, &event);
			    XSelectInput (pen_display, pen_win, 0);
			}
	    }
	    
	    XCloseDisplay (pen_display);
    }
}

void xerase (int command)
/*< erase >*/
{
    unsigned long   erase_color;
    unsigned long	prev_color;
    XGCValues look_em_up;

    XGetGCValues(pen_display, pen_gc, GCForeground|GCBackground,&look_em_up);
    
    switch (command)
    {
	case ERASE_START:
	case ERASE_MIDDLE:
	case ERASE_BREAK:
	    if (!dev.num_col)
		erase_color = 0;
	    else
		erase_color = look_em_up.background;
	    prev_color = look_em_up.foreground;
	    XSetForeground(pen_display, pen_gc, erase_color);
	    XFillRectangle (pen_display, pen_win, pen_gc, 0, 0,
			    dev.xmax, dev.ymax);
	    XSetForeground(pen_display, pen_gc, prev_color);
	    break;
    }
}

int xgetpoint (int *x, int *y)
/*< getpoint >*/
{
    XEvent event;

    XSync (pen_display, True);
    XSelectInput (pen_display, pen_win, ButtonPressMask);
    XNextEvent(pen_display,&event);
    XSelectInput (pen_display, pen_win, 0);
    *x = XCORD (event.xbutton.x);
    *y = YCORD (event.xbutton.y);
    return (0);
}

/* pen colors */
/* black, blue, red, magenta, green, cyan, yellow, white */
unsigned char   red[256] = {0, 0, 255, 255, 0, 0, 255, 255};
unsigned char   green[256] = {0, 0, 0, 0, 255, 255, 255, 255};
unsigned char   blue[256] = {0, 255, 0, 255, 0, 255, 0, 255};

Display        *pen_display;
int		pen_screen;
Window          pen_win;
GC		pen_gc;
Visual 	       *pen_visual;
Colormap	pen_colormap;
int		own_colormap;
float           win_xfrac = .77;
float           win_yfrac = .56;
char           *align;

void opendev (int argc, char* argv[])
/*< open >*/
{
    int             win_x, win_y, i;
    char            title[50], *ap, *option;
    Display        *OpenDisplay ();
    XEvent		event;

    dev.erase = xerase;
    dev.close = xclose;
    dev.area = xarea;
    dev.raster = xraster;
    dev.point = xpoint;
    dev.attributes = xattributes;
    dev.getpoint = xgetpoint;
    dev.plot = xplot;

    dev.xmin = 0;
    dev.ymin = 0;
    dev.pixels_per_inch = 72.;
    dev.aspect_ratio = 1.;

    /* device capabilities */
    dev.need_end_erase = true;

    /* find a display server */
    pen_display = OpenDisplay ();
    if (pen_display == NULL)
	sf_error("Can't find a display");
    pen_screen = DefaultScreen(pen_display);
    
    /* pen window dimensions */
    if ((option = XGetDefault (pen_display, name, "Width")) != NULL)
	win_xfrac = (float) atof (option);
    if ((option = XGetDefault (pen_display, name, "Height")) != NULL)
	win_yfrac = (float) atof (option);
    sf_getfloat("width",&win_xfrac);
    sf_getfloat ("height",&win_yfrac);
    dev.xmax = win_xfrac * (DisplayWidth (pen_display,pen_screen) - 2 * BORDER);
    dev.ymax = win_yfrac * (DisplayHeight (pen_display,pen_screen) - TBORDER - BORDER);

    /* pen window alignment- win_x & win_y */
    align = sf_getstring("align");
    if (NULL == align) {
	align = sf_charalloc(5);

	if ((option = XGetDefault (pen_display, name, "Align")) != NULL) {
	    strcpy (align, option);
	} else {
	    strcpy (align, "tr");
	}
    }
    
    win_x = DisplayWidth (pen_display,pen_screen) - dev.xmax - 2 * BORDER;	/* right default */
    win_y = 0;			/* top default */
    for (ap = align; *ap != '\0'; ap++)
	switch (*ap)
	{
	    case 't':		/* top */
		win_y = 0;
		break;
	    case 'm':		/* middle vertical */
		win_y = (DisplayHeight (pen_display,pen_screen) - dev.ymax - TBORDER - BORDER) / 2;
		break;
	    case 'b':		/* bottom */
		win_y = DisplayHeight (pen_display,pen_screen) - dev.ymax - TBORDER - BORDER;
		break;
	    case 'l':		/* left */
		win_x = 0;
		break;
	    case 'c':		/* center horizontal */
		win_x = (DisplayWidth (pen_display,pen_screen) - dev.xmax - 2 * BORDER) / 2;
		break;
	    case 'r':		/* right horizontal */
		win_x = DisplayWidth (pen_display,pen_screen) - dev.xmax - 2 * BORDER;
		break;
	}
    
    pen_visual = DefaultVisual(pen_display, pen_screen);
    
    /* color support */
    dev.num_col = DisplayCells(pen_display, pen_screen);
    if( dev.num_col <= 2 ) mono=true;

#if defined(__cplusplus) || defined(c_plusplus)
    if ( pen_visual->c_class == PseudoColor ){
#else
    if ( pen_visual->class == PseudoColor ){
#endif	
	pen_colormap = XCreateColormap(pen_display,
				       RootWindow(pen_display, pen_screen), pen_visual,  AllocAll);
	own_colormap = 1;
    }else{
        pen_colormap = DefaultColormap(pen_display, pen_screen);
	own_colormap = 0;
    }
    
    
    /* spawn window */
    sprintf (title, "XPEN %d x %d", dev.xmax, dev.ymax);
    pen_win = XCreateSimpleWindow (pen_display, RootWindow(pen_display,pen_screen),
				   win_x, win_y, dev.xmax, dev.ymax, BORDER, None, None);
    XSetWindowBackground(pen_display, pen_win,
			 BlackPixel(pen_display,pen_screen));
    XSetWindowColormap( pen_display, pen_win, pen_colormap );
    XStoreName(pen_display, pen_win, title);
    XMapWindow (pen_display, pen_win);
    
    
    /* graphics context */
    pen_gc = XCreateGC(pen_display, pen_win, 0, NULL);
    XSetLineAttributes(pen_display,pen_gc, 0, LineSolid, CapButt, JoinRound);
    XSetPlaneMask(pen_display, pen_gc, AllPlanes);
    XSetFunction(pen_display, pen_gc, GXcopy);
    XSetFillRule(pen_display, pen_gc, WindingRule);
    XSetBackground(pen_display,pen_gc,BlackPixel(pen_display, pen_screen));
    
    /* color support */
    xattributes (SET_COLOR, 1, 0, 0, 0);
    
    if( dev.num_col > 2){
    	for (i = 0; i < 8; i++)
	    xattributes (SET_COLOR_TABLE, i, red[i], green[i], blue[i]);
    }
    
    /* ignore all events for now except expose events*/
    XSelectInput(pen_display,pen_win,ExposureMask);
    
    /* wait for exposure event */
    XNextEvent(pen_display,&event);
    switch(event.type) {
	case Expose:
	    break;
	default:
	    break;
    }
}

/* routine to find display server; check getpar */
Display  *OpenDisplay ()
{
    char            name[50], *host;
    Display		*display;
    
    display = NULL;
    
    /* check getpar first */
    host = sf_getstring("server"); /* X server */
    if (NULL != host)
    {
	display = XOpenDisplay (host);
	if (display != NULL)
	    sf_warning ("X-server=%s", DisplayString(display));
	else			/* add :0 to getpar name */
	{
	    sprintf (name, "%s:0", host);
	    display = XOpenDisplay (name);
	    if (display != NULL)
		fprintf (stderr, "X-server=%s\n", DisplayString(display));
	}
	return (display);
    }

    display = XOpenDisplay(NULL);
    if (display != NULL)
	sf_warning ("X-server=%s", DisplayString(display));
    return(display);
}

void xplot (int x, int y, int draw)
/*< plot >*/
{
    static int oldx = 0, oldy = 0;

    if (draw)
    {
	XDrawLine (pen_display, pen_win, pen_gc,
		   XCORD (oldx), YCORD (oldy), XCORD (x), YCORD (y));
    }
    oldx = x;
    oldy = y;
}

void xpoint (int x, int y)
/*< point >*/
{
    XDrawPoint (pen_display, pen_win, pen_gc, XCORD (x), YCORD (y));
}

void xraster (int count, int out_of, int xpos, int ypos, int length, int orient, 
	      unsigned char **raster, int dummy1, int dummy2)
/*< raster >*/
{
    XImage	*rasrect;
    char *data;
    int ii,h,w,xxpos,yypos;
    int xpix,xinc,ypix,yinc;
    
    switch ( orient ) { 
	case 0 :
	default:
	    xxpos = xpos; yypos = ypos; w = length; h = 1;
	    xpix =0 ; xinc =1; ypix=0, yinc=0;
	    break;
	case 1 :
	    xxpos = xpos; yypos = ypos; w = 1; h = length;
	    xpix =0 ; xinc =0; ypix=0, yinc=1;
	    break;
	case 2 :
	    xxpos = xpos-length ; yypos = ypos; w = length; h = 1;
	    xpix =length-1 ; xinc = -1; ypix=0, yinc=0;
	    break;
	case 3 :
	    xxpos = xpos; yypos = ypos+length; w = 1; h = length;
	    xpix =0 ; xinc =0; ypix=length-1, yinc= -1;
	    break;
    }
    
    if( mono ){
    	data = sf_charalloc( (w+7)/8*h ); 
    	rasrect = XCreateImage( pen_display, pen_visual,
				1, ZPixmap, 0,data, w, h, 8, 0 );
    }else{
    	data = sf_charalloc( w*h ); 
    	rasrect = XCreateImage( pen_display, pen_visual,
				8, ZPixmap, 0,data, w, h, 8, 0 );
    }
    
    if( own_colormap || mono ){
	
	/* raster values really are pixels if we have our own colormap or if
	   it is a mono display */
	
	for( ii=0 ; ii <length ; ii++ ){
	    XPutPixel( rasrect, xpix, ypix, (unsigned long)raster[0][ii] );
	    xpix += xinc; ypix += yinc;
	}
    }else{
	
	/* if we are using a public colormap we need to translate vplot-raster
	   values into X-pixel values using the map array */
	
	for( ii=0 ; ii <length ; ii++ ){
	    XPutPixel( rasrect, xpix, ypix, (unsigned long)map[raster[0][ii]] );
	    xpix += xinc; ypix += yinc;
	}
    }
    
    XPutImage(pen_display,pen_win,pen_gc,rasrect,0,0,
	      XCORD(xxpos),YCORD(yypos),w,h);
    XDestroyImage(rasrect);
}
