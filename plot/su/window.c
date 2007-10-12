/* Copyright (c) Colorado School of Mines, 2007.*/
/* All rights reserved.                       */

#include <X11/Xlib.h>
#include <X11/Xutil.h>

Window
xNewWindow (Display *dpy   /* display pointer */, 
	    int x          /* x in pixels of upper left corner */, 
	    int y          /* y in pixels of upper left corner */, 
	    int width      /* width in pixels */, 
	    int height     /* height in pixels */,
	    int border     /* border pixel */, 
	    int background /* background pixel */, 
	    char *name     /* name of window (also used for icon) */)
/*< Create a new window and return the window ID 
Notes:
The parent window is the root window.
The border_width is 4 pixels. >*/
/******************************************************************************
Author:		Dave Hale, Colorado School of Mines, 01/06/90
*****************************************************************************/
{
	Window root,win;
	XSizeHints size_hints;
	XWMHints wm_hints;
	int scr,border_w=4;

	/* get screen and root window */
	scr = DefaultScreen(dpy);
	root = RootWindow(dpy,scr);

	/* create window */
	win = XCreateSimpleWindow(dpy,root,x,y,width,height,
		border_w,border,background);

	/* set window properties for window manager */
	size_hints.flags = USPosition|USSize;
	size_hints.x = x;
	size_hints.y = y;
	size_hints.width = width;
	size_hints.height = height;
	XSetStandardProperties(dpy,win,name,name,None,0,0,&size_hints);
	wm_hints.flags = InputHint;
	wm_hints.input = True;
	XSetWMHints(dpy,win,&wm_hints);

	/* return window ID */
	return win;
}
