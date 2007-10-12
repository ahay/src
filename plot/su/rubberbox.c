/* Copyright (c) Colorado School of Mines, 2007.*/
/* All rights reserved.                       */

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "colormap.h"

void 
xRubberBox (Display *dpy /* display pointer */, 
	    Window win   /* window ID */, 
	    XEvent event /* event of type ButtonPress */,
	    int *x       /* x of upper left hand corner of box in pixels */, 
	    int *y       /* y of upper left hand corner of box in pixels */, 
	    int *width   /* width of box in pixels */, 
	    int *height  /* height of box in pixels */)
/*< Track pointer with rubber box 

xRubberBox assumes that event is a ButtonPress event for the 1st button;
i.e., it tracks motion of the pointer while the 1st button is down, and
it sets x, y, w, and h and returns after a ButtonRelease event for the
1st button.

Before calling xRubberBox, both ButtonRelease and Button1Motion events 
must be enabled. >*/
/******************************************************************************
Author:		Dave Hale, Colorado School of Mines, 01/27/90
*****************************************************************************/
{
	GC gc;
	XGCValues *values=NULL;
	XEvent eventb;
	XStandardColormap scmap;
	int scr=DefaultScreen(dpy);
	int xb,yb,w,h,x1,x2,y1,y2,xorig,yorig,xold,yold;
	unsigned long background;

	/* determine typical background color */
	/* +1 added by John Stockwell 23 Jun 1993 */
	/* to shift xwigb rubberbox from light green to red */
	if (xCreateRGBDefaultMap(dpy,&scmap))
		background = (xGetFirstPixel(dpy)+xGetLastPixel(dpy) + 1)/2;
	else
		background = WhitePixel(dpy,scr);


	/* make graphics context */
	gc = XCreateGC(dpy,win,0,values);
  	XSetFunction(dpy,gc,GXxor);
  	XSetForeground(dpy,gc,BlackPixel(dpy,scr)^background);

	/* track pointer */
	xorig = event.xbutton.x;
	yorig = event.xbutton.y;
	xold = xorig;
	yold = yorig;
	x1 = xorig;
	y1 = yorig;
	w = 0;
	h = 0;
	while(h|(~h)/*True*/) {
		XNextEvent(dpy,&eventb);
		if (eventb.type==ButtonRelease) {
			xb = eventb.xbutton.x;
			yb = eventb.xbutton.y;
			break;
		} else if (eventb.type==MotionNotify) {
			xb = eventb.xmotion.x;
			yb = eventb.xmotion.y;

			/* if box is the same, continue */
			if (xb==xold && yb==yold) 
				continue;

			/* erase old box */
			x1 = (xold<xorig)?xold:xorig;
			y1 = (yold<yorig)?yold:yorig;
			x2 = (xold>xorig)?xold:xorig;
			y2 = (yold>yorig)?yold:yorig;
			w = x2-x1;
			h = y2-y1;
			XDrawRectangle(dpy,win,gc,x1,y1,w,h);

			/* draw current box */
			x1 = (xb<xorig)?xb:xorig;
			y1 = (yb<yorig)?yb:yorig;
			x2 = (xb>xorig)?xb:xorig;
			y2 = (yb>yorig)?yb:yorig;
			w = x2-x1;
			h = y2-y1;
			XDrawRectangle(dpy,win,gc,x1,y1,w,h);

			/* remember current pointer position */
			xold = xb;
			yold = yb;
		}
	}

	/* erase rubber box */
	XDrawRectangle(dpy,win,gc,x1,y1,w,h);

	/* free graphics context */
	XFreeGC(dpy,gc);

	/* set output parameters */
	*x = x1;
	*y = y1;
	*width = w;
	*height = h;
}

