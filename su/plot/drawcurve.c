/* Copyright (c) Colorado School of Mines, 2007.*/
/* All rights reserved.                       */

/* DRAWCURVE: $Revision: 1.2 $ ; $Date: 1999/01/14 22:50:24 $	*/

/*********************** self documentation **********************/
/*****************************************************************************
DRAWCURVE - Functions to draw a curve from a set of points

xDrawCurve	draw a curve from a set of points
*****************************************************************************
Function Prototypes:
void xDrawCurve(Display *dpy, Window win,
		int x, int y, int width, int height,
		float x1beg, float x1end, float p1beg, float p1end,
		float x2beg, float x2end, float p2beg, float p2end,
		float *x1curve, float *x2curve, int ncurve,
		char *curvecolor, int style);
*****************************************************************************
xDrawCurve:
Input:
dpy		display pointer
win		window
x		x coordinate of upper left corner of box
y		y coordinate of upper left corner of box
width		width of box
height		height of box
x1beg		axis value at beginning of axis 1
x1end		axis value at end of axis 1
p1beg		pad value at beginning of axis 1
p1end		pad value at end of axis 1
x2beg		axis value at beginning of axis 2
x2end		axis value at end of axis 2
p2beg		pad value at beginning of axis 2
p2end		pad value at end of axis 2
x1curve		vector of x1 coordinates for points along curve
x2curve		vector of x2 coordinates for points along curve
ncurve		number of points along curve
curvecolor	name of color to use for axes
int style	NORMAL (axis 1 on bottom, axis 2 on left)
		SEISMIC (axis 1 on left, axis 2 on top)

******************************************************************************
Author:		Brian Macy, Phillips Petroleum Co., 11/14/98
		(Adapted after Dave Hale's xDrawAxesBox routine)
*****************************************************************************/
/**************** end self doc ********************************/

#include <X11/Xlib.h>

#include <rsf.h>

#include "axesbox.h"

void xDrawCurve(Display *dpy, Window win,
		int x, int y, int width, int height,
		float x1beg, float x1end, float p1beg, float p1end,
		float x2beg, float x2end, float p2beg, float p2end,
		float *x1curve, float *x2curve, int ncurve,
		char *curvecolor, int style)
/*< draw a curve from a set of points >*/
/******************************************************************************
Input:
dpy		display pointer
win		window
x		x coordinate of upper left corner of box
y		y coordinate of upper left corner of box
width		width of box
height		height of box
x1beg		axis value at beginning of axis 1
x1end		axis value at end of axis 1
p1beg		pad value at beginning of axis 1
p1end		pad value at end of axis 1
x2beg		axis value at beginning of axis 2
x2end		axis value at end of axis 2
p2beg		pad value at beginning of axis 2
p2end		pad value at end of axis 2
x1curve		vector of x1 coordinates for points along curve
x2curve		vector of x2 coordinates for points along curve
ncurve		number of points along curve
curvecolor	name of color to use for axes
int style	NORMAL (axis 1 on bottom, axis 2 on left)
		SEISMIC (axis 1 on left, axis 2 on top)
******************************************************************************
Author:		Brian Macy, Phillips Petroleum Co., 11/14/98
		(Adapted after Dave Hale's xDrawAxesBox routine)
*****************************************************************************/
{
	GC gcc;
	XGCValues *values=NULL;
	XColor scolor,ecolor;
	XWindowAttributes wa;
	Colormap cmap;
	float xbase,ybase,xscale,yscale;
	/* float xamin,xamax,yamin,yamax; */
	float *xcurve,*ycurve;
	XPoint *lpoints;   /* points for drawing line */
	XRectangle rectclip;
	int i;

	/* allocate memory for lpoints */
	lpoints=(XPoint *) sf_alloc(ncurve,sizeof(XPoint));

	/* create graphics contexts */
	gcc = XCreateGC(dpy,win,0,values);

	/* determine window's current colormap */
	XGetWindowAttributes(dpy,win,&wa);
	cmap = wa.colormap;

	/* get and set colors */
	if (XAllocNamedColor(dpy,cmap,curvecolor,&scolor,&ecolor))
		XSetForeground(dpy,gcc,ecolor.pixel);
	else
		XSetForeground(dpy,gcc,1L);

	if (style==NORMAL) {
	    /* xamin = (x1beg<x1end)?x1beg:x1end;
	       xamax = (x1beg>x1end)?x1beg:x1end; */
		xscale = width/(x1end+p1end-x1beg-p1beg);
		xbase = x-xscale*(x1beg+p1beg);
/*		yamin = (x2beg<x2end)?x2beg:x2end;
		yamax = (x2beg>x2end)?x2beg:x2end; */
		yscale = -height/(x2end+p2end-x2beg-p2beg);
		ybase = y+height-yscale*(x2beg+p2beg);
		xcurve=x1curve;
		ycurve=x2curve;
	} else {
/*		xamin = (x2beg<x2end)?x2beg:x2end;
		xamax = (x2beg>x2end)?x2beg:x2end; */
		xscale = width/(x2end+p2end-x2beg-p2beg);
		xbase = x-xscale*(x2beg+p2beg);
/*		yamin = (x1beg<x1end)?x1beg:x1end;
		yamax = (x1beg>x1end)?x1beg:x1end; */
		yscale = height/(x1end+p1end-x1beg-p1beg);
		ybase = y-yscale*(x1beg+p1beg);
		ycurve=x1curve;
		xcurve=x2curve;
	}

	/* Set up clip rectangle (only allows drawing inside rect.) */
	rectclip.x = (short) x;
	rectclip.y = (short) y;
	rectclip.width  = (unsigned short) width;
	rectclip.height = (unsigned short) height;
	XSetClipRectangles(dpy,gcc,0,0,&rectclip,1,Unsorted);

	/*
	 * Draw a curve from the input data xcurve,ycurve.
	 */
	for (i=0; i<ncurve; ++i) {
		lpoints[i].x=(short) (xbase+xscale*xcurve[i]);
		lpoints[i].y=(short) (ybase+yscale*ycurve[i]);
	}
	if (ncurve==1)
		XDrawPoints(dpy,win,gcc,lpoints,ncurve,CoordModeOrigin);
	else if (ncurve > 1)
		XDrawLines(dpy,win,gcc,lpoints,ncurve,CoordModeOrigin);

	/* free resources before returning */
	XFreeGC(dpy,gcc);
	free(lpoints);
}
