/*
  Copyright © 2007, Colorado School of Mines,
  All rights reserved.
  
  
  Redistribution and use in source and binary forms, with or 
  without modification, are permitted provided that the following 
  conditions are met:
  
  *  Redistributions of source code must retain the above copyright 
  notice, this list of conditions and the following disclaimer.
  *  Redistributions in binary form must reproduce the above 
  copyright notice, this list of conditions and the following 
  disclaimer in the documentation and/or other materials provided 
  with the distribution.
  *  Neither the name of the Colorado School of Mines nor the names of
  its contributors may be used to endorse or promote products 
  derived from this software without specific prior written permission.
  
  Warranty Disclaimer:
  THIS SOFTWARE IS PROVIDED BY THE COLORADO SCHOOL OF MINES AND CONTRIBUTORS 
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
  COLORADO SCHOOL OF MINES OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
  POSSIBILITY OF SUCH DAMAGE.
  
  
  Export Restriction Disclaimer:
  We believe that CWP/SU: Seismic Un*x is a low technology product that does
  not appear on the Department of Commerce CCL list of restricted exports.
  Accordingly, we believe that our product meets the qualifications of
  an ECCN (export control classification number) of EAR99 and we believe
  it fits the qualifications of NRR (no restrictions required), and
  is thus not subject to export restrictions of any variety.
  
  Approved Reference Format:
  In publications, please refer to SU as per the following example:
  Cohen, J. K. and Stockwell, Jr. J. W., (200_), CWP/SU: Seismic Un*x 
  Release No. __: an open source software  package for seismic 
  research and processing, 
  Center for Wave Phenomena, Colorado School of Mines.
  
  Articles about SU in peer-reviewed journals:
  Saeki, T., (1999), A guide to Seismic Un*x (SU)(2)---examples of data processing (part 1), data input and preparation of headers, Butsuri-Tansa (Geophysical Exploration), vol. 52, no. 5, 465-477.
  Stockwell, Jr. J. W. (1999), The CWP/SU: Seismic Un*x Package, Computers and Geosciences, May 1999.
  Stockwell, Jr. J. W. (1997), Free Software in Education: A case study of CWP/SU: Seismic Un*x, The Leading Edge, July 1997.
  Templeton, M. E., Gough, C.A., (1998), Web Seismic Un*x: Making seismic reflection processing more accessible, Computers and Geosciences.
  
  Acknowledgements:
  SU stands for CWP/SU:Seismic Un*x, a processing line developed at Colorado 
  School of Mines, partially based on Stanford Exploration Project (SEP) 
  software.
*/

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

