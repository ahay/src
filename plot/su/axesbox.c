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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <X11/Xlib.h>

#include "scaxis.h"

typedef enum {
    NONE,
    DOT,
    DASH,
    SOLID
} gridcode;
/*^*/

enum {
    NORMAL,
    SEISMIC
};
/*^*/

void
xDrawAxesBox (Display *dpy     /* display pointer */, 
	      Window win       /* window */,
	      int x            /* x coordinate of upper left corner of box */, 
	      int y            /* y coordinate of upper left corner of box */, 
	      int width        /* width of box */, 
	      int height       /* height of box */,
	      float x1beg      /* axis value at beginning of axis 1 */, 
	      float x1end      /* axis value at end of axis 1 */, 
	      float p1beg      /* pad value at beginning of axis 1 */, 
	      float p1end      /* pad value at end of axis 1 */,
	      float d1num      /* numbered tic increment for axis 1 (0.0 for automatic) */, 
	      float f1num      /* first numbered tic for axis 1 */, 
	      int n1tic        /* number of tics per numbered tic for axis 1 */, 
	      gridcode grid1   /* grid code for axis 1 */, 
	      char *label1     /* label for axis 1 */,
	      float x2beg      /* axis value at beginning of axis 2 */, 
	      float x2end      /* axis value at end of axis 2 */, 
	      float p2beg      /* pad value at beginning of axis 2 */, 
	      float p2end      /* pad value at end of axis 2 */,
	      float d2num      /* numbered tic increment for axis 2 (0.0 for automatic) */, 
	      float f2num      /* first numbered tic for axis 2 */, 
	      int n2tic        /* number of tics per numbered tic for axis 2 */, 
	      gridcode grid2   /* grid code for axis 2 */, 
	      char *label2     /* label for axis 2 */,
	      char *labelfont  /* name of font to use for axes labels */, 
	      char *title      /* axes box title */, 
	      char *titlefont  /* name of font to use for title */, 
	      char *axescolor  /* name of color to use for axes */, 
	      char *titlecolor /* name of color to use for title */, 
	      char *gridcolor  /* name of color to use for grid */,
	      int style        /* NORMAL (axis 1 on bottom, axis 2 on left)
				  SEISMIC (axis 1 on left, axis 2 on top) */)
/*< draw a labeled axes box 

xDrawAxesBox will determine the numbered tic incremenet and first
numbered tic automatically, if the specified increment is zero.

Pad values must be specified in the same units as the corresponding
axes values.  These pads are useful when the contents of the axes box
requires more space than implied by the axes values.  For example,
the first and last seismic wiggle traces plotted inside an axes box
will typically extend beyond the axes values corresponding to the
first and last traces.  However, all tics will lie within the limits
specified in the axes values (x1beg, x1end, x2beg, x2end) >*/
/******************************************************************************
Author:		Dave Hale, Colorado School of Mines, 01/27/90
*****************************************************************************/
{
	GC gca,gct,gcg;
	XGCValues *values=NULL;
	XColor scolor,ecolor;
	XFontStruct *fa,*ft;
	XWindowAttributes wa;
	Colormap cmap;
	int labelca,labelcd,labelch,labelcw,titleca,
		ntic,xa,ya,tw,ticsize,ticb,numb,labelb,lstr,grided,grid,
		n1num,n2num;
	float dnum,fnum,dtic,amin,amax,base,scale,anum,atic,azero;
	char str[256],dash[2],*label;

	/* create graphics contexts */
	gca = XCreateGC(dpy,win,0,values);
	gct = XCreateGC(dpy,win,0,values);
	gcg = XCreateGC(dpy,win,0,values);

	/* get and set fonts and determine character dimensions */
	fa = XLoadQueryFont(dpy,labelfont);
	if (fa==NULL) fa = XLoadQueryFont(dpy,"fixed");
	if (fa==NULL) {
		fprintf(stderr,"Cannot load/query labelfont=%s\n",labelfont);
		exit(-1);
	}
	XSetFont(dpy,gca,fa->fid);
	labelca = fa->max_bounds.ascent;
	labelcd = fa->max_bounds.descent;
	labelch = fa->max_bounds.ascent+fa->max_bounds.descent;
	labelcw = fa->max_bounds.lbearing+fa->max_bounds.rbearing;
	ft = XLoadQueryFont(dpy,titlefont);
	if (ft==NULL) ft = XLoadQueryFont(dpy,"fixed");
	if (ft==NULL) {
		fprintf(stderr,"Cannot load/query titlefont=%s\n",titlefont);
		exit(-1);
	}
	XSetFont(dpy,gct,ft->fid);
	titleca = ft->max_bounds.ascent;

	/* determine window's current colormap */
	XGetWindowAttributes(dpy,win,&wa);
	cmap = wa.colormap;

	/* get and set colors */
	if (XAllocNamedColor(dpy,cmap,axescolor,&scolor,&ecolor))
		XSetForeground(dpy,gca,ecolor.pixel);
	else
		XSetForeground(dpy,gca,1L);
	if (XAllocNamedColor(dpy,cmap,titlecolor,&scolor,&ecolor))
		XSetForeground(dpy,gct,ecolor.pixel);
	else
		XSetForeground(dpy,gct,1L);
	if (XAllocNamedColor(dpy,cmap,gridcolor,&scolor,&ecolor))
		XSetForeground(dpy,gcg,ecolor.pixel);
	else
		XSetForeground(dpy,gcg,1L);

	/* determine tic size */
	ticsize = labelcw;

	/* determine numbered tic intervals */
	if (d1num==0.0) {
		n1num = (style==NORMAL ? width : height)/(8*labelcw);
		scaxis(x1beg,x1end,&n1num,&d1num,&f1num);
	}
	if (d2num==0.0) {
		n2num = (style==NORMAL ? height : width)/(8*labelcw);
		scaxis(x2beg,x2end,&n2num,&d2num,&f2num);
	}

	/* draw horizontal axis */
	if (style==NORMAL) {
		amin = (x1beg<x1end)?x1beg:x1end;
		amax = (x1beg>x1end)?x1beg:x1end;
		dnum = d1num;  fnum = f1num;  ntic = n1tic;
		scale = width/(x1end+p1end-x1beg-p1beg);
		base = x-scale*(x1beg+p1beg);
		ya = y+height;
		ticb = ticsize;
		numb = ticb+labelca;
		labelb = numb+labelch;
		grid = grid1;
		label = label1;
	} else {
		amin = (x2beg<x2end)?x2beg:x2end;
		amax = (x2beg>x2end)?x2beg:x2end;
		dnum = d2num;  fnum = f2num;  ntic = n2tic;
		scale = width/(x2end+p2end-x2beg-p2beg);
		base = x-scale*(x2beg+p2beg);
		ya = y;
		ticb = -ticsize;
		numb = ticb-labelcd;
		labelb = numb-labelch;
		grid = grid2;
		label = label2;
	}
	if (grid==SOLID)
		grided = True;
	else if (grid==DASH) {
		grided = True;
		XSetLineAttributes(dpy,gcg,1L,LineOnOffDash,CapButt,JoinMiter);
		dash[0] = 8;  dash[1] = 4;
		XSetDashes(dpy,gcg,0,dash,2);
	} else if (grid==DOT) {
		grided = True;
		XSetLineAttributes(dpy,gcg,1L,LineOnOffDash,CapButt,JoinMiter);
		dash[0] = 1;  dash[1] = 4;
		XSetDashes(dpy,gcg,0,dash,2);
	} else
		grided = False;
	azero = 0.0001*(amax-amin);
	for (anum=fnum; anum<=amax; anum+=dnum) {
		if (anum<amin) continue;
		xa = base+scale*anum;
		if (grided) XDrawLine(dpy,win,gcg,xa,y,xa,y+height);
		XDrawLine(dpy,win,gca,xa,ya,xa,ya+ticb);
		if (anum>-azero && anum<azero)
			sprintf(str,"%1.5g",0.0);
		else
			sprintf(str,"%1.5g",anum);
		lstr = (int) strlen(str);
		tw = XTextWidth(fa,str,lstr);
		XDrawString(dpy,win,gca,xa-tw/2,ya+numb,str,lstr);
	}
	dtic = dnum/ntic;
	for (atic=fnum-ntic*dtic-dtic; atic<=amax; atic+=dtic) {
		if (atic<amin) continue;
		xa = base+scale*atic;
		XDrawLine(dpy,win,gca,xa,ya,xa,ya+ticb/2);
	}
	lstr = (int) strlen(label);
	tw = XTextWidth(fa,label,lstr);
	XDrawString(dpy,win,gca,x+width-tw,ya+labelb,label,lstr);

	/* draw vertical axis */
	if (style==NORMAL) {
		amin = (x2beg<x2end)?x2beg:x2end;
		amax = (x2beg>x2end)?x2beg:x2end;
		dnum = d2num;  fnum = f2num;  ntic = n2tic;
		scale = -height/(x2end+p2end-x2beg-p2beg);
		base = y+height-scale*(x2beg+p2beg);
		grid = grid2;
		label = label2;
	} else {
		amin = (x1beg<x1end)?x1beg:x1end;
		amax = (x1beg>x1end)?x1beg:x1end;
		dnum = d1num;  fnum = f1num;  ntic = n1tic;
		scale = height/(x1end+p1end-x1beg-p1beg);
		base = y-scale*(x1beg+p1beg);
		grid = grid1;
		label = label1;
	}
	xa = x;
	ticb = -ticsize;
	numb = ticb-ticsize/4;
	if (grid==SOLID)
		grided = True;
	else if (grid==DASH) {
		grided = True;
		XSetLineAttributes(dpy,gcg,1L,LineOnOffDash,CapButt,JoinMiter);
		dash[0] = 8;  dash[1] = 4;
		XSetDashes(dpy,gcg,0,dash,2);
	} else if (grid==DOT) {
		grided = True;
		XSetLineAttributes(dpy,gcg,1L,LineOnOffDash,CapButt,JoinMiter);
		dash[0] = 1;  dash[1] = 4;
		XSetDashes(dpy,gcg,0,dash,2);
	} else
		grided = False;
	azero = 0.0001*(amax-amin);
	for (anum=fnum; anum<=amax; anum+=dnum) {
		if (anum<amin) continue;
		ya = base+scale*anum;
		if (grided) XDrawLine(dpy,win,gcg,x,ya,x+width,ya);
		XDrawLine(dpy,win,gca,xa,ya,xa+ticb,ya);
		if (anum>-azero && anum<azero)
			sprintf(str,"%1.5g",0.0);
		else
			sprintf(str,"%1.5g",anum);
		lstr = (int) strlen(str);
		tw = XTextWidth(fa,str,lstr);
		XDrawString(dpy,win,gca,xa+numb-tw,ya+labelca/4,str,lstr);
	}
	dtic = dnum/ntic;
	for (atic=fnum-ntic*dtic-dtic; atic<=amax; atic+=dtic) {
		if (atic<amin) continue;
		ya = base+scale*atic;
		XDrawLine(dpy,win,gca,xa,ya,xa+ticb/2,ya);
	}
	lstr = (int) strlen(label);
	if (style==NORMAL)
		XDrawString(dpy,win,gca,
			x+ticb-9*labelcw,
			y+labelca/4-labelch,label,lstr);
	else
		XDrawString(dpy,win,gca,
			x+ticb-9*labelcw,
			y+height+labelca/4+labelch,label,lstr);
	
	/* draw title */
	lstr = (int) strlen(title);
	tw = XTextWidth(ft,title,lstr);
	if (style==NORMAL)
		XDrawString(dpy,win,gct,
			x+width/2-tw/2,
			y+labelca/4-labelch-labelch,title,lstr);
	else
		XDrawString(dpy,win,gct,
			x+width/2-tw/2,
			y+height+labelca/4+labelch+titleca,title,lstr);

	/* draw axes box */
	XDrawRectangle(dpy,win,gca,x,y,width,height);

	/* free resources before returning */
	XFreeGC(dpy,gca);
	XFreeGC(dpy,gct);
	XFreeGC(dpy,gcg);
	XFreeFont(dpy,fa);
	XFreeFont(dpy,ft);
}

void
xSizeAxesBox (Display *dpy    /* display pointer */, 
	      Window win      /* window */, 
	      char *labelfont /* name of font to use for axes labels */, 
	      char *titlefont /* name of font to use for title */, 
	      int style       /* NORMAL (axis 1 on bottom, axis 2 on left)
				 SEISMIC (axis 1 on left, axis 2 on top) */,
	      int *x          /* x coordinate of upper left corner of box */, 
	      int *y          /* y coordinate of upper left corner of box */, 
	      int *width      /* width of box */, 
	      int *height     /* height of box */)
/*< determine optimal origin and size for a labeled axes box 

xSizeAxesBox is intended to be used prior to xDrawAxesBox.

An "optimal" axes box is one that more or less fills the window, 
with little wasted space around the edges of the window. >*/
/******************************************************************************
Author:		Dave Hale, Colorado School of Mines, 01/27/90
*****************************************************************************/
{
	XFontStruct *fa,*ft;
	XWindowAttributes attr;
	int labelch,labelcw,titlech,bl,bt,br,bb;

	/* get fonts and determine character dimensions */
	fa = XLoadQueryFont(dpy,labelfont);
	if (fa==NULL) fa = XLoadQueryFont(dpy,"fixed");
	if (fa==NULL) {
		fprintf(stderr,"Cannot load/query labelfont=%s\n",labelfont);
		exit(-1);
	}
	labelch = fa->max_bounds.ascent+fa->max_bounds.descent;
	labelcw = fa->max_bounds.lbearing+fa->max_bounds.rbearing;
	ft = XLoadQueryFont(dpy,titlefont);
	if (ft==NULL) ft = XLoadQueryFont(dpy,"fixed");
	if (ft==NULL) {
		fprintf(stderr,"Cannot load/query titlefont=%s\n",titlefont);
		exit(-1);
	}
	titlech = ft->max_bounds.ascent+ft->max_bounds.descent;

	/* determine axes box origin and size */
	XGetWindowAttributes(dpy,win,&attr);
	bl = 10*labelcw;
	br = attr.width-5*labelcw;
	while (br<=bl) {
		br += labelcw;
		bl -= labelcw;
	}
	if (bl<0) bl = 0;
	if (br>attr.width) br = attr.width;
	if (style==NORMAL) {
		bt = labelch+labelch/2+titlech;
		bb = attr.height-3*labelch;
	} else {
		bt = 3*labelch;
		bb = attr.height-labelch-labelch/2-titlech;
	}
	while (bb<=bt) {
		bb += labelch;
		bt -= labelch;
	}
	if (bt<0) bt = 0;
	if (bb>attr.height) bb = attr.height;
	
	*x = bl;
	*y = bt;
	*width = br-bl;
	*height = bb-bt;

	XFreeFont(dpy,fa);
	XFreeFont(dpy,ft);
}
