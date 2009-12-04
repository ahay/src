/* X WIGgle-trace plot of f(x1,x2) via Bitmap.

X Functionality:							
 Button 1	Zoom with rubberband box				
 Button 2	Show mouse (x1,x2) coordinates while pressed		
 q or Q key	Quit							
 s key		Save current mouse (x1,x2) location to file		
 p or P key	Plot current window with pswigb (only from disk files)	
 a or page up keys		enhance clipping by 10%			
 c or page down keys		reduce clipping by 10%			
 up,down,left,right keys	move zoom window by half width/height	
 i or +(keypad) 		zoom in by factor 2 			
 o or -(keypad) 		zoom out by factor 2 			
    l 				lock the zoom while moving the coursor	
    u 				unlock the zoom 			
 1,2,...,9	Zoom/Move factor of the window size			
									
 Notes:								
	Reaching the window limits while moving within changes the zoom	
	factor in this direction. The use of zoom locking(l) disables it
*/

/* XWIGB: $Revision: 1.46 $ ; $Date: 2007/05/04 18:08:27 $	*/

/*
 * AUTHOR:  Dave Hale, Colorado School of Mines, 08/09/90
 *
 * Endian stuff by: 
 *    Morten Wendell Pedersen, Aarhus University (visiting CSM, June 1995)
 *  & John Stockwell, Colorado School of Mines, 5 June 1995
 *
 * Stewart A. Levin, Mobil - Added ps print option
 * John Stockwell - Added optional sinc interpolation
 * Stewart A. Levin, Mobil - protect title, labels in pswigb call
 *
 * Brian J. Zook, SwRI - Added style=normal and wigclip flag
 *
 * Brian K. Macy, Phillips Petroleum, 11/27/98, added curve plotting option
 * Curve plotting notes:
 * MODIFIED:  P. Michaels, Boise State Univeristy  29 December 2000
 *            Added solid/grey color scheme for peaks/troughs
 * 
 * G.Klein, IFG Kiel University, 2002-09-29, added cursor scrolling and
 *            interactive change of zoom and clipping.
 *          IFM-GEOMAR Kiel, 2004-03-12, added zoom locking 
 *          IFM-GEOMAR Kiel, 2004-03-25, interactive plotting fixed 
 *
 * Modified by Sergey Fomel for including with Madagascar.
 */
/*
  Copyright � 2007, Colorado School of Mines,
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
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xatom.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>

#include <rsf.h>

#include "image.h"
#include "window.h"
#include "axesbox.h"
#include "rubberbox.h"
#include "rfwtva.h"
#include "rfwtvaint.h"

/* functions defined and used internally */
static void zoomBox (int x, int y, int w, int h, 
	int xb, int yb, int wb, int hb,
	float x1, float x2,
	float y1, float y2,
	float *x1b, float *x2b,
	float *y1b, float *y2b,
        int style);
static XImage *newBitmapCWP (Display *dpy, int width, int height,
	int n1, float d1, float f1, int n2, float *x2, float *z,
	float x1beg, float x1end, float x2beg, float x2end,
	float xcur, float clip, int wt, int va,
	float *p2begp, float *p2endp, int endian, int interp,
	int wigclip, int style);
void xMouseLoc(Display *dpy, Window win, XEvent event, int style, Bool show,
	int x, int y, int width, int height,
	float x1begb, float x1endb, float x2begb, float x2endb,
	float p2beg, float p2end);
void xMousePrint(XEvent event, int style, FILE *mpicksfp,
		 int x, int y, int width, int height,
	float x1begb, float x1endb, float x2begb, float x2endb,
	float p2beg, float p2end);
static XImage *RotImage90(Display *dpy, XImage *oldImage, int bitmap_pad);

int main (int argc, char *argv[])
{
    bool verbose, interp;
    int n1,n2,n1tic,n2tic,wt,va,
	i2,style=0,
	nz,iz,endian,
	xbox,ybox,wbox,hbox,
	xb,yb,wb,hb,
	x,y,width,height,
	imageOutOfDate,winwidth=-1,winheight=-1,
	showloc=0,wigclip;
    gridcode grid1,grid2;
    float labelsize,titlesize,perc,clip,xcur,bias,
	d1,f1,d2,f2,*z,*temp,*x2,
	x1beg,x1end,x2beg,x2end,
	x1min,x1max,x2min,x2max,
	d1num,f1num,d2num,f2num,
	x1begb,x1endb,x2begb,x2endb,p2beg=0,p2end=0;
    char *label1,*label2,*title,*windowtitle,
	*labelfont,*titlefont,
	*styles,*grid1s,*grid2s,
	*labelcolor,*titlecolor,
	*gridcolor,keybuf[256],*mpicks;
    FILE *mpicksfp;
    Display *dpy;
    Window win;
    XEvent event;
    KeySym keysym;
    XComposeStatus keystat;
    XImage *image=NULL;
    GC gci;
    int scr;
    unsigned long black,white;
    
    int lock=0;		/* lock/unlock zoom while scrolling */
    float mve;		/* distance for scrolling */
    float mvefac=8.;	/* window factor for scrolldistance 
			 * 2=half window size; 
			 * 8=one eighths of the window size */
    char  *msg="";		/* message on screen */
    
    sf_file in=NULL;

    /* initialize getpar */
    sf_init(argc,argv);
    in = sf_input("in");
    
    /* get parameters describing 1st dimension sampling */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input\n");
    if (!sf_histfloat(in,"d1",&d1)) d1 = 1.0;
    if (!sf_histfloat(in,"o1",&f1)) f1 = 0.0;
    
    x1min = (d1>0.0)?f1:f1+(n1-1)*d1;
    x1max = (d1<0.0)?f1:f1+(n1-1)*d1;
    
    /* get parameters describing 2nd dimension sampling */
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input\n");
    
    x2 = sf_floatalloc(n2);
    if (!sf_getfloats("x2",x2,n2)) {
	/* array of sampled values in 2nd dimension */
	if (!sf_histfloat(in,"d2",&d2)) d2 = 1.0;
	if (!sf_histfloat(in,"o2",&f2)) f2 = 0.0;
	for (i2=0; i2<n2; i2++)
	    x2[i2] = f2+i2*d2;
    }
    
    for (i2=1,x2min=x2max=x2[0]; i2<n2; i2++) {
	x2min = SF_MIN(x2min,x2[i2]);
	x2max = SF_MAX(x2max,x2[i2]);
    }
    
    if (NULL == (mpicks = sf_getstring("mpicks"))) mpicks = "/dev/tty";
    /* file to save mouse picks in */
    if (NULL == (mpicksfp = fopen(mpicks, "w"))) sf_error("Error opening %s:",mpicks);
    
    /* read binary data to be plotted */
    nz = n1*n2;
    z = sf_floatalloc(nz);
    sf_floatread(z,nz,in);
    
    /* if necessary, subtract bias */
    if (!sf_getfloat("bias",&bias)) bias=0.0;
    /* data value corresponding to location along axis 2 */
    
    if (bias!=0.0)
	for (iz=0; iz<nz; iz++)
	    z[iz] -= bias;
    
    /* if necessary, determine clip from percentile */
    if (!sf_getfloat("clip",&clip)) {
	/* data values < bias+clip and > bias-clip are clipped */
	
	if (!sf_getfloat("perc",&perc)) perc=100.0;
	/* percentile for determining clip */
	
	temp = sf_floatalloc(nz);
	for (iz=0; iz<nz; iz++)
	    temp[iz] = fabsf(z[iz]);
	iz = (nz*perc/100.0);
	if (iz<0) iz = 0;
	if (iz>nz-1) iz = nz-1;
	clip = sf_quantile(iz,nz,temp);
	free(temp);
    }
    
    if (!sf_getbool("verbose",&verbose)) verbose = true;
    /* y for info printed on stderr (n for no info) */
    if (verbose) sf_warning("clip=%g\n",clip);

    if (!sf_getint("wt",&wt)) wt = 1;
    /* =0 for no wiggle-trace; =1 for wiggle-trace */
    if (!sf_getint("va",&va)) va = 1;
    /* =0 for no variable-area; 
       =1 for variable-area fill;
       =2 for variable area, solid/grey fill */

    /* set wt=va for va with solid/grey coloring  */
    if (va>=2) {  wt=va; va=1; } 

    if (!sf_getfloat("xcur",&xcur)) xcur = 1.0;
    /* wiggle excursion in traces corresponding to clip */

    if (!sf_getint("wigclip",&wigclip)) wigclip = 0;
    /* If 0, the plot box is expanded to accommodate	
       the larger wiggles created by xcur>1. If this 
       flag is non-zero, the extra-large wiggles are	
       are clipped at the boundary of the plot box. */

    if (!sf_getint("xbox",&xbox)) xbox = 50;
    /*  x in pixels of upper left corner of window */
    if (!sf_getint("ybox",&ybox)) ybox = 50;
    /* y in pixels of upper left corner of window */
    if (!sf_getint("wbox",&wbox)) wbox = 550;
    /* width in pixels of window */
    if (!sf_getint("hbox",&hbox)) hbox = 700;
    /* height in pixels of window */

    if (!sf_getfloat("x1beg",&x1beg)) x1beg = x1min;
    /* value at which axis 1 begins */
    if (!sf_getfloat("x1end",&x1end)) x1end = x1max;
    /* value at which axis 1 ends */
    if (!sf_getfloat("d1num",&d1num)) d1num = 0.0;
    /* numbered tic interval on axis 1 (0.0 for automatic) */
    if (!sf_getfloat("f1num",&f1num)) f1num = x1min;
    /* first numbered tic on axis 1 (used if d1num not 0.0) */
    if (!sf_getint("n1tic",&n1tic)) n1tic = 1;
    /* number of tics per numbered tic on axis 1 */
    if (NULL == (grid1s = sf_getstring("grid1"))) grid1s="none";
    /* grid lines on axis 1 - none, dot, dash, or solid */
    if (0 == strcmp("dot",grid1s)) {
	grid1 = DOT;
    } else if (0 == strcmp("dash",grid1s)) {
	grid1 = DASH;
    } else if (0 == strcmp("solid",grid1s)) {
	grid1 = SOLID;
    } else {
	grid1 = NONE;
    }
    if (NULL == (label1 = sf_getstring("label1")) &&
	NULL == (label1 = sf_histstring(in,"label1"))) label1="";
    
    if (!sf_getfloat("x2beg",&x2beg)) x2beg = x2min;
    /* value at which axis 2 begins */
    if (!sf_getfloat("x2end",&x2end)) x2end = x2max;
    /* value at which axis 2 ends */
    if (!sf_getfloat("d2num",&d2num)) d2num = 0.0;
    /* numbered tic interval on axis 2 (0.0 for automatic) */
    if (!sf_getfloat("f2num",&f2num)) f2num = x2min;
    /* first numbered tic on axis 2 (used if d2num not 0.0) */
    if (!sf_getint("n2tic",&n2tic)) n2tic = 1;
    /* number of tics per numbered tic on axis 2 */
    if (NULL == (grid2s = sf_getstring("grid2"))) grid2s="none";
    /* grid lines on axis 2 - none, dot, dash, or solid */
    if (0 == strcmp("dot",grid2s)) {
	grid2 = DOT;
    } else if (0 == strcmp("dash",grid2s)) {
	grid2 = DASH;
    } else if (0 == strcmp("solid",grid2s)) {
	grid2 = SOLID;
    } else {
	grid2 = NONE;
    }
    if (NULL == (label2 = sf_getstring("label2")) &&
	NULL == (label2 = sf_histstring(in,"label2"))) label2="";

    if (NULL == (labelfont = sf_getstring("labelfont"))) labelfont="Erg14";
    /* font name for axes labels */
    if (!sf_getfloat("labelsize",&labelsize)) labelsize = 18.0; 
    if (NULL == (title = sf_getstring("title"))) title="";
    if (NULL == (titlefont = sf_getstring("titlefont"))) titlefont="Rom22";
    if (!sf_getfloat("titlesize",&titlesize)) titlesize = 24.0; 

    if (NULL == (styles = sf_getstring("style"))) styles="seismic";
    if (0 == strcmp("seismic",styles)) {
	style = SEISMIC;
    } else if (0 == strcmp("normal",styles)) {
	style = NORMAL;
    } else {
	sf_error("Unknown style='%s'\n",styles);
    }

    if (NULL == (titlecolor  = sf_getstring("titlecolor")))  titlecolor  = "red";
    if (NULL == (labelcolor  = sf_getstring("labelcolor")))  labelcolor  = "blue";
    if (NULL == (gridcolor   = sf_getstring("gridcolor")))    gridcolor  = "blue";
    if (NULL == (windowtitle = sf_getstring("windowtitle"))) windowtitle = "sfwigb";

    /* initialize zoom box parameters */
    x1begb = x1beg;	 x1endb = x1end;
    x2begb = x2beg;	 x2endb = x2end;

    /* connect to X server */
    if ((dpy=XOpenDisplay(NULL))==NULL)
	sf_error("Cannot connect to display %s!\n",XDisplayName(NULL));
    scr = DefaultScreen(dpy);
    black = BlackPixel(dpy,scr);
    white = WhitePixel(dpy,scr);

    if (!sf_getint("endian",&endian)){
	/* endian for display =0 little endian =1 big endian */
	if(BitmapBitOrder(dpy)==LSBFirst)
	    endian=0;
	else if(BitmapBitOrder(dpy)==MSBFirst)
	    endian=1;
	else if (sf_endian())
	    endian=0;
	else
	    endian=1;
    }
 
    if (!sf_getbool("interp",&interp)) interp = false;
    /* if y, use interpolation */
 
    /* create window */
    win = xNewWindow(dpy,xbox,ybox,wbox,hbox,
		     (int) black,(int) white,windowtitle);
		
    /* make GC for image */
    gci = XCreateGC(dpy,win,0,NULL);

    /* make sure foreground/background are black/white */
    XSetForeground(dpy,gci,black);
    XSetBackground(dpy,gci,white);

	/* set normal event mask */
	XSelectInput(dpy,win,
		StructureNotifyMask |
		ExposureMask |
		KeyPressMask |
		PointerMotionMask |
		ButtonPressMask |
		ButtonReleaseMask |
		Button1MotionMask |
		Button2MotionMask);
	
	/* map window */
	XMapWindow(dpy,win);
	
	/* clear the window */
	XClearWindow(dpy,win);
					
	/* determine good size for axes box */
	xSizeAxesBox(dpy,win,
		labelfont,titlefont,style,
		&x,&y,&width,&height);
	
	/* note that image is out of date */
	imageOutOfDate = 1;

	/* main event loop */
	while(imageOutOfDate|(~imageOutOfDate)/*True*/) {
		XNextEvent(dpy,&event);

		/* if window was resized */
		if (event.type==ConfigureNotify &&
			(event.xconfigure.width!=winwidth ||
			 event.xconfigure.height!=winheight)) {
			winwidth = event.xconfigure.width;
			winheight = event.xconfigure.height;
							
			/* determine good size for axes box */
			xSizeAxesBox(dpy,win,
				labelfont,titlefont,style,
				&x,&y,&width,&height);
			
			/* clear the window */
			XClearWindow(dpy,win);
			
			/* note that image is out of date */
			imageOutOfDate = 1;

		/* else if window exposed */
		} else if (event.type==Expose) {
			
			/* clear all expose events from queue */
			while (XCheckTypedEvent(dpy,Expose,&event));
			
			/* if necessary, make new image */
			if (imageOutOfDate) {
				if (image!=NULL) {
				/*	free1(image->data); */
					XDestroyImage(image);
				}
				image = newBitmapCWP(dpy,width,height,
					n1,d1,f1,n2,x2,z,
					x1begb,x1endb,x2begb,x2endb,
					xcur,clip,wt,va,
					&p2beg,&p2end,endian,interp,
					wigclip,style);
				imageOutOfDate = 0;
			}
	
			/* draw image (before axes so grid lines visible) */
			XPutImage(dpy,win,gci,image,0,0,x,y,
				image->width,image->height);

			/* draw axes on top of image */
			xDrawAxesBox(dpy,win,
				x,y,width,height,
				x1begb,x1endb,0.0,0.0,
				     d1num,f1num,n1tic,grid1,label1,
				x2begb,x2endb,p2beg,p2end,
				     d2num,f2num,n2tic,grid2,label2,
				labelfont,title,titlefont,
				labelcolor,titlecolor,gridcolor,
				style);

		/* else if key down */
		} else if (event.type==KeyPress) {

		    XLookupString(&(event.xkey),keybuf,0,&keysym,&keystat);
		    if (keysym==XK_s) {
			xMousePrint(event,style, mpicksfp,
				    x,y,width,height,
				    x1begb,x1endb,x2begb,x2endb,
				    p2beg, p2end);
		    } else if (keysym==XK_l ) {
			/* set lock */		  
			lock = 1 ;
			if (verbose) sf_warning("zoom lock set  %d\n",lock);

		    } else if (keysym==XK_u ) {
			/* unset lock */		  
			lock = 0 ;
			if (verbose) sf_warning("zoom lock released %d\n",lock);

		    } else if (keysym==XK_Shift_L ) { 
			/* if (verbose) 
			   fprintf(stderr,"Shift Left pressed \n");*/
		    } else if (keysym==XK_KP_1 || keysym==XK_1 ) { 
			mvefac=1.;
			sf_warning("Zoom/Move factor = 1");
		    } else if (keysym==XK_KP_2 || keysym==XK_2 ) { 
			mvefac=2.;
			sf_warning("Zoom/Move factor = 2");
		    } else if (keysym==XK_KP_3 || keysym==XK_3 ) { 
			mvefac=3.;
			if (verbose) 
			    sf_warning("Zoom/Move factor = 3");
		    } else if (keysym==XK_KP_4 || keysym==XK_4 ) { 
			mvefac=4.;
			if (verbose) 
			    sf_warning("Zoom/Move factor = 4");
		    } else if (keysym==XK_KP_5 || keysym==XK_5 ) { 
			mvefac=5.;
			if (verbose) 
			    sf_warning("Zoom/Move factor = 5");
		    } else if (keysym==XK_KP_6 || keysym==XK_6 ) { 
			mvefac=6.;
			if (verbose) 
			    sf_warning("Zoom/Move factor = 6");
		    } else if (keysym==XK_KP_7 || keysym==XK_7 ) { 
			mvefac=7.;
			if (verbose) 
			    sf_warning("Zoom/Move factor = 7");
		    } else if (keysym==XK_KP_8 || keysym==XK_8 ) { 
			mvefac=8.;
			if (verbose) 
			    sf_warning("Zoom/Move factor = 8");
		    } else if (keysym==XK_KP_9 || keysym==XK_9 ) { 
			mvefac=9.;
			if (verbose) 
			    sf_warning("Zoom/Move factor = 9");
		    } else if (keysym==XK_Left || keysym==XK_KP_Left ) {
			
			/* move zoom box to left by half window width */
			mve = (x2endb - x2begb)/mvefac ;
			x2begb = x2begb - mve ;
			x2endb = x2endb - mve ;
			msg="move "; 
			/* check for bounds of full window */
			if (x2begb < x2beg) {
			    if ( lock ) { x2begb = x2begb + mve ;
				x2endb = x2endb + mve ;
				msg="limit ";
				mve=0;
			    } else { x2begb = x2beg ;}
			}
			
			if (verbose) sf_warning("%s %g",msg,mve);
			
			/* clear area and force an expose event */
			XClearArea(dpy,win,0,0,0,0,True);
			/* note that image is out of date */
			imageOutOfDate = 1;
			
		    } else if (keysym==XK_Right || keysym==XK_KP_Right ) {
			/* move zoom box to right by half window width*/
			mve = (x2endb - x2begb)/mvefac ;
			x2begb = x2begb + mve ;
			x2endb = x2endb + mve ;
			msg="move "; 
			/* check for bounds of full window */
			if (x2endb > x2end) {
			    if ( lock ) { x2begb = x2begb - mve ;
				x2endb = x2endb - mve ;
				msg="limit ";
				mve=0;
			    } else { x2endb = x2end ;}
			}
			if (verbose) sf_warning("%s %g",msg,mve);
			
			/* clear area and force an expose event */
			XClearArea(dpy,win,0,0,0,0,True);
			/* note that image is out of date */
			imageOutOfDate = 1;
			
		    } else if (keysym==XK_Down || keysym==XK_KP_Down  ) {
			/* move zoom box down by half window height */
			mve = (x1endb - x1begb)/mvefac ;
			x1begb = x1begb + mve ;
			x1endb = x1endb + mve ;
			msg="move "; 
			/* check for bounds of full window */
			if (x1endb > x1end) {
			    if ( lock ) {  x1begb = x1begb - mve ;
				x1endb = x1endb - mve ;
				msg="limit ";
				mve=0;
			    } else { x1endb = x1end ;}
			}
			if (verbose) sf_warning("%s %g",msg,mve);
			
			/* clear area and force an expose event */
			XClearArea(dpy,win,0,0,0,0,True);
			/* note that image is out of date */
			imageOutOfDate = 1;
			
		    } else if (keysym==XK_Up || keysym==XK_KP_Up ) {
			/* move zoom box down by half window height */
			mve = (x1endb - x1begb)/mvefac ;
			x1begb = x1begb - mve ;
			x1endb = x1endb - mve ;
			msg="move "; 
			/* check for bounds of full window */
			if (x1begb < x1beg) {
			    if ( lock ) { x1begb = x1begb + mve ;
				x1endb = x1endb + mve ;
				msg="limit ";
				mve=0;
			    } else { x1begb = x1beg ;}
			}
			if (verbose) sf_warning("%s %g",msg,mve);
				
			/* clear area and force an expose event */
			XClearArea(dpy,win,0,0,0,0,True);
			
			/* note that image is out of date */
			imageOutOfDate = 1;
								
		    } else if (keysym==XK_o || keysym==XK_KP_Subtract ) {
			/* zoom out .... vertical*/
			mve = (x1endb - x1begb)/mvefac ;
			x1begb = x1begb - mve ;
			x1endb = x1endb + mve ;
			/* check for bounds of full window */
			if (x1begb < x1beg) x1begb = x1beg ;
			if (x1endb > x1end) x1endb = x1end ;
			/*   .... and horizontal */
			mve = (x2endb - x2begb)/mvefac ;
			x2begb = x2begb - mve ;
			x2endb = x2endb + mve ;
			/* check bounds of original image */
			if (x2begb < x2beg) x2begb = x2beg ;
			if (x2endb > x2end) x2endb = x2end ;
			
			/* clear area and force an expose event */
			XClearArea(dpy,win,0,0,0,0,True);
			
			/* note that image is out of date */
			imageOutOfDate = 1;
			
		    } else if (keysym==XK_i || keysym==XK_KP_Add ) {
			/* zoom in .... vertical*/
			mve = (x1endb - x1begb)/(2.*mvefac) ;
			x1begb = x1begb + mve ;
			x1endb = x1endb - mve ;
			/*   .... and horizontal */
			mve = (x2endb - x2begb)/(2.*mvefac) ;
			x2begb = x2begb + mve ;
			x2endb = x2endb - mve ;
			
			/* clear area and force an expose event */
			XClearArea(dpy,win,0,0,0,0,True);
			
			/* note that image is out of date */
			imageOutOfDate = 1;
			
		    } else if (keysym==XK_c || keysym==XK_Page_Down) {
		  	
			/* Change clip for image */
			clip += clip/10. ;
			/* if (verbose) warn("clip=%g\n",clip);*/
			if (verbose) fprintf(stderr,"clip=%g\n",clip);
			/* note that image is out of date */
			imageOutOfDate = 1;				
			
		    } else if (keysym==XK_a || keysym==XK_Page_Up) {
			
			/* Change clip for image */
			clip -= clip/10. ;
			/* if (verbose) warn("clip=%g\n",clip);*/
			if (verbose) fprintf(stderr,"clip=%g\n",clip);
			/* note that image is out of date */
			imageOutOfDate = 1;				
			
		    } else if (keysym==XK_q || keysym==XK_Q) {
			/* This is the exit from the event loop */
			break;
		    }

		    /* else if button down (1 == zoom, 2 == mouse tracking */
		} else if (event.type==ButtonPress) {
		    /* if 1st button: zoom */
		    if (event.xbutton.button==Button1) {
			
			/* track pointer and get new box */
			xRubberBox(dpy,win,event,&xb,&yb,&wb,&hb);
			
			/* if new box has tiny width or height */
			if (wb<4 || hb<4) {
			    
			    /* reset box to initial values */
			    x1begb = x1beg;
			    x1endb = x1end;
			    x2begb = x2beg;
			    x2endb = x2end;
			    
			    /* else, if new box has non-zero width */
				/* if new box has zero width or height */
				} else {
			
					/* calculate new box parameters */
					zoomBox(x,y,width,height,
						xb,yb,wb,hb,
						x2begb,x2endb,
						x1begb,x1endb,
						&x2begb,&x2endb,
						&x1begb,&x1endb,
                                                style);
				}

				/* clear area and force an expose event */
				XClearArea(dpy,win,0,0,0,0,True);
			
				/* note that image is out of date */
				imageOutOfDate = 1;
			
			/* else if 2nd button down: display mouse coords */
			} else if (event.xbutton.button==Button2) {

				showloc = 1;
				xMouseLoc(dpy,win,event,style,showloc,
					  x,y,width,height,x1begb,x1endb,
					  x2begb,x2endb,p2beg,p2end);

			} else {
				continue;
			}

		/* else if pointer has moved */
		} else if (event.type==MotionNotify) {
			
			/* if button2 down, show mouse location */
			if (showloc)
				xMouseLoc(dpy,win,event,style,True,
					x,y,width,height,x1begb,x1endb,
					x2begb,x2endb,p2beg,p2end);

		/* else if button2 released, stop tracking */
		} else if (event.type==ButtonRelease &&
			   event.xbutton.button==Button2) {
			showloc = 0;
		}

	} /* end of event loop */

	/* close connection to X server */
	XCloseDisplay(dpy);

	free(z);
    
	exit(0);
}
			
/* update parameters associated with zoom box */
static void zoomBox (int x, int y, int w, int h, 
	int xb, int yb, int wb, int hb,
	float x1, float x2,
	float y1, float y2,
	float *x1b, float *x2b,
	float *y1b, float *y2b,
        int style)
{
	/* if width and/or height of box are zero, just copy values */
	if (wb==0 || hb==0) {
		*x1b = x1; *x2b = x2;
		*y1b = y1; *y2b = y2;
		return;		
	} 
	
	/* clip box */
	if (xb<x) {
		wb -= x-xb;
		xb = x;
	}
	if (yb<y) {
		hb -= y-yb;
		yb = y;
	}
	if (xb+wb>x+w) wb = x-xb+w;
	if (yb+hb>y+h) hb = y-yb+h;	
	
	/* determine box limits */
        if (style == SEISMIC) {
                *x1b = x1+(xb-x)*(x2-x1)/w;
                *x2b = x1+(xb+wb-x)*(x2-x1)/w;
                *y1b = y1+(yb-y)*(y2-y1)/h;
                *y2b = y1+(yb+hb-y)*(y2-y1)/h;
        } else {
                *x2b = x2+(yb-y)*(x1-x2)/h;
                *x1b = x2+(yb+hb-y)*(x1-x2)/h;
                *y1b = y1+(xb-x)*(y2-y1)/w;
                *y2b = y1+(xb+wb-x)*(y2-y1)/w;
        }
}

/* return pointer to new image bitmap of rasterized wiggles */
static XImage *newBitmapCWP (Display *dpy, int width, int height,
	int n1, float d1, float f1, int n2, float *x2, float *z,
	float x1beg, float x1end, float x2beg, float x2end,
	float xcur, float clip, int wt, int va,
	float *p2begp, float *p2endp, int endian, int interp,
	int wigclip, int style)
{
	int widthpad,nbpr,i1beg,i1end,if1r,n1r,b1fz,b1lz,i2,i,n2in;
	float x2min,x2max,p2beg,p2end,bscale,boffset,bxcur,bx2;
	unsigned char *bits;
	int scr=DefaultScreen(dpy);
	XImage *image,*image2;
	float	x2margin,clip1,clip2;
	int	bx1max,bx2min,bx2max,b2f,b2l;
	int	width1,height1;
	int bitmap_pad=0;

	/* Kludge to fix problem with XCreateImage introduced in */
	/* Xorg 7.0 update for security */
	if (BitmapPad(dpy)>16) {
		bitmap_pad = 16;
	} else if (BitmapPad(dpy) < 16) {
		bitmap_pad = 8;
	}
	

	/* determine bitmap dimensions and allocate space for bitmap */
	width1 =  (style==SEISMIC) ? width : height;
	height1 = (style==SEISMIC) ? height : width;
	widthpad = (1+(width1-1)/bitmap_pad)*bitmap_pad;
	nbpr = widthpad-1;
	bits = sf_ucharalloc(nbpr*height1);
	for (i=0; i<nbpr*height1; ++i) bits[i] = 0;

	/* determine number of traces that fall within axis 2 bounds */
	x2min = SF_MIN(x2beg,x2end);
	x2max = SF_MAX(x2beg,x2end);
	for (i2=0,n2in=0; i2<n2; i2++)
		if (x2[i2]>=x2min && x2[i2]<=x2max) n2in++;

	/* determine pads for wiggle excursion along axis 2 */
	xcur = fabs(xcur);
	if (n2in>1) xcur *= (x2max-x2min)/(n2in-1);
	x2margin = (wigclip && n2in>1) ? (x2max-x2min)/(2*(n2in-1)) : xcur;
	p2beg = (x2end>=x2beg)?-x2margin:x2margin;
	p2end = -p2beg;

	bx2min = 0;
	bx2max = width1 - 1;
	bx1max = height1 - 1;

	/* determine scale and offset to map x2 units to bitmap units */
	bscale = bx2max/(x2end+p2end-x2beg-p2beg);
	boffset = -(x2beg+p2beg)*bscale;
	bxcur = xcur*bscale;

	/* adjust x1beg and x1end to fall on sampled values */
	i1beg = SF_NINT((x1beg-f1)/d1);
	i1beg = SF_MAX(0,SF_MIN(n1-1,i1beg));
	x1beg = f1+i1beg*d1;
	i1end = SF_NINT((x1end-f1)/d1);
	i1end = SF_MAX(0,SF_MIN(n1-1,i1end));
	x1end = f1+i1end*d1;

	/* determine first sample and number of samples to rasterize */
	if1r = SF_MIN(i1beg,i1end);
	n1r = SF_MAX(i1beg,i1end)-if1r+1;

	/* determine bits corresponding to first and last samples */
	b1fz = (x1end > x1beg) ? 0 : bx1max;
	b1lz = (x1end > x1beg) ? bx1max : 0;

	/* rasterize traces */
	for (i2=0; i2<n2; i2++,z+=n1) {

		/* skip traces not in bounds */
		if (x2[i2]<x2min || x2[i2]>x2max) continue;

		/* determine bitmap coordinate of trace */
		bx2 = boffset+x2[i2]*bscale;
		b2f = (int)(bx2-bxcur);
		b2l = (int)(bx2+bxcur);
		clip1 = -clip;
		clip2 = clip;
		if (b2f < bx2min) {
			clip1 *= ((bx2-bx2min) / bxcur);
			b2f = bx2min;
		}
		if (b2l > bx2max) {
			clip2 *= ((bx2max-bx2) / bxcur);
			b2l = bx2max;
		}

		/* rasterize one trace */
		if (interp==0) { /* don't use interpolation */
			rfwtva(n1r,&z[if1r],clip1,clip2,va?0:clip2,
				b2f,b2l,b1fz,b1lz,
				wt,nbpr,bits,endian);
		} else { /* use 8 point sinc interpolation */
			rfwtvaint(n1r,&z[if1r],clip1,clip2,va?0:clip2,
				b2f,b2l,b1fz,b1lz,
				wt,nbpr,bits,endian);
		}
		
	}
	
	/* return axis 2 pads */
	*p2begp = p2beg;  *p2endp = p2end;


	/* get pointer to image */
	image = XCreateImage(	(Display *) dpy,
				(Visual *) DefaultVisual(dpy,scr),
				(unsigned int) 1,
				(int) XYBitmap,
				(int) 0,
				(char *) bits,
				(unsigned int) widthpad,
				(unsigned int) height1,
				(int) bitmap_pad,
				(int) nbpr);

	if (style == NORMAL) {
		image2 = RotImage90(dpy,image,bitmap_pad);
		XDestroyImage(image);
		image = image2;
	}

	return image;
}	

void xMouseLoc(Display *dpy, Window win, XEvent event, int style, Bool show,
	int x, int y, int width, int height,
	float x1begb, float x1endb, float x2begb, float x2endb,
	float p2beg, float p2end)
{
	static XFontStruct *fs=NULL;
	static XCharStruct overall;
	static GC gc;
	int dummy,xoffset=5,yoffset=5;
	float x1,x2;
	char string[256];

	/* if first time, get font attributes and make gc */
	if (fs==NULL) {
		fs = XLoadQueryFont(dpy,"fixed");
		gc = XCreateGC(dpy,win,0,NULL);

		/* make sure foreground/background are black/white */
		XSetForeground(dpy,gc,BlackPixel(dpy,DefaultScreen(dpy)));
		XSetBackground(dpy,gc,WhitePixel(dpy,DefaultScreen(dpy)));

		XSetFont(dpy,gc,fs->fid);
		overall.width = 1;
		overall.ascent = 1;
		overall.descent = 1;
	}

	/* erase previous string */
	XClearArea(dpy,win,xoffset,yoffset,
		overall.width,overall.ascent+overall.descent,False);

	/* if not showing, then return */
	if (!show) return;

	/* convert mouse location to (x1,x2) coordinates */
	if (style==NORMAL) {
		x1 = x1begb+(x1endb-x1begb)*(event.xmotion.x-x)/width;
		x2 = p2end+x2endb+(p2beg+x2begb-x2endb-p2end)*
			(event.xmotion.y-y)/height;
	} else {
		x1 = x1begb+(x1endb-x1begb)*(event.xmotion.y-y)/height;
		x2 = p2beg+x2begb+(p2end+x2endb-x2begb-p2beg)*
			(event.xmotion.x-x)/width;
	}

	/* draw string indicating mouse location */
	sprintf(string,"(%0.6g,%0.6g)",x1,x2);
	XTextExtents(fs,string,(int)strlen(string),&dummy,&dummy,&dummy,&overall);
	XDrawString(dpy,win,gc,xoffset,yoffset+overall.ascent,
		string,(int)strlen(string));
}

void xMousePrint(XEvent event, int style, FILE *mpicksfp,
		 int x, int y, int width, int height,
		 float x1begb, float x1endb, float x2begb, float x2endb,
		 float p2beg, float p2end)
{
	float x1;
	int x2;

	/* convert mouse location to (x1,x2) coordinates */
	if (style==NORMAL) {
		x1 = x1begb+(x1endb-x1begb)*(event.xmotion.x-x)/width;
		x2 = SF_NINT(p2end+x2endb+(p2beg+x2begb-x2endb-p2end)*
			     (event.xmotion.y-y)/height);
	} else {
		x1 = x1begb+(x1endb-x1begb)*(event.xmotion.y-y)/height;
		x2 = SF_NINT(p2beg+x2begb+(p2end+x2endb-x2begb-p2beg)*
			     (event.xmotion.x-x)/width);
	}

	/* write string indicating mouse location */
	fprintf(mpicksfp, "%0.6g  %d\n", x1, x2);
}


/*********************************************************/
static XImage *RotImage90(Display *dpy, XImage *oldImage, int bitmap_pad)
{
	int	widthpad,x1,y1,x2,y2,i,nbpr;
	unsigned char	*bits;
	XImage	*image;
	int	width1 =		oldImage->width;
	int	width2 =		oldImage->height;
	int	height2 =	oldImage->width;
	int	scr =			DefaultScreen(dpy);

	widthpad = (1 + (width2-1)/bitmap_pad)*bitmap_pad;
	nbpr =  widthpad - 1;
	bits = sf_ucharalloc(nbpr*height2);
	image = XCreateImage(	(Display *) dpy,
				(Visual *) DefaultVisual(dpy,scr),
				(unsigned int) 1,
				(int) XYBitmap,
				(int) 0,
				(char *) bits,
				(unsigned int) widthpad,
				(unsigned int) height2,
				(int) bitmap_pad,
				(int) nbpr);

	for (i = 0; i < nbpr*height2; i++)	bits[i]=0;
	for (x2 = 0; x2 < width2; x2++) {
		y1 = x2;
		for (y2 = 0; y2 < height2; y2++) {
			x1 = width1 - 1 - y2;
			XPutPixel(image,x2,y2,XGetPixel(oldImage,x1,y1));
		}
	}
	return image;
}
