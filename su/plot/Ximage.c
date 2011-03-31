/* X IMAGE plot of a uniformly-sampled function f(x1,x2).
   X Functionality:							
   Button 1	Zoom with rubberband box				
   Button 2	Show mouse (x1,x2) coordinates while pressed		
   q or Q key	Quit							
   s key		Save current mouse (x1,x2) location to file		
   a or page up keys		enhance clipping by 10%			
   c or page down keys		reduce clipping by 10%			
   up,down,left,right keys	move zoom window by half width/height	
   i or +(keypad) 		zoom in by factor 2 			
   o or -(keypad) 		zoom out by factor 2 			
									
   ... change colormap interactively					
   r	     install next RGB - colormap				
   R	     install previous RGB - colormap				
   h	     install next HSV - colormap				
   H	     install previous HSV - colormap				
   H	     install previous HSV - colormap				
   (Move mouse cursor out and back into window for r,R,h,H to take effect)
*/
/* XIMAGE: $Revision: 1.43 $ ; $Date: 2005/10/04 17:28:06 $	*/

/*
 * AUTHOR:  Dave Hale, Colorado School of Mines, 08/09/90
 *
 * Stewart A. Levin, Mobil - Added ps print option
 *
 * Brian Zook, Southwest Research Institute, 6/27/96, added blank option
 *
 * Toralf Foerster, Baltic Sea Research Institute, 9/15/96, new colormaps
 *
 * Berend Scheffers, Delft, colorbar (legend)
 *
 * Brian K. Macy, Phillips Petroleum, 11/27/98, added curve plotting option
 * 
 * G.Klein, GEOMAR Kiel, 2004-03-12, added cursor scrolling and
 *                                   interactive change of zoom and clipping.
 * 
 * Zhaobo Meng, ConocoPhillips, 12/02/04, added amplitude display
 * 
 * Garry Perratt, Geocon, 08/04/05, modified perc handling to center colorbar if balance==1.
 */

/* Copyright (c) Colorado School of Mines, 2007.
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
#include <math.h> /* for roundf */

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/keysym.h>

#include <rsf.h>

#include "axesbox.h"
#include "rubberbox.h"
#include "colormap.h"
#include "image.h"
#include "window.h"
#include "legendbox.h"
#include "drawcurve.h"
#include "intl2b.h"

#define   STREQ(s,t) (strcmp(s,t) == 0)
#define NINT(x) ((int) roundf(x))

/* functions defined and used internally */
/* ZM: interpolate the amplitude from dataset */
static float getamp(float *zz,float f1,float d1,int n1,
		    float f2,float d2,int n2,float x1,float x2,int verbose);
static void zoomBox (int x, int y, int w, int h, 
		     int xb, int yb, int wb, int hb,
		     int nx, int ix, float x1, float x2,
		     int ny, int iy, float y1, float y2,
		     int *nxb, int *ixb, float *x1b, float *x2b,
		     int *nyb, int *iby, float *y1b, float *y2b);
static unsigned char *newInterpBytes (int n1in, int n2in, unsigned char *bin,
				      int n1out, int n2out, int newInterpBytes);
void xMouseLoc(Display *dpy, Window win, XEvent event, int style, Bool show,
	       int x, int y, int width, int height,
	       float x1begb, float x1endb, float x2begb, float x2endb,
	       float *z, float f1, float d1, int n1,float f2, float d2, 
	       int n2, int verbose);
void xMousePrint(XEvent event, int style, FILE *mpicksfp,
		 int x, int y, int width, int height,
		 float x1begb, float x1endb, float x2begb, float x2endb);

int
main (int argc,char **argv)
{
    int n1,n2,n1tic,n2tic,
	i1,i2,style,
	n1c,n2c,i1beg,i1end,i2beg,i2end,i1c,i2c,
	nz,iz,i1step,i2step,
	xbox,ybox,wbox,hbox,
	xb,yb,wb,hb,
	x,y,width,height,
	i,j,nx,ny,nxb,nyb,ixb,iyb,
	imageOutOfDate,winwidth=-1,winheight=-1,
	showloc=0,
	lwidth,lheight,lx,ly; /* BEREND */
    gridcode grid1, grid2;
    bool balance, verbose, blockinterp, legend;

    int base;
    float fact;
    unsigned char* ptr;

    float labelsize,titlesize,perc,clip,bperc,wperc,bclip,wclip,
	d1,f1,d2,f2,*z,*temp,zscale,zoffset,zi,dx,dy,
	x1beg,x1end,x2beg,x2end,
	x1min,x1max,x2min,x2max,
	d1num,f1num,d2num,f2num,
	x1begb,x1endb,x2begb,x2endb,blank; /* dx,dy added GK */

    unsigned char *cz,*czp,*czb,*czbp,*czbi=NULL;
    char *label1,*label2,*title,*windowtitle,
	*units, *legendfont,
	*labelfont,*titlefont,
	*styles,*grid1s,*grid2s,
	*labelcolor,*titlecolor,
	*gridcolor,*cmap,keybuf[256],*mpicks;
    FILE *mpicksfp;
    Display *dpy;
    Window win;
    XEvent event;
    KeySym keysym;
    XComposeStatus keystat;
    XImage *image=NULL;
    XImage *image_legend=NULL; /* BEREND */
    unsigned char *data_legend; /* BEREND */
    GC gci;
    int scr;
    unsigned long black,white,pmin,pmax;

    float **x1curve,**x2curve;
    int curve,*npair,ncurvecolor;
    char **curvefile,**curvecolor=NULL;
    FILE *curvefp;
	
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
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) d1 = 1.0;
    if (!sf_histfloat(in,"o1",&f1)) f1 = 0.0;

    x1min = (d1>0.0)?f1:f1+(n1-1)*d1;
    x1max = (d1<0.0)?f1:f1+(n1-1)*d1;

    /* get parameters describing 2nd dimension sampling */
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&d2)) d2 = 1.0;
    if (!sf_histfloat(in,"o2",&f2)) f2 = 0.0;

    x2min = (d2>0.0)?f2:f2+(n2-1)*d2;
    x2max = (d2<0.0)?f2:f2+(n2-1)*d2;

    if (NULL == (mpicks = sf_getstring("mpicks"))) mpicks = "/dev/tty";
    /* file to save mouse picks */

    mpicksfp = fopen(mpicks, "w");
    if (NULL == mpicksfp) sf_error("Cannot open \"%s\" for writing:",mpicks);

    /* set up curve plotting */

    if (!sf_getint("ncurve",&curve)) curve=0;
    /* number of curves to draw */

    if (curve > 0) {
	curvefile = (char**) sf_alloc(curve,sizeof(char*));
	if (!sf_getstrings("curve",curvefile,curve)) sf_error("Need curve=");
	/* file(s) containing points to draw curve(s) */

	x1curve = (float**) sf_alloc(curve,sizeof(float*));
	x2curve = (float**) sf_alloc(curve,sizeof(float*));

	npair= sf_intalloc(curve);
	if (sf_getints("npair",npair,curve)) sf_error("Need npair=");
	/* number(s) of pairs in each curve file */	
    } else {
	npair=(int *)NULL;
	curvefile=(char **)NULL;
	x1curve=(float **)NULL;
	x2curve=(float **)NULL;
    }

    if (!sf_getint("ncurvecolor",&ncurvecolor)) ncurvecolor=0;
    /* number of curve colors */

    if (ncurvecolor < curve) {
	curvecolor=(char**) sf_alloc(curve,sizeof(char*));
	if (!sf_getstrings("curvecolor",curvecolor,ncurvecolor)) {
	    curvecolor[0] = strdup("blue\0");
	    ncurvecolor=1;
	}
	for (i=ncurvecolor; i<curve; i++)
	    curvecolor[i] = strdup(curvecolor[ncurvecolor-1]);
    } else if( ncurvecolor ) {
	curvecolor=(char**) sf_alloc(ncurvecolor,sizeof(char*));
	if (!sf_getstrings("curvecolor",curvecolor,ncurvecolor)) sf_error("Need curvecolor=");
	/* color(s) for curve(s) */
    }

    for (j=0; j<curve; j++) {
	curvefp=fopen(curvefile[j],"r");
	if (NULL == curvefp) sf_error("Cannot open \"%s\" for reading:",curvefile[j]);
	x1curve[j]=sf_floatalloc(npair[j]);
	x2curve[j]=sf_floatalloc(npair[j]);
	for (i=0; i<npair[j]; i++) {
	    if (EOF == fscanf(curvefp,"%f",&x1curve[j][i])) sf_error("%s: scan error:",__FILE__);
	    if (EOF == fscanf(curvefp,"%f",&x2curve[j][i])) sf_error("%s: scan error:",__FILE__);
	}
	fclose(curvefp);
    }

    if (!sf_histfloat(in,"d2",&d2)) d2 = 1.0;
    if (!sf_histfloat(in,"o2",&f2)) f2 = 0.0;

    x2min = (d2>0.0)?f2:f2+(n2-1)*d2;
    x2max = (d2<0.0)?f2:f2+(n2-1)*d2;

    /* read binary data to be plotted */
    nz = n1*n2;
    z = sf_floatalloc(nz);
    sf_floatread(z,nz,in);

    /* if necessary, determine clips from percentiles */
    if (sf_getfloat("clip",&clip)) {
	bclip = clip;
	wclip = -clip;
    }
    if ((!sf_getfloat("bclip",&bclip) || 
	 !sf_getfloat("wclip",&wclip)) &&
	!sf_getfloat("clip",&clip)) {
	if (!sf_getfloat("perc",&perc)) perc = 100.0;
	/* percentile used to determine clip */
	if (!sf_getbool("balance",&balance)) balance=false; 
	/* if false, determine bclip & wclip individually; 
	   else set them to the same abs value */
	temp = sf_floatalloc(nz);
	/* Modded by GCP to balance bclip & wclip */

	if (!balance)
	    for (iz=0; iz<nz; ++iz) {
		temp[iz] = z[iz];
	    } else { 
	    for (iz=0; iz<nz; ++iz) temp[iz] = abs(z[iz]);
	    perc=100.0;
	}

	/* End of modded code */
	if (!sf_getfloat("bclip",&bclip)) {
	    if (!sf_getfloat("bperc",&bperc)) bperc = perc;
	    /* percentile for determining black clip value */
	    iz = (nz*bperc/100.0);
	    if (iz<0) iz = 0;
	    if (iz>nz-1) iz = nz-1;
	    bclip = sf_quantile(iz,nz,temp);
	}
	if (!sf_getfloat("wclip",&wclip)) {
	    if (!sf_getfloat("wperc",&wperc)) wperc = 100.0-perc;
	    /* percentile for determining white clip value */
	    iz = (nz*wperc/100.0);
	    if (iz<0) iz = 0;
	    if (iz>nz-1) iz = nz-1;
	    /* Modded by GCP to balance bclip & wclip */
	    if (balance==0) wclip = sf_quantile(iz,nz,temp);
	    else wclip = -1*bclip;
	    /* End of modded code */
	}
	free(temp);
    }
    if (!sf_getbool("verbose",&verbose)) verbose=true;
    /* verbose mode */
    if (verbose) sf_warning("bclip=%g wclip=%g",bclip,wclip);

    if (NULL == (cmap = sf_getstring("cmap"))) {
	/* colormap specification (hsv# or rgb#) */
	cmap = sf_charalloc(5);
	strcpy(cmap,"gray");
    }

    if (!sf_getbool("blockinterp", &blockinterp)) blockinterp=false;
    /* whether to use block interpolation */

    /* get legend specs BEREND */

    if (!sf_getbool("legend", &legend)) legend = false; 
    /* if display the color scale */

    if (NULL == (units = sf_getstring("units"))) units=""; 
    /* unit label for legend */
    if (NULL == (legendfont = sf_getstring("legendfont"))) legendfont="times_roman10"; 
    /* font name for legend */

    if (!sf_getfloat("blank",&blank)) blank = 0.0f;
    /* portion of the lower range to blank out */

    /* get axes parameters */
    if (!sf_getint("xbox",&xbox)) xbox = 50; /* x in pixels of upper left corner of window */
    if (!sf_getint("ybox",&ybox)) ybox = 50; /* y in pixels of upper left corner of window */
    if (!sf_getint("wbox",&wbox)) wbox = 550; /* width in pixels of window */
    if (!sf_getint("hbox",&hbox)) hbox = 700; /* height in pixels of window */

    /* legend dimensions */
    if (!sf_getint("lwidth",&lwidth))	lwidth = 16; /* colorscale (legend) width in pixels */
    if (!sf_getint("lheight",&lheight))	lheight = hbox/3; /* colorscale (legend) height in pixels */
    if (!sf_getint("lx",&lx))	lx = 3; /* colorscale (legend) x-position in pixels */
    if (!sf_getint("ly",&ly))	ly = (hbox-lheight)/3; /* colorscale (legend) y-position in pixels */

    if (!sf_getfloat("x1beg",&x1beg)) x1beg = x1min; /* value at which axis 1 begins */
    if (!sf_getfloat("x1end",&x1end)) x1end = x1max; /* value at which axis 1 ends */
    if (!sf_getfloat("d1num",&d1num)) d1num = 0.0; /* numbered tic interval on axis 1 (0.0 for automatic) */
    if (!sf_getfloat("f1num",&f1num)) f1num = x1min; /* first numbered tic on axis 1 (used if d1num not 0.0) */
    if (!sf_getint("n1tic",&n1tic)) n1tic = 1; /* number of tics per numbered tic on axis 1 */

    if (NULL == (grid1s = sf_getstring("grid1"))) grid1s="none"; /* grid lines on axis 1 (none, dot, dash, or solid) */ 
    if (STREQ("dot",grid1s)) grid1 = DOT;
    else if (STREQ("dash",grid1s)) grid1 = DASH;
    else if (STREQ("solid",grid1s)) grid1 = SOLID;
    else grid1 = NONE;
    if (NULL == (label1 = sf_getstring("label1"))) label1=""; /* label on axis 1 */

    if (!sf_getfloat("x2beg",&x2beg)) x2beg = x2min; /* value at which axis 2 begins */
    if (!sf_getfloat("x2end",&x2end)) x2end = x2max; /* value at which axis 2 ends */
    if (!sf_getfloat("d2num",&d2num)) d2num = 0.0; /* numbered tic interval on axis 2 (0.0 for automatic) */
    if (!sf_getfloat("f2num",&f2num)) f2num = 0.0; /* first numbered tic on axis 2 (used if d2num not 0.0) */
    if (!sf_getint("n2tic",&n2tic)) n2tic = 1; /* number of tics per numbered tic on axis 2 */

    if (NULL == (grid2s = sf_getstring("grid2"))) grid2s="none"; /* grid lines on axis 2 (none, dot, dash, or solid) */
    if (STREQ("dot",grid2s)) grid2 = DOT;
    else if (STREQ("dash",grid2s)) grid2 = DASH;
    else if (STREQ("solid",grid2s)) grid2 = SOLID;
    else grid2 = NONE;
    if (NULL == (label2 = sf_getstring("label2"))) label2=""; /* label on axis 2 */

    if (NULL == (labelfont = sf_getstring("labelfont"))) labelfont = "Erg14"; /* font name for axes labels */

    if (!sf_getfloat("labelsize",&labelsize)) labelsize = 18.0; 
    if (NULL == (title = sf_getstring("title"))) title=""; /* title of plot */
    if (NULL == (titlefont = sf_getstring("titlefont"))) titlefont="Rom22"; /* font name for title */
    if (!sf_getfloat("titlesize",&titlesize)) titlesize = 24.0;

    if (NULL == (styles = sf_getstring("style"))) styles="seismic"; 
    /*  normal (axis 1 horizontal, axis 2 vertical) 
	or  seismic (axis 1 vertical, axis 2 horizontal) */
    if (STREQ("normal",styles)) style = NORMAL;
    else style = SEISMIC;

    if (NULL == (titlecolor = sf_getstring("titlecolor"))) titlecolor="red"; /* color for title */
    if (NULL == (labelcolor = sf_getstring("labelcolor"))) labelcolor="blue"; /* color for axes labels */
    if (NULL == (gridcolor = sf_getstring("gridcolor"))) gridcolor="blue"; /* color for grid lines */
    if (NULL == (windowtitle = sf_getstring("windowtitle"))) windowtitle="sfimage"; /* title on window */

    /* adjust x1beg and x1end to fall on sampled values */
    i1beg = NINT((x1beg-f1)/d1);
    i1beg = SF_MAX(0,SF_MIN(n1-1,i1beg));
    x1beg = f1+i1beg*d1;
    i1end = NINT((x1end-f1)/d1);
    i1end = SF_MAX(0,SF_MIN(n1-1,i1end));
    x1end = f1+i1end*d1;

    /* adjust x2beg and x2end to fall on sampled values */
    i2beg = NINT((x2beg-f2)/d2);
    i2beg = SF_MAX(0,SF_MIN(n2-1,i2beg));
    x2beg = f2+i2beg*d2;
    i2end = NINT((x2end-f2)/d2);
    i2end = SF_MAX(0,SF_MIN(n2-1,i2end));
    x2end = f2+i2end*d2;

    /* allocate space for image bytes */
    n1c = 1+abs(i1end-i1beg);
    n2c = 1+abs(i2end-i2beg);
    cz = sf_ucharalloc(n1c*n2c);

    /* convert data to be imaged into signed characters */
    zscale = (wclip!=bclip)?255.0/(wclip-bclip):1.0e10;
    zoffset = -bclip*zscale;
    i1step = (i1end>i1beg)?1:-1;
    i2step = (i2end>i2beg)?1:-1;
    if (style==NORMAL) {
	for (i2c=0,i2=i2beg; i2c<n2c; i2c++,i2+=i2step) {
	    czp = cz+n1c*n2c-(i2c+1)*n1c;
	    for (i1c=0,i1=i1beg; i1c<n1c; i1c++,i1+=i1step) {
		zi = zoffset+z[i1+i2*n1]*zscale;
		if (zi<0.0) zi = 0.0;
		if (zi>255.0) zi = 255.0;
		*czp++ = (unsigned char)zi;
	    }
	}
    } else {
	czp = cz;
	for (i1c=0,i1=i1beg; i1c<n1c; i1c++,i1+=i1step) {
	    for (i2c=0,i2=i2beg; i2c<n2c; i2c++,i2+=i2step) {
		zi = zoffset+z[i1+i2*n1]*zscale;
		if (zi<0.0) zi = 0.0;
		if (zi>255.0) zi = 255.0;
		*czp++ = (unsigned char)zi;
	    }
	}
    }
/*	free1float(z);      keep data for plotting GK */
	
    /* initialize zoom box parameters */
    dx = (style==NORMAL ? d1 : d2);
    dy = (style==NORMAL ? d2 : d1);
    nxb = nx = (style==NORMAL ? n1c : n2c);
    nyb = ny = (style==NORMAL ? n2c : n1c);
    ixb = iyb = 0;
    czb = cz;
    x1begb = x1beg;	 x1endb = x1end;
    x2begb = x2beg;	 x2endb = x2end;

    /* connect to X server */
    if ((dpy=XOpenDisplay(NULL))==NULL)
	sf_error("Cannot connect to display %s!",XDisplayName(NULL));
    scr = DefaultScreen(dpy);
    black = BlackPixel(dpy,scr);
    white = WhitePixel(dpy,scr);
	
    /* create window */
    win = xNewWindow(dpy,xbox,ybox,wbox,hbox,(int) black,(int) white,windowtitle);

    /* backwards compatibility */
    if (STREQ(cmap,"gray")) {
	strcpy(cmap,"rgb0");
    } else if (STREQ(cmap,"hue")) {
	free(cmap);
	cmap = sf_charalloc(5);
	strcpy(cmap,"hsv1");
    } else  if ((strncmp(cmap,"hsv",3)) && (strncmp(cmap,"rgb",3))){
	if (verbose) sf_warning ("cmap=%s using cmap=gray", cmap);

	free(cmap);
	cmap = sf_charalloc(5);
	strcpy(cmap,"rgb0");
    } 
	
    /* here are the new colormaps				*/
    if (strncmp(cmap, "rgb", 3) == 0)
	XSetWindowColormap(dpy,win,
			   xCreateRGBColormap(dpy,win, cmap, verbose));
    else if (strncmp (cmap, "hsv", 3) == 0)
	XSetWindowColormap(dpy,win,
			   xCreateHSVColormap(dpy,win, cmap, verbose));
	
    /* determine min and max pixels from standard colormap */
    pmin = xGetFirstPixel(dpy);
    pmax = xGetLastPixel(dpy);
    if(pmax==0L)pmax=255L;

    if (verbose) sf_warning("pmin=%x pmax=%x",pmin,pmax);
    data_legend = sf_ucharalloc(lwidth * lheight);

    if( bclip < wclip ){
	base=256;
	fact=-256.0;
    }else{
	base=0;
	fact=256.0;
    }
    ptr = data_legend;
    for (i=0; i<lheight; i++){
	for( j=0; j<lwidth; j++ ){
	    *ptr++ = (unsigned char) 
		(base + (fact*i)/lheight);
	}
	/* fprintf(stderr," %d ",*(ptr-1) ); */
    }
		
    /* make GC for image */
    gci = XCreateGC(dpy,win,0,NULL);
	
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
					
    /* determine good size for axes box */
    xSizeAxesBox(dpy,win,
		 labelfont,titlefont,style,
		 &x,&y,&width,&height);
	
    /* clear the window */
    XClearWindow(dpy,win);
	
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
		czbi = newInterpBytes(nxb,nyb,czb,
				      width,height,blockinterp);

		if (image!=NULL) XDestroyImage(image);
		image = xNewImage(dpy,pmin,pmax,
				  width,height,blank,czbi);

		/* BEREND create image */
		if (legend) {
		    if (image_legend!=NULL) XDestroyImage(image_legend);
		    image_legend = xNewImage(dpy,pmin,pmax,lwidth,lheight,0,data_legend);
		}

		imageOutOfDate = 0;
	    }
	
	    /* draw image (before axes so grid lines visible) */
	    XPutImage(dpy,win,gci,image,0,0,x,y,
		      image->width,image->height);

	    /* BEREND display image */
	    if (legend)
		XPutImage(dpy,win,gci,image_legend,
			  0,0,lx,y+ly,lwidth,lheight);

	    /* BEREND draw legend axes on top of image */
	    if (legend)
		xDrawLegendBox(dpy,win,
			       lx,y+ly,lwidth,lheight,
			       bclip,wclip,units,legendfont,
			       labelfont,title,titlefont,
			       labelcolor,titlecolor,gridcolor,
			       style);

	    /* draw curve on top of image */
	    for (i=0; i<curve; i++)
		xDrawCurve(dpy,win,
			   x,y,width,height,
			   x1begb,x1endb,0.0,0.0,
			   x2begb,x2endb,0.0,0.0,
			   x1curve[i],x2curve[i],npair[i],
			   curvecolor[i],style);

	    /* draw axes on top of image */
	    xDrawAxesBox(dpy,win,
			 x,y,width,height,
			 x1begb,x1endb,0.0,0.0,
			 d1num,f1num,n1tic,grid1,label1,
			 x2begb,x2endb,0.0,0.0,
			 d2num,f2num,n2tic,grid2,label2,
			 labelfont,title,titlefont,
			 labelcolor,titlecolor,gridcolor,
			 style);

	    /* else if key down */
	} else if (event.type==KeyPress) {

	    XLookupString(&(event.xkey),keybuf,0,&keysym,&keystat);

	    /*  added moving, clipping and zooming GK */
	    if (keysym==XK_s) {
		xMousePrint(event,style, mpicksfp,
			    x,y,width,height,
			    x1begb,x1endb,x2begb,x2endb);

	    } else if (keysym==XK_l ) {
		/* set lock */		  
		lock = 1 ;
		if (verbose) sf_warning("zoom lock set  %d",lock);

	    } else if (keysym==XK_u ) {
		/* unset lock */		  
		lock = 0 ;
		if (verbose) sf_warning("zoom lock released %d",lock);

	    } else if (keysym==XK_Shift_L ) { 
		/* if (verbose) 
		   sf_warning("Shift Left pressed");*/
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
	    } else if (keysym==XK_Left ) {
		/* move zoom box to left by half window width */
		mve = (x2endb - x2begb)/mvefac ;
		x2begb = x2begb - mve ;
		x2endb = x2endb - mve ;
		msg="move "; 
		/* check for bounds of full window */
		if (x2begb < x2beg){
		    if ( lock ) { x2begb = x2begb + mve ;
			x2endb = x2endb + mve ;
			msg="limit ";
			mve=0;
		    } else {  x2begb = x2beg ;
			nxb=(int)((x2endb-x2begb)/dx);
		    }
		}

		if (verbose) sf_warning("%s %g",msg,mve);

		ixb+=-(int)(mve/dx);
		if ( (ixb<0) || 
		     ((ixb+nxb)>nx) || 
		     (nxb<2) || 
		     (nxb>nx)) {ixb=0;nxb=nx;
		    x2begb=x2beg;
		    x2endb=x2end;}

		if (czb!=cz) free(czb);
		czb = (unsigned char*) sf_charalloc(nxb*nyb);
		for (i=0,czbp=czb; i<nyb; i++) {
		    czp = cz+(iyb+i)*nx+ixb;
		    for (j=0; j<nxb; j++)
			*czbp++ = *czp++; 
		}						
			  
		/* clear area and force an expose event */
		XClearArea(dpy,win,0,0,0,0,True);
		/* note that image is out of date */
		imageOutOfDate = 1;
								
	    } else if (keysym==XK_Right ) {
		/* move zoom box to right by half window width*/
		mve = (x2endb - x2begb)/mvefac ;
		x2begb = x2begb + mve ;
		x2endb = x2endb + mve ;
		msg="move "; 
		/* check for bounds of full window */
		if (x2endb > x2end){
		    if ( lock ) { x2begb = x2begb - mve ;
			x2endb = x2endb - mve ;
			msg="limit ";
			mve=0;
		    } else { x2endb = x2end;
			nxb=(int)((x2endb-x2begb)/dx);
		    }
		}
		if (verbose) sf_warning("%s %g",msg,mve);
			  
		/* for replot require 
		 * ixb,iyb   start samples of image
		 * nxb,nyb   number of samples of image */
			   
		ixb+=(int)(mve/dx);
		if ( (ixb<0) || 
		     ((ixb+nxb)>nx) || 
		     (nxb<2) || 
		     (nxb>nx)) {ixb=0;nxb=nx;
		    x2begb=x2beg;
		    x2endb=x2end;}


		if (czb!=cz) free(czb);
		czb = (unsigned char*) sf_charalloc(nxb*nyb);
		for (i=0,czbp=czb; i<nyb; i++) {
		    czp = cz+(iyb+i)*nx+ixb;
		    for (j=0; j<nxb; j++)
			*czbp++ = *czp++; 
		}						
			  
	
		/* clear area and force an expose event */
		XClearArea(dpy,win,0,0,0,0,True);
		/* note that image is out of date */
		imageOutOfDate = 1;
								
	    } else if (keysym==XK_Down ) {
		/* move zoom box down by half window height */
		mve = (x1endb - x1begb)/mvefac ;
		x1begb = x1begb + mve ;
		x1endb = x1endb + mve ;
		msg="move "; 
		/* check for bounds of full window */
		if (x1endb > x1end){
		    if ( lock ) { x1begb = x1begb - mve ;
			x1endb = x1endb - mve ;
			msg="limit ";
			mve=0;
		    } else { x1endb = x1end;
			nyb=(int)((x1endb-x1begb)/dy);
		    }
		}
		if (verbose) sf_warning("%s %g",msg,mve);

		iyb+=(int)(mve/dy);
			   
		/* reset to original if out of range */
		if ( (iyb<0) || 
		     ((iyb+nyb)>ny) || 
		     (nyb<2) || 
		     (nyb>ny)) {iyb=0;nyb=ny;
		    x1begb=x1beg;
		    x1endb=x1end;}


		/* clear area and force an expose event */
		XClearArea(dpy,win,0,0,0,0,True);
		/* note that image is out of date */
		imageOutOfDate = 1;

		if (czb!=cz) free(czb);
		czb = (unsigned char*) sf_charalloc(nxb*nyb);
		for (i=0,czbp=czb; i<nyb; i++) {
		    czp = cz+(iyb+i)*nx+ixb;
		    for (j=0; j<nxb; j++)
			*czbp++ = *czp++; 
		}						

	    } else if (keysym==XK_Up || keysym==XK_KP_Up ) {
		/*********** 
		 * move zoom box up in .... vertical* 
		 ***********                          */
		mve = (x1endb - x1begb)/mvefac ;
		x1begb = x1begb - mve ;
		x1endb = x1endb - mve ;
		msg="move "; 
		/* check for bounds of full window */
		if (x1begb < x1beg){
		    if ( lock ) { x1begb = x1begb + mve ;
			x1endb = x1endb + mve ;
			msg="limit ";
			mve=0;
		    } else { x1begb = x1beg ;
			nyb=(int)((x1endb-x1begb)/dy);
		    }
		}
		if (verbose) sf_warning("%s %g",msg,mve);

		iyb+=-(int)(mve/dy);
			  
		/* reset to original if out of range */
		if ( (iyb<0) || 
		     (nyb<2) || 
		     (nyb>ny)) {iyb=0;nyb=ny;
		    x1begb=x1beg;
		    x1endb=x1end;}

		if (czb!=cz) free(czb);
		czb = (unsigned char*) sf_charalloc(nxb*nyb);
		for (i=0,czbp=czb; i<nyb; i++) {
		    czp = cz+(iyb+i)*nx+ixb;
		    for (j=0; j<nxb; j++)
			*czbp++ = *czp++; 
		}						
				
		/* clear area and force an expose event */
		XClearArea(dpy,win,0,0,0,0,True);
			
		/* note that image is out of date */
		imageOutOfDate = 1;
									
	    } else if (keysym==XK_o || keysym==XK_KP_Subtract ) {
		/*********** 
		 *zoom out .... vertical* 
		 ***********            */
		mve = (x1endb - x1begb)/mvefac ;
		x1begb = x1begb - mve ;
		x1endb = x1endb + mve ;
		/* check for bounds of full window */
		if (x1begb < x1beg){
		    if ( lock ) { x1begb = x1begb + mve ;
			msg="limit ";
			mve=0;
		    } else { x1begb = x1beg ;}
		}
		if (x1endb > x1end){
		    if ( lock ) { x1endb = x1endb - mve ;
			msg="limit ";
			mve=0;
		    } else { x1endb = x1end ;}
		}
		nyb=(int)((x1endb-x1begb)/dy);
		iyb+=-(int)(mve/dy);
		if ( (iyb<0) || (nyb>ny)) {iyb=0;nyb=ny;}
			  
		/*   .... and horizontal */
		mve = (x2endb - x2begb)/mvefac ;
		x2begb = x2begb - mve ;
		x2endb = x2endb + mve ;
		/* check bounds of original image */
		if (x2begb < x2beg){
		    if ( lock ) { x2begb = x2begb + mve ;
			msg="limit ";
			mve=0;
		    } else { x2begb = x2beg ;}
		}
		if (x2endb > x2end){
		    if ( lock ) { x2endb = x2endb - mve ;
			msg="limit ";
			mve=0;
		    } else { x2endb = x2end ;}
		}
		nxb=(int)((x2endb-x2begb)/dx);
		ixb+=-(int)(mve/dx);
		if ( (ixb<0)        || 
		     ((ixb+nxb)>nx) || 
		     (nxb<0)        || 
		     (nxb>nx))  { ixb=0;nxb=nx;
		    x2begb=x2beg;
		    x2endb=x2end;}
			   
		if (czb!=cz) free(czb);
		czb = (unsigned char*) sf_charalloc(nxb*nyb);
		for (i=0,czbp=czb; i<nyb; i++) {
		    czp = cz+(iyb+i)*nx+ixb;
		    for (j=0; j<nxb; j++)
			*czbp++ = *czp++; 
		}			 	
		/* clear area and force an expose event */
		XClearArea(dpy,win,0,0,0,0,True);
			 
		/* note that image is out of date */
		imageOutOfDate = 1;
								
	    } else if (keysym==XK_i || keysym==XK_KP_Add ) {
		/*********** 
		 *zoom in .... vertical* 
		 ***********           */
		mve = (x1endb - x1begb)/(2.*mvefac) ;
		x1begb = x1begb + mve ;
		x1endb = x1endb - mve ;
		iyb+=(int)(mve/dy);

		/*   .... and horizontal */
		mve = (x2endb - x2begb)/(2.*mvefac) ;
		x2begb = x2begb + mve ;
		x2endb = x2endb - mve ;
		ixb+=(int)(mve/dx);

		nxb=(int)((x2endb-x2begb)/dx);
		nyb=(int)((x1endb-x1begb)/dy);
		if ( (ixb<0) || 
		     (nxb>nx)||
		     (ixb>nx)||
		     (nxb<0) ) {ixb=0;nxb=nx;
		    x2begb=x2beg;
		    x2endb=x2end;}
		if ( (iyb<0) || 
		     (nyb>ny)||
		     (iyb>ny)||
		     (nyb<0) ) {iyb=0;nyb=ny;
		    x1begb=x1beg;
		    x1endb=x1end;}

		/* clear area and force an expose event */
		XClearArea(dpy,win,0,0,0,0,True);
			 
		/* note that image is out of date */
		imageOutOfDate = 1;
		if (czb!=cz) free(czb);
		czb = (unsigned char*) sf_charalloc(nxb*nyb);
		for (i=0,czbp=czb; i<nyb; i++) {
		    czp = cz+(iyb+i)*nx+ixb;
		    for (j=0; j<nxb; j++)
			*czbp++ = *czp++; 
		}					
	    } else if (keysym==XK_c || keysym==XK_Page_Down) {
		  		
		/* Change clip for image */
		clip += clip/10. ;
		if (verbose) sf_warning("clip=%g",clip);
		/* note that image is out of date */
		imageOutOfDate = 1;				
				 
	    } else if (keysym==XK_a || keysym==XK_Page_Up) {

		/* Change clip for image */
		clip -= clip/10. ;
		if (verbose) sf_warning("clip=%g",clip);
		/* note that image is out of date */
		imageOutOfDate = 1;
				
		if (czb!=cz) free(czb);
		czb = (unsigned char*) sf_charalloc(nxb*nyb);
		for (i=0,czbp=czb; i<nyb; i++) {
		    czp = cz+(iyb+i)*nx+ixb;
		    for (j=0; j<nxb; j++)
			*czbp++ = *czp++; 
		}				
		/* end of section for moving clipping and zooming GK */		    

	    } else if (keysym==XK_q || keysym==XK_Q) {
		/* This is the exit from the event loop */
		break;

	    } else if (keysym==XK_r) {
		Colormap mycp=xCreateRGBColormap(dpy,win,"rgb_up",verbose);

		XSetWindowColormap(dpy,win,mycp);
		XInstallColormap(dpy,mycp);
		/* clear area and force an expose event */
		XClearArea(dpy,win,0,0,0,0,True);
		/* note that image is out of date */
		imageOutOfDate = 1;


	    } else if (keysym==XK_R) {
		Colormap mycp=xCreateRGBColormap(dpy,win,"rgb_down",verbose);

		XSetWindowColormap(dpy,win,mycp);
		XInstallColormap(dpy,mycp);

		/* clear area and force an expose event */
		XClearArea(dpy,win,0,0,0,0,True);
		/* note that image is out of date */
		imageOutOfDate = 1;

	    } else if (keysym==XK_h) {
		Colormap mycp=xCreateHSVColormap(dpy,win,"hsv_up",verbose);

		XSetWindowColormap(dpy,win,mycp);
		XInstallColormap(dpy,mycp);

		/* clear area and force an expose event */
		XClearArea(dpy,win,0,0,0,0,True);
		/* note that image is out of date */
		imageOutOfDate = 1;


	    } else if (keysym==XK_H) {

		Colormap mycp=xCreateHSVColormap(dpy,win,"hsv_down",verbose);

		XSetWindowColormap(dpy,win,mycp);
		XInstallColormap(dpy,mycp);

		/* clear area and force an expose event */
		XClearArea(dpy,win,0,0,0,0,True);
		/* note that image is out of date */
		imageOutOfDate = 1;

	    } else {
		continue;
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
		    nxb = nx;
		    nyb = ny;
		    ixb = iyb = 0;
		    if (czb!=cz) free(czb);
		    czb = cz;
			
		    /* else, if new box has non-zero width */
		    /* and height */
		} else {
			
		    /* calculate new box parameters */
		    if (style==NORMAL) {
			zoomBox(x,y,width,height,
				xb,yb,wb,hb,
				nxb,ixb,x1begb,x1endb,
				nyb,iyb,x2endb,x2begb,
				&nxb,&ixb,&x1begb,&x1endb,
				&nyb,&iyb,&x2endb,&x2begb);
		    } else {
			zoomBox(x,y,width,height,
				xb,yb,wb,hb,
				nxb,ixb,x2begb,x2endb,
				nyb,iyb,x1begb,x1endb,
				&nxb,&ixb,&x2begb,&x2endb,
				&nyb,&iyb,&x1begb,&x1endb);
		    }
			
		    /* make new bytes in zoombox */
		    if (czb!=cz) free(czb);
		    czb = (unsigned char*) sf_charalloc(nxb*nyb);
		    for (i=0,czbp=czb; i<nyb; i++) {
			czp = cz+(iyb+i)*nx+ixb;
			for (j=0; j<nxb; j++)
			    *czbp++ = *czp++; 
		    }
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
			  x2begb,x2endb,z,f1,d1,n1,f2,d2,
			  n2,verbose);
	    } else if (event.xbutton.button==Button3) {
		/* ZM: Mouse-Button-3 shows the amplitude constantly */

		showloc = 1;
		xMouseLoc(dpy,win,event,style,showloc,
			  x,y,width,height,x1begb,x1endb,
			  x2begb,x2endb,z,f1,d1,n1,f2,d2,
			  n2,verbose);

	    } else {
		continue;
	    }

	    /* else if pointer has moved */
	} else if (event.type==MotionNotify) {
			
	    /* if button2 down, show mouse location */
	    if (showloc)
		xMouseLoc(dpy,win,event,style,True,
			  x,y,width,height,x1begb,x1endb,
			  x2begb,x2endb,z,f1,d1,n1,f2,d2,
			  n2,verbose);

	    /* else if button2 released, stop tracking */
	} else if (event.type==ButtonRelease &&
		   event.xbutton.button==Button2) {
	    showloc = 0;
	}

    } /* end of event loop */

    /* close connection to X server */
    XCloseDisplay(dpy);

    exit(0);
}

/* update parameters associated with zoom box */
static void zoomBox (int x, int y, int w, int h, 
		     int xb, int yb, int wb, int hb,
		     int nx, int ix, float x1, float x2,
		     int ny, int iy, float y1, float y2,
		     int *nxb, int *ixb, float *x1b, float *x2b,
		     int *nyb, int *iyb, float *y1b, float *y2b)
{
    /* if width and/or height of box are zero, just copy values */
    if (wb==0 || hb==0) {
	*nxb = nx; *ixb = ix; *x1b = x1; *x2b = x2;
	*nyb = ny; *iyb = iy; *y1b = y1; *y2b = y2;
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
	
    /* determine number of samples in rubber box (at least 2) */
    *nxb = SF_MAX(nx*wb/w,2);
    *nyb = SF_MAX(ny*hb/h,2);
	
    /* determine indices of first samples in box */
    *ixb = ix+(xb-x)*(nx-1)/w;
    *ixb = SF_MIN(*ixb,ix+nx-*nxb);
    *iyb = iy+(yb-y)*(ny-1)/h;
    *iyb = SF_MIN(*iyb,iy+ny-*nyb);
	
	
    /* determine box limits to nearest samples */
    *x1b = x1+(*ixb-ix)*(x2-x1)/(nx-1);
    *x2b = x1+(*ixb+*nxb-1-ix)*(x2-x1)/(nx-1);
    *y1b = y1+(*iyb-iy)*(y2-y1)/(ny-1);
    *y2b = y1+(*iyb+*nyb-1-iy)*(y2-y1)/(ny-1);
}

/* return pointer to new interpolated array of bytes */
static unsigned char *newInterpBytes (int n1in, int n2in, unsigned char *bin,
				      int n1out, int n2out, int useBlockInterp) /* JG */
{
    unsigned char *bout;
    float d1in,d2in,d1out,d2out,f1in,f2in,f1out,f2out;
	
    f1in = f2in = f1out = f2out = 0.0;
    d1in = d2in = 1.0;
    d1out = d1in*(float)(n1in-1)/(float)(n1out-1);
    d2out = d2in*(float)(n2in-1)/(float)(n2out-1);
    bout = sf_ucharalloc(n1out*n2out);
    /* JG .... */
    if (!useBlockInterp)
    {
	intl2b(n1in,d1in,f1in,n2in,d2in,f2in,bin,
	       n1out,d1out,f1out,n2out,d2out,f2out,bout);
    }
    else
    {
	intl2b_block(n1in,d1in,f1in,n2in,d2in,f2in,bin,
		     n1out,d1out,f1out,n2out,d2out,f2out,bout);
    }
    /* .... JG */
    return bout;
}

void xMouseLoc(Display *dpy, Window win, XEvent event, int style, Bool show,
	       int x, int y, int width, int height,
	       float x1begb, float x1endb, float x2begb, float x2endb,
	       float *z, float f1, float d1, int n1,float f2, float d2, 
	       int n2, int verbose)
{
    static XFontStruct *fs=NULL;
    static XCharStruct overall;
    static GC gc;
    int dummy,xoffset=5,yoffset=5;
    float x1,x2,amp;
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
	x2 = x2endb+(x2begb-x2endb)*(event.xmotion.y-y)/height;
    } else {
	x1 = x1begb+(x1endb-x1begb)*(event.xmotion.y-y)/height;
	x2 = x2begb+(x2endb-x2begb)*(event.xmotion.x-x)/width;
    }

    /* draw string indicating mouse location */
       
    /* ZM: computing amplitude at the poked point from dataset */
    amp = getamp(z,f1,d1,n1,f2,d2,n2,x1,x2,verbose);
    sprintf(string,"(%0.6g,%0.6g,%0.6g)",x2,x1,amp); /* ZM */

    XTextExtents(fs,string,(int)strlen(string),&dummy,&dummy,&dummy,&overall);
    XDrawString(dpy,win,gc,xoffset,yoffset+overall.ascent,
		string,(int) strlen(string));
}

void xMousePrint(XEvent event, int style, FILE *mpicksfp,
		 int x, int y, int width, int height,
		 float x1begb, float x1endb, float x2begb, float x2endb)
{
    float x1,x2;

    /* convert mouse location to (x1,x2) coordinates */
    if (style==NORMAL) {
	x1 = x1begb+(x1endb-x1begb)*(event.xmotion.x-x)/width;
	x2 = x2endb+(x2begb-x2endb)*(event.xmotion.y-y)/height;
    } else {
	x1 = x1begb+(x1endb-x1begb)*(event.xmotion.y-y)/height;
	x2 = x2begb+(x2endb-x2begb)*(event.xmotion.x-x)/width;
    }

    /* write string indicating mouse location */
    fprintf(mpicksfp, "%0.6g  %0.6g\n", x1, x2);
}

static float getamp(float *zz,float f1,float d1,int n1,
		    float f2,float d2,int n2,float x1,float x2,int verbose)

/*****************************************************************************
return the amplitude value at x1,x2
******************************************************************************
Input:
zz		zz(n1*nz) is the data
f1		coordinate in first x
d1		x coordinate increment
n1		number samples in x
f2		coordinate in first y
d2		y coordinate increment
n2		number samples in y
x1              x1 coordinate of the probed point
x2              x2 coordinate of the probed point
******************************************************************************
Author: Zhaobo Meng, ConocoPhillips, Feb. 03,2004
*****************************************************************************/
{
    float x1last,x2last,x1i,x2i,xfrac,zfrac,xfrac0,zfrac0,temp;
    int ix1,ix2,i00;

    if (d1==0.0) sf_error("d1 can not be 0.0");
    if (d2==0.0) sf_error("d2 can not be 0.0");

    x1last = f1 + (n1-1)*d1;
    if (x1<f1 || x1>x1last) return -999.0; 	
    x2last = f2 + (n2-1)*d2;
    if (x2<f2 || x2>x2last) return -999.0; 	

    x1i = (x1-f1)/d1;
    ix1 = SF_MIN(x1i,n1-1);
    xfrac = x1i-ix1;
    xfrac0 = 1.0-xfrac;

    x2i = (x2-f2)/d2;
    ix2 = SF_MIN(x2i,n2-1);
    zfrac = x2i-ix2;
    zfrac0 = 1.0-zfrac;

    i00 = ix1 + n1*ix2;
    temp = zfrac *( xfrac*zz[i00+1+n1] + xfrac0*zz[i00+n1])
	+ zfrac0*( xfrac*zz[i00+1   ] + xfrac0*zz[i00   ]);

    return(temp);
}
