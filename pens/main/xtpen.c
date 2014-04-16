/*  Vplot filter for X windows using the X Toolkit (Xt). 

    This pen takes all the standard X-toolkit options
    E.g. -geometry 500x400 -font fixed

    The pen has two display modes

    RUNNING MODE: Runs through all the frames in a loop
    Active buttons are:
    QUIT : quits the program
    STOP : enter frame mode

    FRAME MODE (pause=-1): Pauses after each frame
    Active buttons are:
    NEXT : next frame
    PREV : previous frame
    QUIT : quits the program
    RESTART : go to the first frame
    RUN  : enter running mode
    STRETCHY/RIGID : make plot fill the frame or preserve aspect ratio
    FORWARDS/BACKWARDS/BOTH-WAYS : change direction of frame flipping
    Note that a backwards run will only show those frames already plotted
    It is advisable to run once through all the frames forwards.

    The following actions are available for binding to keystrokes;
    xt_quit(): quit program   xt_next(): next frame   xt_prev(): prev frame
    xt_run(): run mode        xt_stop(): frame mode   xt_restart(): first frame
    xt_faster(): reduce pause between frames in run mode
    xt_slower(): increase pause between frames in run mode
    xt_stretchy(): toggle between stretchy and rigid modes
    xt_number(digit): enter a digit in the current number
    xt_reset_number(): reset the current number
    xt_goto_frame(): goto the frame number given by the current number
    xt_print_coord(): Print mouse coords in the file given by interact=

    The default key bindings are:
    <Btn1Down>:            xt_print_coord()  \n\
    <KeyPress>n:           xt_stop() xt_reset_number() xt_next()  \n\
    <KeyPress>m:           xt_stop() xt_reset_number() xt_prev()  \n\
    <KeyPress>r:           xt_run()  \n\
    <KeyPress>q:           xt_quit()  \n\
    <KeyPress>.:           xt_stop()  \n\
    <KeyPress>f:           xt_faster()  \n\
    <KeyPress>s:           xt_slower()  \n\
    <KeyPress>t:           xt_stretchy()  \n\
    <KeyPress>0:           xt_number(0)  \n\
    ......                  .......
    <KeyPress>9:           xt_number(9)  \n\
    <KeyPress>Return:      xt_goto_frame() xt_reset_number()  \n\
    <KeyPress>Escape:      xt_reset_number()

    Here is an example of overriding these in your ~/.Xdefaults file
    this binds the keypad number 1 to skip to the first frame
    xtpen*translations: #override\n\
    <KeyPress>Q:       xt_quit() \n\
    <KeyPress>KP_1:       xt_number(1) xt_goto_frame() xt_reset_number()

    N.B using prev when you are at the first frame takes you to the last
    frame plotted so far; this is not necessarily the last frame!
    You can only jump to a frame if it has already been plotted once.
    If you give an invalid frame number it will jump to the next frame.

    Many parameters may have their defaults set in the Xresource database
    Here are the equivalent names:
    option name          X-Resource name         Type
    ===========          ===============         ====
    stretchy              XTpen.stretchy         Boolean
    images                XTpen.useImages        Boolean
    pixmaps               XTpen.usePixmaps       Boolean
    buttons               XTpen.showButtons      Boolean
    labels                XTpen.showLabels       Boolean
    want_text             XTpen.showText         Boolean
    numcol                XTpen.numCol           int
    pause                 XTpen.pause            int
 
    E.g. If you want xtpen to come up in stretchy mode as a default
    put this line in your Xdefaults file:
    XTpen.stretchy: True
*/

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
 *  source file:   ./filters/xtlib/appdata.c
 *
 * Dave Nichosl (SEP), April 13 1993
 *      Inserted this sample edit history entry.
 * Dave Nichols (SEP), Jan 14 1994
 *      app_data.pause is a float not an int.
 * Dave Nichols (SEP), Nov 16 1999
 *      Add visual depth as a resource
 */

/*
 *
 *  source file:   ./filters/xtlib/popup.c
 *
 * Dave Nichols (SEP), July 18 1992
 *      Inserted this sample edit history entry.
 *
 * Stewart A. Levin (Mobil), Feb 19, 1993
 *      Removed erroneous ; after actionPointPopupConfirm function body.
 *
 * Stewart A. Levin (SEP), May 22, 1994
 *      Add explicit XtNinput specification to allow use with OpenWindows
 */

/*
 *
 *  source file:   ./filters/xtlib/xtcolors.c
 *
 * Dave Nichols (SEP), December 14 1993
 *      Inserted this sample edit history entry.
 */
/*
 * Joe Dellinger (AMOCO), June 9 1995
 *	This code had a check in the xt_set_color_table routine
 *	which disallowed setting colors 0 through 7 to black.
 *	This was to get around a bug elsewhere in the code
 *	that has been fixed, so the check (and the bug it
 *	created) can now safely be removed.
 * Dave Nichols (Schlumberger) Oct 10 1998
 *      Add support for true color visuals
 * Dave Nichols (Schlumberger) Nov 16 1999
 *      Restore colormaps all at once (for speed)
 */

/*
 *
 *  source file:   ./filters/xtlib/xtcommands.c
 *
 * Steve Cole (SEP), August 4 1991
 *      Inserted this sample edit history entry.
 *
 * Dave Nichols (SEP), February 14 1992
 *      Changed to new control structure, added new buttons
 *      for run, stop, prev.
 * Dave Nichols (SEP), April 30 1992
 *      Final version of the new xtpen, added message window
 * Dave Nichols (SEP), July 18 1992
 *      Added boxy mode for interact, it pops up a dialog to enter
 *	information for the label and it outputs a format suitable
 *	for the command line of Box
 */

/*
 *
 *  source file:   ./filters/xtlib/xt_dovplot.c
 *
 * Steve Cole (SEP), August 4 1991
 *      Inserted this sample edit history entry.
 *
 * Dave Nichols (SEP), February 14 1992
 *      Rewrote control structure, dovplot is now invoked only when
 *      needed and on a per-frame basis.
 * Dave Nichols (SEP), April 30 1992
 *      Final version of the new xtpen, added message window.
 * Dave Nichols (SEP), April 12 1993
 *      Moved xt_size_n_scale into this routine so that user_scales are
 *      set before it is called.
 * Dave Nichols (SEP), April 13 1993
 *      Use app_data to get defaults from resource database
 */

/* 
 * source file:   ./filters/xtlib/xtpaint.c 
 *
 * Steve Cole (SEP), August 4 1991
 *      Inserted this sample edit history entry.
 *
 * Dave Nicols (SEP), February 14 1992
 *      Rewrote using new xtFrame interface that only plots one
 *      frame at a time.
 * Dave Nichols (SEP), April 30 1992
 *      Final version of the new xtpen, added message window
 * Dave Nichols (SEP), Nov 19 1993
 *      Added proper handling of repaints that involve breaks.
 *      Frame structure now contains info on how frame is terminated.
 */

/*
 *
 *  source file:   ./filters/xtlib/xtpixmap.c
 *
 * Steve Cole (SEP), August 4 1991
 *      Inserted this sample edit history entry.
 * Dave Nicols (SEP), February 14 1992
 *      Rewrote using new xtFrame interface that only plots one
 *      frame at a time.
 * Dave Nichols (SEP), April 30 1992
 *      Final version of the new xtpen, added message window
 * Dave Nichols (SEP), May 6 1992
 *      see_progress controls offscreen pixmap.
 */

#include <string.h>
#include <stdio.h>
#include <strings.h>

extern void bcopy (const void *src, void *dst, size_t len);

#include <sys/ioctl.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include <X11/IntrinsicP.h>
#include <X11/Xatom.h>
#include <X11/Object.h>
#include <X11/ObjectP.h>

#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>

#include <X11/cursorfont.h>
#include <X11/Shell.h>

#include <X11/Xaw/Box.h>
#include <X11/Xaw/Paned.h>
#include <X11/Xaw/Form.h>
#include <X11/Xaw/Command.h>
#include <X11/Xaw/Label.h>
#include <X11/Xaw/Dialog.h>
#include <X11/Xaw/Cardinals.h>
#include <X11/Xaw/AsciiText.h>

#include <rsfplot.h>

#include "../include/attrcom.h"
#include "../include/closestat.h"
#include "../include/err.h"
#include "../include/erasecom.h"
#include "../include/intcom.h"
#include "../include/mesgcom.h"

#include "dovplot.h"
#include "init_vplot.h"

#include "../utilities/util.h"
#include "../genlib/genpen.h"

#include "_xt.h"
#include "xtdoc.h"

#include "xtpen.h"

#include "_device.h"

AppData app_data;

#define XtNstretchy "stretchy"
#define XtCStretchy "Stretchy"
#define XtNusePixmaps "usePixmaps"
#define XtCUsePixmaps "UsePixmaps"
#define XtNuseImages    "useImages"
#define XtCUseImages    "UseImages"
#define XtNshowButtons "showButtons"
#define XtCShowButtons "ShowButtons"
#define XtNshowLabels "showLabels"
#define XtCShowLabels "ShowLabels"
#define XtNshowText "showText"
#define XtCShowText "ShowText"
#define XtNnumCol "numCol"
#define XtCNumCol "NumCol"
#define XtNpause "pause"
#define XtCPause "Pause"
#define XtNvisDepth "visDepth"
#define XtCVisDepth "VisDepth"

static XtResource resources[] = {
    {
	XtNstretchy,
	XtCStretchy,
	XtRBoolean,
	sizeof(Boolean),
	XtOffset( AppDataPtr, stretchy ),
	XtRString,
	(XtPointer) "False",
    },
    {
	XtNusePixmaps,
	XtCUsePixmaps,
	XtRBoolean,
	sizeof(Boolean),
	XtOffset( AppDataPtr, pixmaps ),
	XtRString,
	(XtPointer) "False",
    },
    {
	XtNuseImages,
	XtCUseImages,
	XtRBoolean,
	sizeof(Boolean),
	XtOffset( AppDataPtr, images ),
	XtRString,
	(XtPointer) "True",
    },
    {
	XtNshowButtons,
	XtCShowButtons,
	XtRBoolean,
	sizeof(Boolean),
	XtOffset( AppDataPtr, buttons ),
	XtRString,
	(XtPointer) "True",
    },
    {
	XtNshowLabels,
	XtCShowLabels,
	XtRBoolean,
	sizeof(Boolean),
	XtOffset( AppDataPtr, labels ),
	XtRString,
	(XtPointer) "True",
    },
    {
	XtNshowText,
	XtCShowText,
	XtRBoolean,
	sizeof(Boolean),
	XtOffset( AppDataPtr, textpane ),
	XtRString,
	(XtPointer) "True",
    },
    {
	XtNnumCol,
	XtCNumCol,
	XtRInt,
	sizeof(int),
	XtOffset( AppDataPtr, num_col ),
	XtRString,
	(XtPointer) "0",
    },
    { 
	XtNvisDepth,
	XtCVisDepth,
	XtRInt,
	sizeof(int),
	XtOffset( AppDataPtr, vis_depth ),
	XtRString,
	(XtPointer) "32",
    },
    {
	XtNpause,
	XtCPause,
	XtRFloat,
	sizeof(float),
	XtOffset( AppDataPtr, pause ),
	XtRString,
	(XtPointer) "0.",
    }
};

int xt_after_break = 0;
int xt_after_erase = 0;

Display		*pen_display;
GC		pen_gc;
Window		pen_window;
Pixmap		pen_pixmap;
Drawable	pen_drawable;
Visual		*pen_visual;
unsigned long	*pen_colors;
Colormap	pen_colormap;
Widget		xtpen, vpane, control_panel, pen_picture; 
Widget		text_region=NULL;

XtAppContext	pen_context;
int		own_colormap;
int		x_num_col;
int		num_col_req = 0;
int		num_col_max = 128;
int		num_col_min = 16;
int             xmono = NO;
int             xtruecol = NO;
int		screen_depth;
int		visual_depth;
long		screen_black, screen_white;
bool             xt_stretchy;
bool 		boxy,want_images,greedy_pixmaps,see_progress;
const char 	*message_text;
int		plotting_started = 0;

static void  create_pixmap(int width, int height );
static void choose_visual(Display* display, int screen, int max_depth, 
			  XVisualInfo* vinfo);
static void clear_pixmap(void);
static void remove_pixmap(void);
static void xt_size_n_scale(Dimension *width, Dimension *height);
static void PenRepaint (Widget w, XEvent *event, String *params, 
			Cardinal *num_params);

String fallback_resources[] = {
    "*text.scrollVertical:            whenNeeded",
    "*text.autoFill:                  True",
    "*text.Wrap:                      word",
    "*text.height:                    50",
    "*Font:                           fixed",
    "*Foreground:                     black",
    "*Background:                     white",
    NULL,
};

/* 	 <ColormapNotify>:	PenRepaint() \n\
         <ConfigureNotify>:	PenRepaint() \n\ */

/* default translation table for pen_picture widget */

static const char trans[] =
    "<Expose>:		PenRepaint() \n\
         <ConfigureNotify>:	PenRepaint() \n\
         <Btn1Down>:            xt_print_coord() \n\
         None<KeyPress>n:       xt_stop() xt_reset_number() xt_next() \n\
         None<KeyPress>m:       xt_stop() xt_reset_number() xt_prev() \n\
         None<KeyPress>r:       xt_run()  \n\
         None<KeyPress>q:       xt_quit()  \n\
         None<KeyPress>.:       xt_stop()  \n\
         None<KeyPress>f:       xt_faster()  \n\
         None<KeyPress>s:       xt_slower()  \n\
         None<KeyPress>t:       xt_stretchy()  \n\
	 None<KeyPress>Escape: 	xt_reset_number()\n\
         None<KeyPress>0: 	xt_number(0)\n\
         None<KeyPress>1: 	xt_number(1)\n\
         None<KeyPress>2: 	xt_number(2)\n\
         None<KeyPress>3: 	xt_number(3)\n\
         None<KeyPress>4: 	xt_number(4)\n\
         None<KeyPress>5: 	xt_number(5)\n\
         None<KeyPress>6: 	xt_number(6)\n\
         None<KeyPress>7: 	xt_number(7)\n\
         None<KeyPress>8: 	xt_number(8)\n\
         None<KeyPress>9: 	xt_number(9)\n\
	 None<KeyPress>Return:	xt_goto_frame() xt_reset_number()";


static XtActionsRec window_actions[] = {
    {"PenRepaint",	  PenRepaint},
    {"xt_quit",           actionQuit},	
    {"xt_next",           actionNext},	
    {"xt_stretchy",       actionStretch},	
    {"xt_prev",           actionPrev},
    {"xt_run",            actionRun},
    {"xt_stop",           actionStop},
    {"xt_restart",        actionRestart},
    {"xt_faster",         actionFaster},
    {"xt_slower",         actionSlower},
    {"xt_number",         actionNumber},
    {"xt_reset_number",   actionNumReset},
    {"xt_goto_frame",     actionGoto},
    {"xt_print_coord",    actionCoord},
    {"xt_run_mode",       actionRunMode},
};

extern int xt_next_num;

char            name[] = "xtpen";

int xt_app_data(Widget app)
{

    XtVaGetApplicationResources( app,
				 &app_data,
				 resources, XtNumber( resources ),
				 NULL );
    return 0;

}


typedef union
{
    int i;
    float f;
} IntFloatUnion;


static Widget PopUpWidget;
static Widget PointLabelWidget;
static Widget xOffset,yOffset;
static Widget xLabel,yLabel; 
static Widget Color,Fat,Boxit,Size,XOVAL,YOVAL,POINTER;

static float vpx, vpy;
static FILE *outfp;

void PopupConfirm(
    Widget     w         /* widget */,
    XtPointer  x1        /* closure */,  /* data the application registered */
    XtPointer  x2        /* call_data */ /* callback specific data */
    )
{
    float xoff,yoff;
    int   col,fat,box,point;
    float size;
    float x_oval,y_oval;
    char boxit[128],arrow[129];
	

    sscanf(XawDialogGetValueString(xOffset),"%f",&xoff );
    sscanf(XawDialogGetValueString(yOffset),"%f",&yoff );
    sscanf(XawDialogGetValueString(Color),"%d",&col );
    sscanf(XawDialogGetValueString(Fat),"%d",&fat );
    sscanf(XawDialogGetValueString(Size),"%f",&size );
    sscanf(XawDialogGetValueString(Boxit),"%s",boxit );
    sscanf(XawDialogGetValueString(YOVAL),"%f",&y_oval );
    sscanf(XawDialogGetValueString(XOVAL),"%f",&x_oval );
    sscanf(XawDialogGetValueString(POINTER),"%s",arrow );
    if(boxit[0]=='n') box=0;
    else box=1;
    if(arrow[0]=='n') point=0;
    else point=1;

    fprintf( outfp,"x0=%f y0=%f label=\"%s\" xt=%f yt=%f lab_fat=%d lab_color=%d boxit=%d size=%f pointer=%d x_oval=%f y_oval=%f \n",vpx,vpy,
	     XawDialogGetValueString(PointLabelWidget), xoff, yoff,fat,col,box,size,point,x_oval,y_oval);
    fflush( outfp );
    XtPopdown( PopUpWidget );
}

void actionPointPopupConfirm(Widget w, XEvent *event, 
			     String *params, Cardinal *num_params)
{ PopupConfirm(w,NULL,NULL); }


void PopupCancel(
    Widget     w         /* widget */,
    XtPointer  x1        /* closure */,  /* data the application registered */
    XtPointer  x2        /* call_data */ /* callback specific data */
    )
{
    XtPopdown( PopUpWidget );
}

void actionPointPopupCancel(Widget w, XEvent *event, 
			    String *params, Cardinal *num_params)
{ PopupCancel(w,NULL,NULL); }

static void set_labels(float x, float y)
{
    Arg arg[1];
    char text[32];

    sprintf( text, "%6.2f ", x );
    XtSetArg( arg[0], XtNlabel, text );
    XtSetValues( xLabel , arg, ONE );

    sprintf( text, "%6.2f ", y );
    XtSetArg( arg[0], XtNlabel, text );
    XtSetValues( yLabel , arg, ONE );
}


void doPointPopup(FILE* fp, float x, float y )
{
    Arg argList[20];
    int args;
  
    Window root, child;
    static Widget cancelWidget, confirmWidget;
    static Widget vpaneWidget, labelBox, offsetBox,otherBox,otherBox2;

    char xoff_string[32];
    char yoff_string[32];
    char text_string[32];
    char box_string[32];
    char size_string[32];
    char color_string[32];
    char fat_string[32];
    char pointer_string[32];
    char xov_string[32];
    char yov_string[32];
  
    static Boolean inited = False;
    int root_x, root_y, child_x, child_y;
    unsigned int buttons;
  
    Dimension width;
    Dimension height;
    Dimension px,py,pw,ph;

    vpx=x; vpy=y; outfp= fp;

    /*
     * Find out where the mouse is, so we can put the confirmation
     * box right there.
     */
    XQueryPointer( pen_display,pen_window,
		   &root, &child,
		   &root_x, &root_y, &child_x, &child_y, &buttons);
    /*
     * If we need to construct the label box do that,
     * otherwise just reset the position and callbacks and
     * put it up again.
     */
    if (!inited) {
	/*
	 * The confirmation box will be a pop-up widget.
	 */
	args = 0;
	XtSetArg(argList[args],XtNinput, True); args++;
	PopUpWidget =
	    XtCreatePopupShell("Point", transientShellWidgetClass,
			       xtpen, argList, args);
    
	args = 0;
	vpaneWidget = XtCreateWidget("vpane",panedWidgetClass,
				     PopUpWidget,
				     NULL, 0);

	/*
	 * Make a box to put the labels in.
	 */

	labelBox = XtCreateManagedWidget("labelBox", boxWidgetClass,
					 vpaneWidget, NULL, 0);

	xLabel = XtCreateManagedWidget("000.00",labelWidgetClass,labelBox,NULL,0);
	yLabel = XtCreateManagedWidget("000.00",labelWidgetClass,labelBox,NULL,0);

	offsetBox = XtCreateManagedWidget("offsetBox", boxWidgetClass,
					  vpaneWidget, NULL, 0);
    
	/* dialogs for the offsets */
	args = 0;
	strcpy(pointer_string, "yes" );
	XtSetArg( argList[0], XtNvalue, (XtArgVal)pointer_string ); args++;
	XtSetArg( argList[1], XtNlabel, (XtArgVal)"Do arrow" ); args++;
	XtSetArg( argList[2], XtNleft, XtChainLeft ); args++;
	POINTER = XtCreateManagedWidget("dialog1", dialogWidgetClass,
					offsetBox, argList, args);
	args = 0;
	sprintf( xoff_string, "%4.2f", 1.0 );
	XtSetArg( argList[0], XtNvalue, (XtArgVal)xoff_string ); args++;
	XtSetArg( argList[1], XtNlabel, (XtArgVal)"x-off" ); args++;
	XtSetArg( argList[2], XtNfromHoriz, POINTER); args++;
	xOffset = XtCreateManagedWidget("dialog1", dialogWidgetClass,
					offsetBox, argList, args);

	args = 0;
	sprintf( yoff_string, "%4.2f", 1.0 );
	XtSetArg( argList[0], XtNvalue, (XtArgVal)yoff_string ); args++;
	XtSetArg( argList[1], XtNlabel, (XtArgVal)"y-off" ); args++;
	XtSetArg( argList[2], XtNfromHoriz, xOffset); args++;
	yOffset = XtCreateManagedWidget("dialog2", dialogWidgetClass,
					offsetBox, argList, args);
	args = 0;
	sprintf( xov_string, "%4.2f", 0.0 );
	XtSetArg( argList[0], XtNvalue, (XtArgVal)xov_string ); args++;
	XtSetArg( argList[1], XtNlabel, (XtArgVal)"Oval-x" ); args++;
	XtSetArg( argList[2], XtNfromHoriz, yOffset); args++;
	XOVAL = XtCreateManagedWidget("dialog2", dialogWidgetClass,
				      offsetBox, argList, args);
	args = 0;
	sprintf( yov_string, "%4.2f", 0.0 );
	XtSetArg( argList[0], XtNvalue, (XtArgVal)yov_string ); args++;
	XtSetArg( argList[1], XtNlabel, (XtArgVal)"Oval-y" ); args++;
	XtSetArg( argList[2], XtNfromHoriz, XOVAL); args++;
	YOVAL = XtCreateManagedWidget("dialog2", dialogWidgetClass,
				      offsetBox, argList, args);

	/*addition features to set */
	otherBox = XtCreateManagedWidget("otherBox", boxWidgetClass,
					 vpaneWidget, NULL, 0);

	args = 0;
	sprintf( fat_string, "%3d", 0);
	XtSetArg( argList[0], XtNvalue, (XtArgVal)fat_string ); args++;
	XtSetArg( argList[1], XtNlabel, (XtArgVal)"fat" ); args++;
	XtSetArg( argList[2], XtNleft, XtChainLeft ); args++;
	Fat = XtCreateManagedWidget("dialog3", dialogWidgetClass,
				    otherBox, argList, args);

	args = 0;
	sprintf( color_string, "%2d", 7);
	XtSetArg( argList[0], XtNvalue, (XtArgVal)color_string ); args++;
	XtSetArg( argList[1], XtNlabel, (XtArgVal)"color" ); args++;
	XtSetArg( argList[2], XtNfromHoriz, Fat); args++;
	Color = XtCreateManagedWidget("dialog4", dialogWidgetClass,
				      otherBox, argList, args);

	otherBox2 = XtCreateManagedWidget("otherBox", boxWidgetClass,
					  vpaneWidget, NULL, 0);

	args = 0;
	sprintf( size_string, "%6.3f", .25);
	XtSetArg( argList[0], XtNvalue, (XtArgVal)size_string ); args++;
	XtSetArg( argList[1], XtNlabel, (XtArgVal)"Size" ); args++;
	XtSetArg( argList[2], XtNfromHoriz, (XtArgVal) Color); args++;
	Size = XtCreateManagedWidget("dialog5", dialogWidgetClass,
				     otherBox2, argList, args);

	args = 0;
	sprintf( box_string, "%s","yes");
	XtSetArg( argList[0], XtNvalue, (XtArgVal)box_string ); args++;
	XtSetArg( argList[1], XtNlabel, (XtArgVal)"Box" ); args++;
	XtSetArg( argList[2], XtNfromHoriz, (XtArgVal) Size); args++;
	Boxit = XtCreateManagedWidget("dialog6", dialogWidgetClass,
				      otherBox2, argList, args);


    
	/*
	 * Dialog for the label.
	 */
	args = 0;
	strcpy(text_string,"");
	XtSetArg( argList[0], XtNvalue, (XtArgVal)text_string); args++;
	XtSetArg( argList[1], XtNlabel, (XtArgVal)"Label for this point" ); args++;
	PointLabelWidget = XtCreateManagedWidget("dialog", dialogWidgetClass,
						 vpaneWidget, argList, args);
    
	/*
	 * Confirmation button.
	 */
	confirmWidget = XtCreateManagedWidget("confirm", commandWidgetClass,
					      PointLabelWidget , NULL, 0);
	XtAddCallback( confirmWidget ,  XtNcallback, PopupConfirm, NULL );
    
	/*
	 * Cancellation button.
	 */
	cancelWidget = XtCreateManagedWidget("cancel", commandWidgetClass,
					     PointLabelWidget , NULL, 0);
	XtAddCallback( cancelWidget ,  XtNcallback, PopupCancel, NULL );

	/*
	 * Let the shell widget know we're here.
	 */
	XtManageChild(vpaneWidget);
	XtRealizeWidget(PopUpWidget);
	inited = True;
    }

    set_labels( vpx,vpy);
  
    /*
     * Take some pains to center the popup on the pointer, but be certain
     * the thing is visible, else they can never exit
     */
    args = 0;
    XtSetArg(argList[args], XtNx, &px); args++;
    XtSetArg(argList[args], XtNy, &py); args++;
    XtSetArg(argList[args], XtNwidth, &pw); args++;
    XtSetArg(argList[args], XtNheight, &ph); args++;
    XtGetValues(xtpen, argList,args);
  
    args = 0;
    XtSetArg(argList[args], XtNwidth, &width); args++;
    XtSetArg(argList[args], XtNheight, &height); args++;
    XtGetValues(PopUpWidget,argList,args);
  
    root_x -= (width/2);
    root_y -= (3*height/4);
  
    if ( root_x < 0 ) root_x = 0;
    if ( root_y < 0 ) root_y = 0;
  
    args = 0;
    XtSetArg(argList[args], XtNx, root_x); args++;
    XtSetArg(argList[args], XtNy, root_y); args++;
    XtSetValues(PopUpWidget,argList,args);
  
    /* reset the string */
    args = 0;
    strcpy(text_string,"");
    XtSetArg( argList[0], XtNvalue, (XtArgVal)text_string); args++;
    XtSetValues(PointLabelWidget, argList, args);
  
    /* pop it up */
    XtPopup(PopUpWidget, XtGrabExclusive);
}

#define MAXVERT 1000

void xtarea (int npts, struct vertex *head)
/*< area >*/
{
    int             i;
    XPoint vlist[MAXVERT];

    /* translate data structures */
    for (i = 0; i < npts && i < MAXVERT; i++)
    {
	vlist[i].x = XCORD (head->x);
	vlist[i].y = YCORD (head->y);
	head = head->next;
    }
    XDrawLines(pen_display, pen_drawable, pen_gc, vlist, npts, CoordModeOrigin);
    XFillPolygon(pen_display, pen_drawable, pen_gc, vlist, npts, Complex, CoordModeOrigin);
}

void xtattributes (int command, int v, int v1, int v2, int v3)
/*< attributes >*/
{
    XRectangle clip_rect;
    int clip_xorigin,clip_yorigin;

    switch (command)
    {
	case SET_COLOR:
	    xt_set_color( v );
	    break;

	case SET_COLOR_TABLE:
	    xt_set_color_table( v, v1, v2, v3 );
	    break;

	case SET_WINDOW:

	    clip_rect.x = XCORD( v  );
	    clip_rect.y = YCORD( v3 );
	    clip_rect.width = XCORD( v2 ) - XCORD( v ) + 1;
	    clip_rect.height = YCORD( v1 ) - YCORD( v3 ) + 1;
	    clip_xorigin = XCORD(dev.xmin); clip_yorigin =  YCORD(dev.ymax);
	    XSetClipRectangles( pen_display, pen_gc, clip_xorigin, clip_yorigin,
				&clip_rect, 1, Unsorted );
	    break;
	
	case BEGIN_GROUP:
	    break;
	case END_GROUP:
	    if( v == 0 ) {
		/* add the group name to the message window (if it exists ) */
		/*   addText( group_name );*/
		/*   addText( "\n" );*/
		switch( v1 ) {
		    case UNSPECIFIED_EGROUP:
		    case USER_EGROUP:
		    case ERASE_EGROUP:
		    case EOF_ERASE_EGROUP:
			xt_endframe = YES;
			break;
		    case BREAK_EGROUP:
			xt_after_break = 1;
			xt_endframe = YES;
			break;
		    case IGNORE_EGROUP:
		    case EOF_IGNORE_EGROUP:
			xt_endframe = NO; /* eof but not the end of a frame */
			break;
		}
	    }
		
	    break;
    }
}

void xtclose (int status)
/*< close >*/
{
    switch (status)
    {
	case CLOSE_PAUSE:
	    break;

	case CLOSE_NORMAL:
	case CLOSE_ERROR:
	case CLOSE_INTERRUPT:
	    xt_clear_images(frame_list);
	    remove_pixmap();
	    break;

	case CLOSE_FLUSH:
	    /*horrible hack to skip the first two flushes after erase */
	    if( xt_after_erase ){
		xt_after_erase--;
	    }else{
		xt_flush_display(); 
	    }
	    break;
    }
}


#ifndef SEP_MAP_SIZE
#define SEP_MAP_SIZE 65536
#endif
unsigned long    color = 0;	/* current color */

/* black, blue, red, purple, green, cyan, yellow, white */
unsigned char   red[SEP_MAP_SIZE] =   {0, 0,   255, 255, 0,   0,   255, 255};
unsigned char   green[SEP_MAP_SIZE] = {0, 0,   0,   0,   255, 255, 255, 255};
unsigned char   blue[SEP_MAP_SIZE] =  {0, 255, 0,   255, 0,   255, 0,   255};
#ifndef UO_SEP
#define U0_SEP (~((unsigned long) 0))
#endif
unsigned long  map[SEP_MAP_SIZE]  ;


/* routine to find a colormap with at most maxcol and at least mincol
 * available color cells in it
 * It returns the number of cells available and fills in the entries in
 * the array "map".
 */

int xtnumcol(int mincol, int maxcol)
{
    Colormap colormap;
    int ncol ;
    int sucess;
    int i;
    unsigned long planeMask;

    for(i=0; i < SEP_MAP_SIZE; i++) map[i]=U0_SEP;


    /* try to allocate these cells in the current colormap */
    ncol = maxcol;
    colormap = pen_colormap;
    do {
	sucess = XAllocColorCells(pen_display,colormap,
				  0,&planeMask, 0,map,ncol ) ;
	if( !sucess ) ncol = ncol/2;

    } while( !sucess && ncol >=mincol );


    if( ! sucess ){
        /*-- It didn't work so ... construct our own XColormap 
	  We copy the old one so that  we still get black and pwhite
	  defined OK --*/
	ERR( COMMENT,name,"allocated new colormap");
	colormap = XCopyColormapAndFree( pen_display,pen_colormap);
	/* try to allocate these cells in the current colormap */
	ncol = maxcol;
	do{
	    sucess = XAllocColorCells(pen_display,colormap,
				      0,&planeMask, 0,map,ncol ) ;
	    if( !sucess ) ncol = ncol/2;

	} while( !sucess && ncol >=mincol );
    }


    if( !sucess ) 
	ERR( FATAL,name,
	     "unable to obtain minimum number of cells in any colormap");
  
    pen_colormap = colormap; 
    return ncol; 
}


void xt_set_color(int col )
/*< set color >*/
{
    if (mono) {
	if (col){
	    color = map[1];
	}else{
	    color = map[0];
	}
    }else{
	color = map[col];
    }
    XSetForeground(pen_display, pen_gc, color);
}

void xt_set_color_table(int col, int ired, int igreen, int iblue )
/*< set color table >*/
{
    XColor           pen_color;

    if ( mono) return;

    red[col] = ired;
    green[col] = igreen;
    blue[col] = iblue;

    /* put the color in the colormap, we specify the pixel value from the
       values in the table "map".
    */
    
    pen_color.red = red[col] << 8;
    pen_color.green = green[col] << 8;
    pen_color.blue = blue[col] << 8;
    pen_color.flags = DoRed | DoGreen | DoBlue;

    if( xtruecol ) {
	/* True-Color we get told the pixel */
	XAllocColor(pen_display,pen_colormap,&pen_color);
	map[col] = pen_color.pixel;
	
    }else{
	/* Pesudo-Color we suppy the pixel */
        pen_color.pixel = map[col];
	XStoreColor(pen_display,pen_colormap,&pen_color);
    }
	
    if(col == 0) { /* change of background color */
    	XSetWindowBackground( pen_display, pen_window, map[0]);
    }
}

/* save the current colormap in the structure associated with this frame */
void xt_save_colormap(xtFrame *frame)
{
    int i;

    if( mono ) return;

    for( i=0; i < dev.num_col ; i++ ){
	frame->cmap.map[i] = map[i];
	frame->cmap.red[i] = red[i];
	frame->cmap.green[i] = green[i];
	frame->cmap.blue[i] = blue[i];
    }

}

/* restore the saved  current colormap from this frame */
void xt_restore_colormap(xtFrame *frame)
{
    XColor*           pen_color;
    int i;

    if( mono ) return;

    pen_color = (XColor*)XtMalloc( dev.num_col*sizeof(XColor) );

    for( i=0; i < dev.num_col ; i++ ){
        map[i] = frame->cmap.map[i];
        red[i] = frame->cmap.red[i] ;
        green[i] = frame->cmap.green[i];
        blue[i] = frame->cmap.blue[i];
	if( ! xtruecol ) {
	    pen_color[i].red = frame->cmap.red[i] << 8;
	    pen_color[i].green = frame->cmap.green[i] << 8;
	    pen_color[i].blue = frame->cmap.blue[i] << 8;
	    pen_color[i].pixel = frame->cmap.map[i];
	    pen_color[i].flags = DoRed | DoGreen | DoBlue;
	}
    }

    if( ! xtruecol ) {
	XStoreColors(pen_display,pen_colormap,pen_color, dev.num_col);
    }
	
    XSetWindowBackground( pen_display, pen_window, map[0]);

    XtFree((char*) pen_color );
}

/* variable that control the existance of buttons and labels */
static int wantButtons = 0;
static int  first_time = YES;
extern char *interact;
static int wantLabels = 0;

/* The buttons that go into the control panel */
Widget quit_button, next_button, prev_button,  restart_button;
Widget run_button, stop_button, stretch_button;
Widget mode_button, frame_label, delay_label;


/* static variables that are set when a button is pressed */
static int didNEXT = NO, didPREV = NO, didQUIT = NO, didTIME=NO; 
static int didFRAM1 = NO,  didRUN = NO, didSTOP=NO , didCHANGE = NO;


/* Function active? flags. We can't use XtIsSensitive because we might
 * not have active buttons */
static int next_on=NO ,prev_on=NO ,quit_on=NO ,restart_on=NO ;
static int run_on=NO , stop_on=NO, size_on=NO; 

static void  dummy_proc () { return; }

void actionDummy(Widget w, XEvent ev, String *p, Cardinal *np)
/*< dummy >*/
{ dummy_proc(); }

static void  next_proc(
    Widget     w         /* widget */,
    XtPointer  x1        /* closure */,  /* data the application registered */
    XtPointer  x2        
    ){  if( next_on == YES ) didNEXT = YES;    return; }

void actionNext(Widget w, XEvent *ev, String *p, Cardinal *np)
/*< next >*/
{ next_proc(w,NULL,NULL); }

static void  prev_proc(
    Widget     w         /* widget */,
    XtPointer  x1        /* closure */,  /* data the application registered */
    XtPointer  x2        
    ){  if( prev_on==YES ) didPREV = YES;    return; }

void actionPrev(Widget w, XEvent *ev, String *p, Cardinal *np)
/*< prev >*/
{ prev_proc(w,NULL,NULL); }

static void quit_proc(
    Widget     w         /* widget */,
    XtPointer  x1        /* closure */,  /* data the application registered */
    XtPointer  x2        
    ){  if( quit_on == YES ) didQUIT = YES;    return; }

void actionQuit(Widget w, XEvent *ev, String *p, Cardinal *np)
/*< quit >*/
{
    if( interact[0] != '\0' ) {
	if (first_time == YES){
	    outfp = fopen (interact, "w");
	    if (outfp == NULL) { 
		ERR (FATAL, name, "Can't open interact output file %s!", interact);
	    }
	    first_time = NO;
	}
/*	 fprintf( outfp,"\n");*/
	fflush( outfp );
    }
    quit_proc(w,NULL,NULL); 
    
}

static void restart_proc(
    Widget     w         /* widget */,
    XtPointer  x1        /* closure */,  /* data the application registered */
    XtPointer  x2        
    ){  if( restart_on == YES) didFRAM1 = YES;    return; }

void actionRestart(Widget w, XEvent *ev, String *p, Cardinal *np)
/*< restart >*/
{ restart_proc(w,NULL,NULL); }

static void run_proc(
    Widget     w         /* widget */,
    XtPointer  x1        /* closure */,  /* data the application registered */
    XtPointer  x2
    ){  if( run_on == YES ) didRUN = YES;    return; }

void actionRun(Widget w, XEvent *ev, String *p, Cardinal *np)
/*< run >*/
{ run_proc(w,NULL,NULL); }

static void stop_proc(
    Widget     w         /* widget */,
    XtPointer  x1        /* closure */,  /* data the application registered */
    XtPointer  x2        
    ){  if( stop_on == YES ) didSTOP = YES;    return; }

void actionStop(Widget  w, XEvent *ev, String *p, Cardinal *np)
/*< stop >*/
{ stop_proc(w,NULL,NULL); }

static Arg rigidArgs[] = { { XtNlabel, (XtArgVal)"Rigid" } };
static Arg stretchArgs[] = { { XtNlabel, (XtArgVal)"Stretchy" } };

static void stretch_proc(
    Widget     w         /* widget */,
    XtPointer  x1        /* closure */,  /* data the application registered */
    XtPointer  x2  
    ){ 
    xt_stretchy = (bool) !xt_stretchy;  /* toggle the "stretchy" attribute */

    if( wantButtons ){

	if( xt_stretchy ){
	    /* if we are now stretchy change the label to "rigid" */
	    XtSetValues( stretch_button, rigidArgs, XtNumber(rigidArgs) );
	}else{
	    /* if we are now rigid change the label to "stretchy" */  
	    XtSetValues( stretch_button, stretchArgs, XtNumber(rigidArgs) );
	}

    }
    didCHANGE = YES; return; 
}

void actionStretch(Widget w, XEvent *ev, String *p, Cardinal *np)
/*< stretch >*/
{ stretch_proc(w,NULL,NULL); }

/* the timeout proc called in run mode by an application timer */
static void timeout_proc(
    XtPointer      x1,     /* closure */
    XtIntervalId*  id1     /* id */ 
    ){    didTIME = YES; return; }

static void set_delay_label(float delay)
{ 
    Arg arg[1];
    char text[32];

    if( wantLabels ){
	sprintf( text, "delay %4.2f ", delay );
	XtSetArg( arg[0], XtNlabel, text );
	XtSetValues( delay_label , arg, ONE );
    }
}

/* procs to set pause between frames, these are not currently set as
 * button callbacks but they may be used in a translation table 
 */
static void slower_proc(){ 
    if( fpause==0. ) fpause =.25; 
    else fpause = fpause*2.; 
    if( fpause >=10. ) fpause =9.99;

    set_delay_label(fpause);
}

void actionSlower(Widget w, XEvent *ev, String *p, Cardinal *np)
/*< slower >*/
{ slower_proc(); }

void faster_proc(void)
{ 
    fpause = fpause/2.; 
    if( fpause < .05 ) fpause =0.05;
    set_delay_label(fpause);
}

void actionFaster(Widget w, XEvent *ev, String *p, Cardinal *np)
/*< faster >*/
{ faster_proc(); }

/* the number that is updated by user keystrokes */
static int inputNumber=0;

void actionNumber(Widget w, XEvent *event, String *params, Cardinal *num_params)
/*< number >*/
{

    if( params && *params){
        int digits = strlen(*params);
        int value = atoi(*params);
        if( inputNumber == -1 ) inputNumber =0;
        while( digits > 0 ){
	    inputNumber *= 10;
	    digits--;
        }
        inputNumber += value;
    }
}

/* this is the global variable that is examained to see what the user typed */
int xt_next_num = -1;

static void numberReset(void)
{
    inputNumber = -1;
}

void actionNumReset(Widget w, XEvent *ev, String *p, Cardinal *np)
/*< numreset >*/
{ numberReset(); }

/* set next frame number to the users number and then do a next frame */
static void gotoFrame(void)
{
    xt_next_num = inputNumber ;
    if( next_on == YES ) didNEXT = YES;    return;
}

void actionGoto(Widget w, XEvent *ev, String *p, Cardinal *np)
/*< goto >*/
{ gotoFrame(); }


/* toggle run mode between Forward, Backward, Both-ways */

/* this routine sets the mode label */
void set_mode_label(const char *newlab)
{ 
    Arg arg[1];
    char text[32];

    if( wantLabels ){
	strncpy( text, newlab, 32 );
	XtSetArg( arg[0], XtNlabel, text );
	XtSetValues( mode_button , arg, ONE );
    }
}
/* this is the global variable that controls run mode */
/* the enumeration is defined in xtpen.h */
int xt_run_mode=XT_FORWARD;

void toggle_run_mode(
    Widget     w         /* widget */,
    XtPointer  x1        /* closure */,  /* data the application registered */
    XtPointer  x2  
    ){
    switch( xt_run_mode ) {
	case XT_FORWARD:
	    xt_run_mode = XT_BACKWARD;
	    set_mode_label("Backwards");
	    break;
	case XT_BACKWARD:
	    xt_run_mode = XT_BOTH_WAYS;
	    set_mode_label("Both Ways");
	    break;
	case XT_BOTH_WAYS:
	    xt_run_mode = XT_FORWARD;
	    set_mode_label("Forwards ");
	    break;
    }
	    
}

void actionRunMode(Widget  w, XEvent *ev, String *p, Cardinal *np)
/*< runmode >*/
{ toggle_run_mode(w,NULL,NULL); }

static void set_frame_label(int frame_num)
{ 
    Arg arg[1];
    char text[10];

    if( wantLabels ){
	sprintf( text, " %3d ", frame_num );
	XtSetArg( arg[0], XtNlabel, text );
	XtSetValues( frame_label , arg, ONE );
    }
}

void actionCoord(Widget w,XEvent *event,String *params,Cardinal *numparams)
/*< this action procedure prints the event location to the file "interact" >*/
{
    int x,y;

    if( interact[0] == '\0' ) return;

    if (first_time == YES){
        outfp = fopen (interact, "w");
        if (outfp == NULL) {
            ERR (FATAL, name, "Can't open interact output file %s!", interact);
        }
        first_time = NO;
    }

    if( event->type == ButtonPress ){
	x = XCORD (event->xbutton.x);
	y = YCORD (event->xbutton.y);
	devtovpxy ( x, y, &x, &y );
	if( boxy ){
	    doPointPopup( outfp, (float) x/RPERIN, (float) y/RPERIN );
	}else{
	    fprintf( outfp,"%f\t%f\n",(float) x/RPERIN,(float) y/RPERIN );
	    fflush( outfp );
	}

    }
}

void create_panel(Widget parent,bool want_buttons,bool want_labels)
{
    wantButtons = want_buttons;
    wantLabels = want_labels;

    if( wantButtons ){
  
        next_button = XtVaCreateManagedWidget("Next",commandWidgetClass,
					      parent, NULL);
        XtAddCallback(next_button, XtNcallback, next_proc, NULL);

        prev_button = XtVaCreateManagedWidget("Prev",commandWidgetClass,
					      parent,NULL);
        XtAddCallback(prev_button, XtNcallback, prev_proc, NULL);

        quit_button = XtVaCreateManagedWidget("Quit",commandWidgetClass,
					      parent,NULL);
        XtAddCallback(quit_button, XtNcallback, quit_proc, NULL);

        restart_button = XtVaCreateManagedWidget("Restart",commandWidgetClass,
						 parent,NULL);
        XtAddCallback(restart_button, XtNcallback, restart_proc, NULL);

        run_button = XtVaCreateManagedWidget("Run",commandWidgetClass,
					     parent,NULL);
        XtAddCallback(run_button, XtNcallback, run_proc, NULL);
 
        stop_button = XtVaCreateManagedWidget("Stop",commandWidgetClass,
					      parent,NULL);
        XtAddCallback(stop_button, XtNcallback, stop_proc, NULL);   

        stretch_button = XtVaCreateManagedWidget("Stretchy",commandWidgetClass,
						 parent,NULL);
        XtAddCallback(stretch_button, XtNcallback, stretch_proc, NULL);   
 	if( xt_stretchy )
	    XtSetValues( stretch_button, rigidArgs, XtNumber(rigidArgs) );

        mode_button = XtVaCreateManagedWidget("Forwards ",commandWidgetClass,
					      parent,NULL);
        XtAddCallback(mode_button, XtNcallback, toggle_run_mode, NULL);   

    }

    if( wantLabels ){
        frame_label = XtCreateManagedWidget("  0 ",labelWidgetClass,parent,
					    NULL,ZERO);

        delay_label = XtCreateManagedWidget("delay      ",
					    labelWidgetClass,parent, NULL,ZERO);

    }

    return;
}
 
void activate_buttons()
{
    if( wantButtons ){
	/* turn on the buttons, if required */
	XtSetSensitive( next_button, ( next_on == YES ) );
	XtSetSensitive( prev_button, ( prev_on == YES ) );
	XtSetSensitive( quit_button, ( quit_on == YES ) );
	XtSetSensitive( restart_button, ( restart_on == YES ) );
	XtSetSensitive( run_button,  ( run_on == YES ) );
	XtSetSensitive( stop_button, ( stop_on == YES ) );
	XtSetSensitive( stretch_button, ( size_on == YES ) );
    }
 
    return;
}

void inactivate_buttons()
{    
    /* turn all the buttons off */

    if( wantButtons ){
	XtSetSensitive( next_button, False );
	XtSetSensitive( prev_button, False );
	XtSetSensitive( quit_button, False );
	XtSetSensitive( restart_button, False );
	XtSetSensitive( run_button,  False );
	XtSetSensitive( stop_button, False );
	XtSetSensitive( stretch_button, False );
    }
   
    return;

}

int xt_pause(int doNEXT, int doPREV, int doREST, int doQUIT, int doRUN, int doSTOP, int doXSIZE)
{

    int retval=0;
    XtIntervalId id;


    next_on=doNEXT; prev_on=doPREV; quit_on=doQUIT; restart_on=doREST;
    run_on=doRUN;stop_on=doSTOP; size_on=doXSIZE;

    if( wantButtons ){
        activate_buttons();
    }

    if( doSTOP == YES ){ 
	unsigned long mpause; /* pause time in milliseconds */
	/* add a timeout event , add two millisecs to prevent zero */
        mpause = fpause*1000. + 2;
	didTIME = NO;
        id = XtAppAddTimeOut( pen_context,mpause, timeout_proc, NULL );

	/* we are in a run loop process events until the timer expires */
	/* only the stop and quit buttons will be active */

	while( !didTIME && !didSTOP && !didQUIT && !didNEXT && !didPREV ){
	    XtAppProcessEvent(pen_context,XtIMXEvent|XtIMTimer);
	}
	XtRemoveTimeOut(id);

	didTIME = NO;
	if( didSTOP ){
	    didSTOP=NO;	    return STOP;
        } else if (didQUIT) {
	    return  QUIT; 
	}else{
	    return TIME;
	}

    }else{
	/* if not in a run loop wait for a user event */
	while (!didNEXT && !didPREV && !didQUIT && 
	       !didFRAM1 && !didRUN && !didSTOP &&!didCHANGE) {
	    XtAppProcessEvent(pen_context,XtIMXEvent);
	}
  
	if (didNEXT) {
	    didNEXT = NO;	retval = NEXT;
        } else if (didPREV) {
	    didPREV = NO;    retval = PREV; 
        } else if (didRUN) {
	    didRUN = NO;     retval= RUN; 
        } else if (didSTOP) {
	    didSTOP = NO;    retval = STOP; 
        } else if (didQUIT) {
	    retval = QUIT; 
        } else if (didFRAM1) {
	    didFRAM1 = NO;	 retval = RESTART;  
        } else if (didCHANGE) {
	    didCHANGE = NO;	 retval = CHANGED;  
        }
	if( wantButtons ){
	    inactivate_buttons();
        }
	return retval;
    }

}

FileInfo       *inFiles;
int		num_files;

xtFrame *redraw_frame=0;
xtFrameList* frame_list; /* the list of frames to plot */

float fpause;

void xt_dovplot (int nn, FILE **inpltin, char **innames)
/*< do vplot >*/
{
    int i;
    int draw_file, draw_pos;
    int next_file, next_pos;
    int what_next;
    xtFrame* draw_frame;
    xtFrame* next_frame;
    int doNEXT,doPREV,doRESTART,doQUIT,doRUN,doSTOP,doSIZE;
    int running,pause;
    int going_forwards;
    Dimension width,height;

    if (nn == 0)
	return;

    xt_size_n_scale( &width,&height);


    num_files = nn;

    going_forwards = YES;

    /* allocate space for array inFiles[] */
    inFiles = (FileInfo *) sf_alloc (nn,sizeof (FileInfo));

    /* fill array inFiles */
    for (i = 0; i < nn; i++)
    {
	inFiles[i] .stream = inpltin[i];
	strcpy (inFiles[i] .name, innames[i]);
	inFiles[i] .fnew = YES;
    }

    /* ready to go! Generate an expose event */
    XFlush(pen_display);


    /* loop over plotting the files until the quit button is pressed */

    /* plot one frame at a time (forced by setting epause to -1 in
       xt_draw_frame, dovplot will return when a INT_PAUSE or
       INT_USER_PAUSE or end of file is found. 
       At this point we record the current frame and start a new frame 
       at the current position of the input stream */

    /* There are two modes for plotting frames:
     *  1) running mode  ( default mode if epause >= 0 )
     *  2) stopping mode ( default mode if epause < 0 ) stops after every frame
     */
    pause=1;
    if( !sf_getint( "pause",&epause) ) {
        pause=0;
        epause=app_data.pause;
    }

    running = ( epause >= 0 ); if( epause <0 ) epause=0;

    fpause = (float)epause; /* initial pause is the user's epause */
    if( pause == 0. ){
	if( !sf_getfloat( "fpause",&fpause) )  
	    if(fpause==0)fpause=.05;
    }

    set_delay_label(fpause);

    /* frame_list is a circular list of the frames so far */
    frame_list = xt_new_list( (xtFrameList*)0, (int*)0, 0 );

    /* start at the beginning of the first file */
    draw_file = 0; draw_pos = 0;
    draw_frame = xt_frame( frame_list, draw_file, draw_pos );

    /* Once this is set the Expose function actually does something */
    plotting_started=1;

    do{

	redraw_frame = draw_frame;

	xt_draw_frame( draw_frame ); /* draw the frame */

	set_frame_label( draw_frame->frame_num ); /*set the frame number label*/	/* keep it in sync every 10 frames */
	if( draw_frame->frame_num % 10 == 0 ){
            XSync(pen_display,0);
	}else{
            XFlush(pen_display);
	}


	if( draw_frame->end_pos == -1 ){  
	    /* end of file we must update next_file */
	    next_file = draw_frame->end_file + 1;
	    next_pos = 0;
	    if( next_file == nn ) {
		next_file = -1; /* at the end */

		/* check for case where we only have one frame */
		if( frame_list->num_frame == 1 ) running=NO;
	    }

	}else{ /* still in the same file */
	    next_file = draw_frame->end_file;
	    next_pos = draw_frame->end_pos;
	}


    	/* do any pausing and/or button reading */

	if( running ){ 

	   
	    /* run in a loop, only the stop & quit buttons are valid */
	    doNEXT = doPREV = doRESTART  = doRUN = doSIZE = NO;
	    doQUIT = doSTOP = YES; 

	    what_next = xt_pause(doNEXT,doPREV,doRESTART,
				 doQUIT,doRUN,doSTOP,doSIZE);
	     
            /* what to do if the timeout occurs */
	    if( what_next == TIME ){
	        switch( xt_run_mode ){
		    case XT_FORWARD:
			what_next = NEXT;
			break;
		    case XT_BACKWARD:
			what_next = PREV;
			break;
		    case XT_BOTH_WAYS:

	                if( draw_frame == xt_last_frame(frame_list) ){
			    going_forwards = NO;
			}

	                if( draw_frame == xt_first_frame(frame_list) ){
			    going_forwards = YES;
			}
 		
		        if( going_forwards ){
			    what_next = NEXT;
			}else{
			    what_next = PREV;
			}		    
		        break;
		}
	    }
	  
	}else{ 
	    /* pause for user interaction */
	    doRESTART = doQUIT = doRUN = doSIZE = YES;
            doPREV = (frame_list->num_frame != 1 );
	    doNEXT = !(frame_list->num_frame == 1 && next_file == -1 );
	    doSTOP = NO;
	    what_next =  xt_pause(doNEXT,doPREV,doRESTART,
				  doQUIT,doRUN,doSTOP,doSIZE);
	 
	}


	/* act on the result of the pause */
	switch(what_next){
	    case QUIT:
		draw_file = -1;
		break;

	    case STOP:
		running = NO;
		break;

	    case RUN:
		running = YES;
		break;
	    
	    case NEXT:

		/* see if the user typed a frame number */
		if( xt_next_num != -1  && 
		    (next_frame=xt_frame_num(frame_list,xt_next_num)) != 0  ){
		    draw_frame = next_frame;
		    draw_file = draw_frame->file_num;
		    draw_pos = draw_frame->file_position;
		    xt_next_num = -1;

		}else{

		    /* they didn't type anything or it was an invalid number */
		    if( (draw_frame == xt_last_frame(frame_list)) 
			&& next_file != -1 ){ 
			/* this is currently our last frame,,
			 * but we are not at the end of the input.
			 * So get next frame to plot */
			draw_file = next_file; draw_pos = next_pos;
			draw_frame = xt_frame( frame_list, draw_file, draw_pos ); 

		    }else{
			/* we already know what the next frame is */
			draw_frame = draw_frame->next;
			draw_file = draw_frame->file_num;
			draw_pos = draw_frame->file_position;
		    }
		}

		break;
				
	    case PREV:
		/* get previous frame to plot */
		draw_frame = draw_frame->prev;
		draw_file = draw_frame->file_num;
		draw_pos = draw_frame->file_position;
		break;

	    case CHANGED:
		/* something changed in the plotting parameters
		 * flush all cached images and replot the current frame */
		xt_clear_images(frame_list);
		xtredraw();
		break;

	    case RESTART:
		/* go back to the first frame */
		draw_frame = xt_first_frame(frame_list); 
		draw_file = draw_frame->file_num;
		draw_pos = draw_frame->file_position;
		break;
	}


    }while( draw_file != -1 );

    /* clean up */
    xt_clear_images(frame_list);
    remove_pixmap();
	
}

void xterase (int command)
/*< erase >*/
{
    /* ignore erase middle after a break */
    if( xt_after_break && command == ERASE_MIDDLE ){
	xt_after_break = 0;
	return;
    }

    xt_after_erase = 2;

    switch (command)
    {
	case ERASE_START:
	case ERASE_MIDDLE:
	    if( have_pixmap ){
		clear_pixmap();
	    }else{
		/* clear window */
		Dimension width, height;
		Position left, top;

		XtVaGetValues(pen_picture,XtNwidth,&width,XtNheight,&height,
			      XtNx,&left,XtNy,&top,NULL);

		XSetForeground(pen_display, pen_gc, map[0] );
		XFillRectangle(pen_display, pen_drawable, pen_gc,
			       0,0,width,height);
		XSetForeground(pen_display, pen_gc, color );
	    }
	    break;
    }
}

xtFrameList* xt_new_list(xtFrameList*  old_list, int *frames, int num_frame)
/*<
 * Maintain a circular list of frame info
 *
 * The info includes the file_name, file_position, frame_number, and image.
 *
 * Given a redraw_file,redraw_position pair we can then return a valid image.
 >*/
{
    xtFrameList* new_list;

    new_list = (xtFrameList*) sf_alloc(1,sizeof(xtFrameList));
   
    if( old_list == 0 ){
	new_list->start=0;
	new_list->end =0;
	new_list->parent=0;
	new_list->num_frame=0;
       
    }else{
	int i;
	for( i=0; i<num_frame; i++ ){
	    xt_add_frame( new_list, xt_frame_num( old_list, frames[i] ) );
	}
	new_list->parent = old_list;
    }
	
    return new_list;
}   
       

/* constructor for an new xtFrame object */       
xtFrame* xt_new_frame(int filenum, int filepos)
{
    xtFrame *newfr;

    newfr = (xtFrame*) sf_alloc(1,sizeof(xtFrame));
    newfr->file_num = filenum;
    newfr->file_position = filepos;
    newfr->end_pos = 0;
    newfr->total_len = 0;
    newfr->frame_num = 0;
    newfr->break_end = 0;
    newfr->has_image = 0; 
    newfr->image = 0;
    newfr->pixmap = 0;
    newfr->next =  0;
    newfr->prev =  0;

    return newfr;
}

void xt_add_frame(xtFrameList*  list, xtFrame* frame)
/*< add frame >*/
{

    xtFrame *end; 
    
    frame->frame_num = (list->num_frame)++;
    
    /* list->start is the first fame in the circulary linked list
       and list->start->prev is the last */
    
    if( list->start == 0 ){
	list->start = frame;
	frame->prev = frame;
	frame->next = frame;
    }else{
	end = list->start->prev;
	end->next = frame;
	frame->prev = end;
	list->start->prev = frame;
	frame->next = list->start;
    }    
}

xtFrame* xt_frame_num(xtFrameList* list, int num)
/*< frame number >*/
{
    xtFrame *curr, *start, *found;

    if( list==0 || list->start == 0 ) return (xtFrame*) 0;

    found = (xtFrame*)0;
    start = curr = list->start;
    
    do{
	if( num == curr->frame_num ) found=curr;
	curr = curr->next;
    }while( (found==(xtFrame*)0) && curr != start );
    
    return found;
}
 
xtFrame* xt_find_frame(xtFrameList *list, int filenum, int filepos)
{
    xtFrame *curr, *start;
    xtFrame *found;

    if( list==0 || list->start == 0 ) return (xtFrame*) 0;

    found = (xtFrame*)0;
    start = curr = list->start;

    do{
	if( filenum == curr->file_num && filepos == curr->file_position) found=curr;
	curr = curr->next;
    }while( (found==(xtFrame*)0 ) && curr != start );

    return found;

}


xtFrame*  xt_frame(xtFrameList*  list, int filenum, int filepos )
/*< frame >*/
{
    xtFrame *curr;
    curr = xt_find_frame( list, filenum, filepos );
    
    if( curr == 0 ) {
	curr = xt_new_frame( filenum, filepos );
	xt_add_frame( list, curr );
    }

    return ( curr );
}

void xt_store_image(xtFrame* fram )
{
    fram->has_image = 0;

/* get the image from the drawable and
 * copy it to either a Pixmap on the server or an XImage on the client
 */

    if( greedy_pixmaps ){

        if( fram->pixmap != (Pixmap)0 ) {
	    XFreePixmap(pen_display,fram->pixmap );
        }
        fram->pixmap = (Pixmap)0;

        fram->pixmap = MyCreatePixmap(pen_display,pen_drawable,
				      pen_width,pen_height,
				      visual_depth );

	if( fram->pixmap == (Pixmap)0 ){
	    /* allocation failed, give up on this method */
	    greedy_pixmaps = false;
	}else{
	    XSync( pen_display,0 );
	    XSetClipMask( pen_display, pen_gc, None );
	    XCopyArea(pen_display,pen_drawable,fram->pixmap,pen_gc,
		      0,0,pen_width,pen_height,0 ,0);
	    fram->has_image = 1;
	}

    }

    if( (!fram->has_image) && want_images ){

	/* We will save the image as an XImage in the client */

	/* saving the image in the client costs us a lot in
	   in data transmission so we try to make sure that the
           effort is worthwhile.

	   The image size = pen_width * pen_height * visual_depth
	   The size of the plot data is frame->total_len.
	   If the plot is a lot smaller than the image size we will
	   just replot it from scratch next time

	   This should probaly be related to an estimated connection speed
	   and the time taken to render the frame.
        */
	if( pen_width * pen_height * visual_depth < 50*fram->total_len ){
            XSync( pen_display,0 );
	    if( fram->image != 0 ){
		XDestroyImage( fram->image );
	    }
	    fram->image=0;
            fram->image = XGetImage( pen_display, pen_drawable, 0,0,
				     pen_width, pen_height, AllPlanes, ZPixmap );
            fram->has_image = 1;
	}
    }
    xt_save_colormap( fram );

}

/* put the image back into the drawable */
void xt_put_image(xtFrame*  fram)
{
    XImage* old_image;

    if( fram->pixmap != (Pixmap)0 ){
	XSetClipMask( pen_display, pen_gc, None );
        XCopyArea(pen_display,fram->pixmap,pen_drawable,pen_gc,
		  0,0,pen_width,pen_height,0 ,0);
    }else{
	old_image = fram->image;
	if( old_image != 0 ){
	    /* put the old image onto the drawable */
	    XPutImage( pen_display, pen_drawable, pen_gc, old_image, 0, 0, 0, 0,
		       old_image->width, old_image->height );
	}
    }
    xt_restore_colormap( fram );
}

void xt_clear_images(xtFrameList* list)
/*< clear images >*/
{
    xtFrame *curr;
    xtFrameList* currlist;
    
    currlist = list;
    if( list == 0 ) return;
    
    do{ 
	if( list->start != 0 ) {
    
	    curr = list->start;
	    do{
		/* free the image associated with this frame */
		if( curr->image != 0 )
		    XDestroyImage( curr->image ); curr->image=0;
		curr->image= 0;

		/* free the pixmap associated with this frame */
		if( curr->pixmap != (Pixmap)0 )
		    XFreePixmap(pen_display,curr->pixmap);
		curr->pixmap = 0;

		/* we don't have a stored image for this frame */
		curr->has_image = 0;

		curr = curr->next;
	    }while( curr != list->start );
	}
	currlist = list->parent;
    }while( currlist != 0 );
}


xtFrame* xt_first_frame(xtFrameList* list )
/*< first frame >*/
{
    return list->start;
}

xtFrame* xt_last_frame(xtFrameList* list)
/*< last frame >*/
{
    return list->start->prev;
}

static Cursor crossCursor = (Cursor)0;
static Cursor pointerCursor = (Cursor)0;
static Cursor currentCursor = (Cursor)-1;

static void cross_cursor(void)
{

    if( crossCursor == (Cursor)0 ){
   	crossCursor =  XCreateFontCursor(pen_display,XC_cross);
    }
    XDefineCursor( pen_display, pen_window, crossCursor );
    currentCursor = crossCursor;
    XFlush(pen_display);
}
	
static void pointer_cursor(void)
{
    if( pointerCursor == (Cursor)0 ){
   	pointerCursor =  XCreateFontCursor(pen_display,XC_left_ptr);
    }
    XDefineCursor( pen_display, pen_window, pointerCursor );
    currentCursor = pointerCursor;
    XFlush(pen_display);
}

int xtgetpoint (int *x, int *y)
{
    XEvent event;

    if( currentCursor != crossCursor ) cross_cursor();
    XSync (pen_display, True);
    XSelectInput (pen_display, pen_window, ButtonPressMask);
    XNextEvent(pen_display,&event);
    XSelectInput (pen_display, pen_window, 0);
    if( event.xbutton.button == Button1 ){
        *x = XCORD (event.xbutton.x);
        *y = YCORD (event.xbutton.y);
        return (0);
    }else{
	pointer_cursor();
	return(1);
    }
}

extern int      group_number;
int	skipit = NO;

int xtinteract (int what, FILE *controltty, char *string)
/*< interact >*/
{
 
    switch (what)
    {
	case INT_PAUSE:
	    if (skipit){
		skipit = NO;
		break;   /* skip the pause that occurs after a restart */
	    }
	case INT_USER_PAUSE:
	    skipit = NO;

	    /*
	     * Pauses caused by the "break" command too hard to handle!
	     */
	    if (group_number != 0){
		return DOVPLOT_CONT;
	    }

	    return DOVPLOT_EXIT; /* exit dovplot routine at end of frame */


	case INT_GET_STRING:
	    /* Not supported */
	    break;
	default:
	    break;
    }

    return DOVPLOT_CONT;

}

void xtaddText(const char *str)
{
  
    if ( text_region ) {
	int len ;
	int rc;
	char *p;
	char *currentString;
    
	Arg argList[20];
	Cardinal args;
    
	XawTextPosition start;
	XawTextBlock tblk, empty_tblk;
    
	len = strlen(str);

	if ( len <= 0 )
	    return;

	p = (char*)malloc( len+1 );
	strcpy( p, str );
	p[len] = '\0';
    
	tblk.firstPos = 0;
	tblk.length = len;
	tblk.ptr = p;
	tblk.format = FMT8BIT;
    
	/* empty text block for deletion */
	empty_tblk.firstPos = 0;
	empty_tblk.length = 0;
	empty_tblk.ptr = "";
	empty_tblk.format = FMT8BIT;

	args = 0;
	XtSetArg(argList[args], XtNstring, &currentString); args++;
	XtGetValues(text_region, argList, args);
    
	start = strlen( currentString );
    
	XawTextDisableRedisplay(text_region); 

	while ( start + len > TEXT_BUFFER_SIZE ) {
      
	    /* find first line and delete it */
      
	    char *eol = strchr( currentString, '\n');
      
	    if ( eol != 0 && eol < &currentString[TEXT_BUFFER_SIZE] ) {
		XawTextReplace( text_region, 0, eol-currentString+1, 
				&empty_tblk );
		start -= eol - currentString + 1;
	    }
	    else {
		XawTextReplace( text_region, 0, start, &empty_tblk );
		start = 0;
		break;
	    }
	}

	rc= XawTextReplace( text_region, start, start, &tblk );
	if (rc != XawEditDone) {
	    fprintf(stderr, "XawTextReplace(text, %d, %d, &textblock): %d\n",
		    (int) start, (int) start, rc);
	}
	XawTextSetInsertionPoint(text_region, start + len);
	XawTextEnableRedisplay(text_region);
	XawTextDisplay(text_region);

    }
}


void xtmessage (int command, const char *text)
/*< message >*/
{
    static int eatit = NO;

    switch (command)
    {
	case MESG_TEXT:
	    if (!eatit){
		if( text_region && plotting_started ){
		    xtaddText( text );
		}else{
		    printf ("%s", text);
		}
	    }
	    break;
	case MESG_HIGHLIGHT_ON:
	    eatit = YES;
	    break;
	case MESG_DONE:
	    eatit = NO;
	    if( !text_region && plotting_started ){
		fflush (stdout);
	    }
	    break;
    }
}

void opendev (int argc, char* argv[])
/*< open >*/
{
    int 		dwidth,dheight, mwidth,mheight;
    int             default_width,default_height;
    bool		want_buttons,want_labels,want_text;
    bool tellme_resolution;
    Visual		*default_visual;
    XVisualInfo	vinfo;
    int		pen_screen;
    Arg		args[20];
    int             cnt,depth;
    int iarg;

#if  XtSpecificationRelease > 4
    int		xtargc;
    int           myxargc;
#else
    Cardinal	xtargc;
    Cardinal        myxargc;
#endif
    char ** myxargv;
    const char *color;
    
    sf_parenv("SFPENOPTS");

    dev.message = xtmessage;
    dev.erase = xterase;
    dev.close = xtclose;

    dev.area = xtarea;
    dev.raster = xtraster;
    dev.point = xtpoint;
    dev.attributes = xtattributes;

    dev.reader = xt_dovplot;
//    dev.interact = xtinteract;
    dev.plot = xtplot;

    if (NULL == (interact = sf_getstring("interact"))) interact="";

    /*
     * save the command line arguments
     */
    
    myxargc = argc;
    myxargv = (char **) XtMalloc (myxargc * sizeof (char *));
    bcopy ((char *) argv, (char *) myxargv, myxargc * sizeof (char *));
    
    /* device capabilities */
    dev.txfont = DEFAULT_HARDCOPY_FONT;
    dev.txprec = DEFAULT_HARDCOPY_PREC;
    dev.brake = BREAK_IGNORE;
    endpause = false;
    dev.smart_clip = false;
    dev.cachepipe = true;

    /* initialize the first app (we will delete this one later ) */
    xtargc = argc;
    xtpen = XtVaAppInitialize(
	&pen_context,
	"XTpen",
	NULL, ZERO,
	&xtargc, argv,
	fallback_resources,
	NULL);
    
    /* get app data from the resource database */
    xt_app_data( xtpen );

    /* get info about the screen */
    pen_display = XtDisplay(xtpen);
    pen_screen = DefaultScreen(pen_display);
    screen_black = BlackPixel( pen_display,pen_screen);
    screen_white = WhitePixel( pen_display,pen_screen);

    if(!sf_getint ("depth",&depth)) depth=app_data.vis_depth;
    /* Choose a visual */
    choose_visual( pen_display, pen_screen, depth, &vinfo );
		
    pen_visual = vinfo.visual;
    screen_depth = visual_depth;
/*    ERR( COMMENT,name,"visual depth = %d",visual_depth);*/


    /* figure out the screen size */
    dwidth = DisplayWidth(pen_display,pen_screen);
    dheight = DisplayHeight(pen_display,pen_screen);

    /* Height and Width of the screen in Millimeters */
    mwidth = DisplayWidthMM( pen_display,pen_screen);
    mheight = DisplayHeightMM( pen_display,pen_screen);
    dev.pixels_per_inch = ((float)dwidth/(float)mwidth)*25.4 ;
    dev.aspect_ratio = (((float)dheight/(float)mheight)*25.4 )/dev.pixels_per_inch ;

    /*
     * Since X rounds values to integer millimeters, the aspect ratio
     * has some error in it. Push aspect ratio to nearest round percent.
     */
    dev.aspect_ratio = ((int)((dev.aspect_ratio * 100.) + .5)) / 100.;

    if (!sf_getbool("x_screen_info", &tellme_resolution)) 
	/* output the default screen information */
	tellme_resolution = false;
    if (tellme_resolution) {
	ERR(COMMENT, name,
	    "display width=%d, height=%d (pixels);  width=%d, height=%d (mm)",
	    dwidth, dheight, mwidth, mheight);
	ERR(COMMENT, name, "pixels_per_inch=%f,  aspect_ratio=%f",
	    dev.pixels_per_inch, dev.aspect_ratio);
    }

    sf_getfloat("aspect",&dev.aspect_ratio);
    /* aspect ratio */
    sf_getfloat("ppi",&dev.pixels_per_inch);
    /* pixels per inch */
 
    if ( !sf_getint("numcol",&num_col_req) ) num_col_req = app_data.num_col;
    /* number of colors (take a default from the resource database) */

    /* colormap support */
    if( !mono ) x_num_col = DisplayCells(pen_display, pen_screen);
    if( x_num_col <= 2 ) mono=true;
    if( x_num_col < num_col_req ) num_col_req = x_num_col;
    /* limit to 254 to get around tektronix bug */
    if( num_col_req > 254 ) num_col_req=254;

    /* if the user specifies the number of colors dont accept anything less */
    if( num_col_req != 0 ){
	num_col_min = num_col_max = num_col_req;
    }

    /* initialze argument list counter to zero */
    cnt=0;

    /* For the default visual set the initial colormap to be 
       the default colormap, otherwise aloocate your own */
    default_visual = DefaultVisual( pen_display, pen_screen);
    if( pen_visual->visualid  == default_visual->visualid ){
	pen_colormap = DefaultColormap(pen_display, pen_screen);
    }else{
	pen_colormap = XCreateColormap( pen_display, 
					DefaultRootWindow(pen_display) , pen_visual, AllocNone );
	XtSetArg (args[cnt], XtNvisual, pen_visual); ++cnt;
	XtSetArg (args[cnt], XtNdepth, visual_depth ); ++cnt;
    }

    if(  mono  ){
	dev.num_col = x_num_col;
    }else{
	if( xtruecol ){
	    dev.num_col = num_col_max;
        }else{
	    /* NB may change pen_colormap ! */
            dev.num_col = xtnumcol( num_col_min, num_col_max );
	}
    }

    XtSetArg (args[cnt], XtNcolormap,pen_colormap ); ++cnt;
    /* end of colormap code */


    /* now destroy the widget we created and make a new one with our 
       chosen visual and colormap etc. */
    XtDestroyWidget( xtpen );
    
    /*
     * Now create the real toplevel widget.
     */
        
/*    XtSetArg (args[cnt], XtNargv, myxargv); ++cnt;*/
/*    XtSetArg (args[cnt], XtNargc, myxargc); ++cnt;*/
    iarg=0;
    while(iarg < myxargc){
	if(myxargv[iarg][0]=='-' && iarg < myxargc-1){
	    XtSetArg(args[cnt],&(myxargv[iarg][1]),myxargv[iarg+1]);++cnt;
	    iarg++;
	}
	iarg++;
    }
    



    xtpen = XtAppCreateShell ( "xtpen", "XTpen",
			       applicationShellWidgetClass,
			       pen_display, args, cnt);


    /* default size for the picture is 1/2 of screen height and width
     * this can be overridden by resources or the command line geometry spec 
     */
    default_width = dwidth*.5;
    default_height = dheight*.5;
  
    vpane = XtVaCreateManagedWidget("vpane",panedWidgetClass,
				    xtpen, 
				    XtNwidth, default_width,
				    NULL);

    if( !sf_getbool("buttons",&want_buttons) ) 
	want_buttons = (bool) app_data.buttons;
    /* if y, display a panel of buttons on top of the plot */
    if( !sf_getbool("labels",&want_labels) ) 
	want_labels = (bool) app_data.labels;
    /* if y, display frame number and inter-frame delay at the top of plot */
    if( !sf_getbool("want_text",&want_text) ) 
	want_text = (bool) app_data.textpane;
    /* if y, display a message window */
    if( !sf_getbool("stretchy",&xt_stretchy) ) 
	xt_stretchy= (bool) app_data.stretchy;
    /* if y, use the stretchy mode and fill the window */

    if( !sf_getbool("boxy",&boxy) ) boxy = false;
    /* output coordinates and labels suitable for sfbox */

    if( !sf_getbool("see_progress",&see_progress) ) see_progress=false;
    /* show progress of each frame, slow */

    if( !sf_getbool("images",&want_images) ) 
	want_images = (bool) app_data.images;
    /* copy the image created by plotting each frame and save it in
       the client program (xtpen). This will increase memory usage in
       the machine that runs the pen command. If you have a fast
       connection to your X-server it will make redisplay of frames
       faster. If you have a slow connection, it may make replotting
       slower. */

    if( !sf_getbool("pixmaps",&greedy_pixmaps)) 
	greedy_pixmaps = (bool) app_data.pixmaps;
    /* Copy the image created by plotting each frame and save it in
       the X-server. This will increase memory usage of the machine that
       displays the window! Redisplay of frames will be very fast and
       the network traffic is very low so this is a suitable option for
       slow connections.  If your X-server is a workstation with plenty
       of memory and swap space then this option should be very useful.
       If your X-server has limited memory, this option may have
       undesirable effects on the response of your terminal. */

    if( want_buttons || want_labels ){
        control_panel = XtCreateManagedWidget("control_panel",boxWidgetClass,
					      vpane, NULL, ZERO);

        create_panel(control_panel,want_buttons,want_labels);
    }

    /* if they want to add a message to the window do it here */

    if( visual_depth != 8 ) want_text=false; /* BUG for 24 bit visuals */
    if( want_text ) {
        if(NULL == (message_text = sf_getstring("message"))) message_text="";

        text_region = XtVaCreateManagedWidget("text",
					      asciiTextWidgetClass,vpane,
					      XtNstring, (XtArgVal)message_text,
					      XtNlength, (XtArgVal)TEXT_BUFFER_SIZE, 
					      XtNeditType, XawtextEdit,
					      NULL );
    }

    pen_picture = XtVaCreateManagedWidget("pen_picture",
					  widgetClass, vpane,
					  XtNwidth, default_width,
					  XtNheight, default_height,
					  XtNcolormap, pen_colormap,
					  NULL);



    XtAppAddActions(pen_context,window_actions, XtNumber(window_actions));

    XtAugmentTranslations(pen_picture,XtParseTranslationTable(trans));

    XtRealizeWidget(xtpen);

    inactivate_buttons();

/* moved into xt_dovplot (should it be xt_reset ? )*/
/*    xt_size_n_scale( &width, &height);*/
  
    /* Pixmap is initially undefined, the drawable is initially the window */
    pen_pixmap = (Pixmap)0;
    have_pixmap = NO;
    pen_window = XtWindow(pen_picture);
    pen_drawable = pen_window;

/*
 * create and initialize gc
 */
    pen_gc = XCreateGC(pen_display,pen_window,0,NULL);
    XSetLineAttributes(pen_display,pen_gc, 0, LineSolid, CapButt, JoinRound);
    XSetFillRule(pen_display, pen_gc, EvenOddRule);
    XSetClipMask(pen_display, pen_gc, None );
 
    /* set standard colors */
    if( mono ){
	if( xmono ){
	    map[1] = screen_white;
	    map[0] = screen_black;
	}else{
	    map[7] = screen_white;
	    map[1] = screen_white;
	    map[0] = screen_black;
        }
    }else{
	if (NULL == (color = sf_getstring("bgcolor"))) color="black";
	/* background color */
	if (color[0] == 'w' || color[0] == 'l') {
	    xtattributes (SET_COLOR_TABLE, 0, MAX_GUN, MAX_GUN, MAX_GUN);
	    xtattributes (SET_COLOR_TABLE, 1, MAX_GUN, MAX_GUN, 0);
	    xtattributes (SET_COLOR_TABLE, 2, 0, MAX_GUN, MAX_GUN);
	    xtattributes (SET_COLOR_TABLE, 3, 0, MAX_GUN, 0);
	    xtattributes (SET_COLOR_TABLE, 4, MAX_GUN, 0, MAX_GUN);
	    xtattributes (SET_COLOR_TABLE, 5, MAX_GUN, 0, 0);
	    xtattributes (SET_COLOR_TABLE, 6, 0, 0, MAX_GUN);
	    xtattributes (SET_COLOR_TABLE, 7, 0, 0, 0);
	} else { 
	    xtattributes (SET_COLOR_TABLE, 0, 0, 0, 0);
	    xtattributes (SET_COLOR_TABLE, 1, 0, 0, MAX_GUN);
	    xtattributes (SET_COLOR_TABLE, 2, MAX_GUN, 0, 0);
	    xtattributes (SET_COLOR_TABLE, 3, MAX_GUN, 0, MAX_GUN);
	    xtattributes (SET_COLOR_TABLE, 4, 0, MAX_GUN, 0);
	    xtattributes (SET_COLOR_TABLE, 5, 0, MAX_GUN, MAX_GUN);
	    xtattributes (SET_COLOR_TABLE, 6, MAX_GUN, MAX_GUN, 0);
	    xtattributes (SET_COLOR_TABLE, 7, MAX_GUN, MAX_GUN, MAX_GUN);
	}
    }

    /* tell which entry in the colormap is the background color */
/*    XSetWindowBackground( pen_display, pen_window, map[0]);*/
    xtattributes (SET_COLOR, 7, 0, 0, 0);;

    /* xtpen controls the file stream itself */
    dev.reader = xt_dovplot;

    /* xtpen controls messages */
    message = dev.message;

 
}

extern float user_xscale;
extern float user_yscale;

static void xt_size_n_scale(Dimension *width, Dimension *height)
{  

    Position left, top;
    float win_ratio; 
    float new_stretch, new_xstretch,new_ystretch;
    static float old_xstretch = 1.;
    static float old_ystretch = 1.;

    /* now get the size of the actual picture and set dev.xmax etc. */

    XtVaGetValues(pen_picture,XtNwidth,width,XtNheight,height,
		  XtNx,&left,XtNy,&top,NULL);
  
    /* set VPLOT plotting window */
    dev.xmin = left;
    dev.xmax = dev.xmin + *width - 1;
    dev.ymin = top;
    dev.ymax = dev.ymin + *height - 1;  

    /* attempt to create pixmap  to speed redraws */
    create_pixmap( *width, *height );

    if( xt_stretchy ){ 
	/* If the pen is in stretch mode we should scale the plots
	 * axis so that the plot fills the frame */
	
	win_ratio = (float)(*height) / (float)(*width);
	
	new_stretch = win_ratio/VP_SCREEN_RATIO;

	if( new_stretch > 1. ){
	    new_ystretch = new_stretch; new_xstretch =1.;
	}else{
	    new_xstretch = 1./new_stretch; new_ystretch =1.;
	}

    }else{
	/* go back to a rigid aspect ratio */
	new_xstretch=1.;
	new_ystretch=1.;
    }


    user_yscale = user_yscale * new_ystretch / old_ystretch;
    user_xscale = user_xscale * new_xstretch / old_xstretch;
	
    old_xstretch = new_xstretch; old_ystretch = new_ystretch;
  
}

static void choose_visual(Display* display, int screen, int max_depth, 
			  XVisualInfo* vinfo)
{

    /* try and figure out a good visual */
  
    /* try for 24bit true color */
    if(  max_depth >=24 &&
	 XMatchVisualInfo(display,screen,24,TrueColor,vinfo) ){
	xmono=NO;
	xtruecol=YES;
	visual_depth = 24;
     
	/* try for 8 bit pseudo color */
    }else if( max_depth >=8 &&
	      XMatchVisualInfo(display,screen,8,PseudoColor,vinfo) ){
	xmono=NO;
	xtruecol=NO;
	visual_depth = 8;
	x_num_col = 256;
     
	/* try for 16bit true color */
    }else if( max_depth >=16 &&
	      XMatchVisualInfo(display,screen,16,TrueColor,vinfo) ){
	xmono=NO;
	xtruecol=YES;
	visual_depth = 16;
     
	/* try for 8 true color */
    }else if( max_depth >=8 &&
	      XMatchVisualInfo(display,screen,8,TrueColor,vinfo) ){
	xmono=NO;
	xtruecol=YES;
	visual_depth = 8;
     
	/* try for 4 bit pseudo color */
    }else if( max_depth >=4 &&
	      XMatchVisualInfo(display,screen,4,PseudoColor,vinfo) ){
	xmono=NO;
	num_col_min = 8;
	x_num_col = 8;
	visual_depth = 4;
	xtruecol=NO;
     
	/* if we can't get any of those try for a monochrome visual */
	/* try for 1 bit static grey */
    }else if( max_depth >=1 &&
	      XMatchVisualInfo(display,screen,1,StaticGray,vinfo) ){
	xmono = YES;
	visual_depth = 1;
	xtruecol=NO;
	mono = true;
	x_num_col = 2;
    }else{
	ERR(FATAL,name,"Could not obtain a suitable visual");
    }   
}

int             in_repaint = NO;

int xt_endframe = NO;
 
/* draw one frame from the file stream, epause is set to -1 to force
 * a USER_INT command at the end of the frame.
 * this will cause dovplot to exit at the end of the frame.
 */
void xt_draw_file(int start_file, 
                  /* index of start file in the array of input files */
		  long start_pos, 
                  /* offset in the file to start reading the plot commands */
		  int *end_file,
                  /* index of the file we finished reading this frame */
		  int *end_pos , 
                  /* position at which we eneded reading this frame */  
		  long *total_len 
		  /* the total length of vplot data for this frame */
    )
{

    int save_pause;
    int file_num;
    int file_pos;

    file_num = start_file;
    file_pos = start_pos;

    save_pause = epause;
    epause = -1;
    xt_endframe=0;
    reset_parameters ();


    *total_len =0;

    do {

        pltin = inFiles[file_num] .stream;
        strcpy (pltname, inFiles[file_num].name);

	if( start_pos==0 && inFiles[file_num].fnew ){
            inFiles[file_num].fnew = NO;
    	    dovplot ();
	}else{
            inFiles[file_num].fnew = NO;

	    if( start_pos < -1 ){
		/* start_pos will be < -1 for a pipe input */
		if( *end_pos < -1 ){
		    /* this is a replot request */
		    ERR (WARN, name, "sorry, unable to draw frame.");
		    xt_endframe = YES ;
		}else{
		    /* just do it! */
		    dovplot ();
		}
	    }else{
	        /* this is for a file input */
                if ( !fseek (pltin, file_pos, 0)) {
    	            dovplot ();
                }else{
                    ERR (WARN, name, "sorry, unable to draw frame.");
	            xt_endframe = YES ;
                }
            }
	}

	/* find the ending position */
        *end_pos = ftell( inFiles[file_num].stream );

	if( *end_pos == -1 ){ 
	    /* when reading from a pipe */
	    *end_pos = start_pos - 2;
	    *total_len += 99999; /* fake it to make it save images */
        }else{
	    /* accumalate a length measure */
	    *total_len += *end_pos - file_pos;
	}
	

        if( feof(inFiles[file_num].stream) ){
	    *end_pos = -1;
	}

	/* now see if this is the end of a frame. If it isn't it means that
         * dovplot returned but the frame has not ended 
         * so we should keep plotting */
	if( ! xt_endframe ){
	    if( *end_pos == -1 ){
	        /* go on to the next file */
	        file_num++; file_pos = 0;
	    }else{
		/* keep plotting from this position */
		file_pos = *end_pos;
	    }

	    /* check to see if we have reached the end of all the files */
	    if( file_num == num_files ) {
		file_num--; xt_endframe = YES ;
	    }
	}

        skipit = YES; /* force skiping of pause at start of next dovplot */

    }while( ! xt_endframe );

    *end_file = file_num;


    epause = save_pause;

}

void xt_draw_frame(xtFrame* frame )
/*< draw frame >*/
{
    int end_file;
    int end_pos;
    long total_len;

    if( frame->has_image ){
	xt_put_image( frame );
	end_pos = frame->end_pos;

    }else{
	if( frame->prev != 0 ){
	    if( frame->prev->break_end ) xt_after_break = 1;
	}
	end_pos = frame->end_pos;
	xt_draw_file(frame->file_num, frame->file_position, &end_file, 
		     &end_pos, &total_len);
	frame->total_len = total_len;
	frame->end_file = end_file;
	frame->end_pos = end_pos;
	if( xt_after_break ) frame->break_end = 1;
	if( want_images || greedy_pixmaps ) xt_store_image( frame );

    }

    xt_flush_display();

}

void xtredraw (void)
/*< paint a file, get the window size, allocate a pixmap and call xt_draw_file 
  save the current file info because this may not be the one we are redrawing
  dovplot can be confused if we change the file pointers in midstream >*/
{
    static FILE     *save_pltin;
    static char     save_pltname[120];
    Dimension width, height;
    xtFrame* replot_frame;
 

    xt_size_n_scale(&width,&height);

    inactivate_buttons();
  
    /* save the current file info */
    /* we might be about to call dovplot to repaint the previous file */
    save_pltin = pltin;
    strcpy (save_pltname, pltname);

    /* to handle breaks correctly we need to backup redplot_frame
     * to the last frame that started with an erase (or the first frame)*/
    replot_frame = redraw_frame ;
    do{
	if( replot_frame->prev == 0 ) break;
	if( !replot_frame->prev->break_end ) break;
	replot_frame = replot_frame->prev;
    }while( replot_frame->frame_num != 0 );
	

    do{	
        xt_draw_frame( replot_frame );
	replot_frame = replot_frame->next;
    }while( replot_frame != redraw_frame->next );

    /* restore file info */
    strcpy (pltname, save_pltname);
    pltin = save_pltin;

    activate_buttons();

}


/*
 * repaint plot by copying in-memory pixmap to canvas
 */

static void
PenRepaint (Widget w, XEvent *event, String *params, Cardinal *num_params)
{
    Dimension width, height;
    Position left, top;


    XtVaGetValues(pen_picture,XtNwidth,&width,XtNheight,&height,
		  XtNx,&left,XtNy,&top,NULL);


    if( have_pixmap && (width == pen_width )&& (height == pen_height)){
	/* lazy man's way! copy the whole pixmap */
        XSetClipMask( pen_display, pen_gc, None );
        XCopyArea(pen_display, pen_pixmap, pen_window, pen_gc, 
		  0, 0, width, height, 0, 0);
        XFlush(pen_display);
    }else{
        if( in_repaint ) return;
        in_repaint = YES;
    	xtredraw() ;
        in_repaint = NO;
    }

}

void xt_flush_display(void)
/*< flush display >*/
{
    if( have_pixmap ){
        /* copy the pixmap to the screen */
        XSetClipMask( pen_display, pen_gc, None );
        if( pen_pixmap != (Pixmap)0) XCopyArea(pen_display,pen_pixmap,
					       pen_window,pen_gc,
					       0,0,pen_width,pen_height,0 ,0);
    }
    XFlush(pen_display);
}

int pen_width = -1;
int pen_height = -1;
int have_pixmap;

/* draw a filled rectangle on the pixmap to clear it */
static void clear_pixmap(void)
{
    XSetForeground(pen_display, pen_gc, map[0] );
    XFillRectangle(pen_display, pen_pixmap, pen_gc,
		   0,0,pen_width,pen_height);
    XSetForeground(pen_display, pen_gc, color );
}

static void remove_pixmap(void)
{
    if( pen_pixmap != (Pixmap)0 ) {
        XFreePixmap(pen_display,pen_pixmap);
	pen_pixmap = (Pixmap)0;
    }
}

/* this is the pixmap we draw onto, it should only be created when we
   first create a plot or when it is resized */

static void  create_pixmap(int width, int height )
{

    /* Make sure we have a pixmap of the correct size */

    if( pen_width == width && pen_height== height ) { return; }

    remove_pixmap();

    xt_clear_images(frame_list);

    if(  see_progress ){ /* we don't want to draw to an offscreen pixmap */
    	pen_height = height; pen_width = width;
	have_pixmap = NO;
	pen_drawable = pen_window;
	return;
    }

    pen_pixmap = MyCreatePixmap(pen_display,pen_window,
				width,height,
				visual_depth );

    /* check if pixmap creation succeeded */
    if( pen_pixmap != (Pixmap)0 ) {
    	pen_height = height; pen_width = width;
    	clear_pixmap();
    	have_pixmap = YES;
    	pen_drawable = pen_pixmap;
    }else{
    	pen_height = height; pen_width = width;
	have_pixmap = NO;
	pen_drawable = pen_window;
    }

}

/* pixmap allocation with an error handler to detect failed allocation */

static int attachfailed; 
 
static int MyHandler(Display *dpy, XErrorEvent *errorevent)
{
    attachfailed = 1;
    return(0);
}

Pixmap MyCreatePixmap(Display *display, Drawable drawable,
		      unsigned int width,
		      unsigned int height, unsigned int depth)
{
    Pixmap pix;

    XErrorHandler handler;
    XSync(display,0);
    attachfailed  = 0;

    handler = XSetErrorHandler(MyHandler);      /* Start critical time */
    pix = XCreatePixmap(display,drawable, width,height, depth );
    XSync(display,0);
    (void) XSetErrorHandler(handler);           /* Critical time over */

    if( attachfailed ) return (Pixmap)0; 
    else return pix;

}

void xtplot (int x, int y, int draw)
/*< plot >*/
{
    static int oldx = 0, oldy = 0;
    
    if (draw)
    {
	XDrawLine (pen_display, pen_drawable, pen_gc, XCORD (oldx), YCORD (oldy), XCORD (x), YCORD (y));
    }
    oldx = x;
    oldy = y;
}

void xtpoint (int x, int y)
/*< point >*/
{
    XDrawPoint (pen_display, pen_drawable, pen_gc, XCORD (x), YCORD (y));
}

#define	MASK0 0200

void xtraster (int count, int out_of, int xpos, int ypos, int length, int orient, 
	       unsigned char **raster, int dummy1, int dummy2)
/*< raster >*/
{
    int             i, j, k;
    XImage	       *pen_image;
    unsigned char  *rp, *re, *rp3;
    unsigned char  *rp2;
    static unsigned char *raster2;
    static unsigned char *raster3;
    static int      xloc, yloc;
    static int      offset, incr1, incr2;
    static int      width, height, widthpad;
    int		format,xoffset,bitmap_pad,bytes_per_line;
    unsigned int	depth,xwidth=0,xheight=0;
    unsigned char     *tip;
    int              do32;
    unsigned long   xcol;
 

/*
 * The following comments were written to go along with the
 * SunView vplot filter "sunpen". I've left them here in the hope
 * that even if they are not completely relevant they will help
 * others understand the code.
 * - S. Cole, 3 Aug 91
 *
 * Comments to help understand the code for different orientations:
 * XCORD(xpos), YCORD(ypos) are the coordinates of the first point in
 * the first line of raster passed to this routine. For the different
 * orientations:         this point is this corner of the raster:
 *      0			upper left
 *      1			upper right
 *      2			lower right
 *      3			lower left
 * The Sunview routines pw_rop and pr_rop want the coordinates of the
 * upper left corner. Keeping in mind the order in which vplot fills
 * the raster (for orient 0, this is TV-style: from left to right and
 * top to bottom; orient 1 is this rotated 90 degrees clockwise,
 * orient 2 180 degrees, orient 3 270 degrees) the computations of
 * xloc, yloc, width, and height should be clear.
 * The mapping of pixels in the input raster to the new rotated raster
 * is done in core, because there is no easy way to do it in Sunview
 * (that I could find). offset is the location where the first pixel
 * of the first row belongs on output (relative to the start of the
 * array), incr1 is the change in output location as you move from one
 * pixel to the next in the input raster line, and incr2 is the change
 * in output address as you move from one line to the next on input.
 * Together, offset incr1 and incr2 provide a simple way to determine
 * the output address from the input address for any pixel.
 * - S. Cole
 */
    if (count == 0)
    {
	if (orient == 0)
	{
	    xloc = XCORD (xpos);
	    yloc = YCORD (ypos);
	    width = length;
	    height = out_of;
	    if (xmono)
		widthpad = ((width + 15) / 16) * 16;
	    else
		widthpad = ((width + 1) / 2) * 2;
	    offset = 0;
	    incr1 = 1;
	    incr2 = widthpad;
	}
	else
	    if (orient == 1)
	    {
		xloc = XCORD (xpos) - out_of + 1;
		yloc = YCORD (ypos);
		width = out_of;
		height = length;
		if (xmono)
		    widthpad = ((width + 15) / 16) * 16;
		else
		    widthpad = ((width + 1) / 2) * 2;
		offset = width - 1;
		incr1 = widthpad;
		incr2 = -1;
	    }
	    else
		if (orient == 2)
		{
		    xloc = XCORD (xpos) - length + 1;
		    yloc = YCORD (ypos) - out_of + 1;
		    width = length;
		    height = out_of;
		    if (xmono)
			widthpad = ((width + 15) / 16) * 16;
		    else
			widthpad = ((width + 1) / 2) * 2;
		    offset = widthpad * height - 1 - widthpad + width;
		    incr1 = -1;
		    incr2 = -widthpad;
		}
		else
		    if (orient == 3)
		    {
			xloc = XCORD (xpos);
			yloc = YCORD (ypos) - length + 1;
			width = out_of;
			height = length;
			if (xmono)
			    widthpad = ((width + 15) / 16) * 16;
			else
			    widthpad = ((width + 1) / 2) * 2;
			offset = widthpad * (height - 1);
			incr1 = -widthpad;
			incr2 = 1;
		    }
	/* raster2 = (unsigned char *) calloc ((widthpad / 2) * 	
	   visual_depth * height / 8, 2); */
	raster3 = (unsigned char *) calloc ((widthpad / 2) * height, 2);

    }

    if (xmono)		/* Just store in correct order AS 2/2/89 */
    {
	rp3 = raster3 + offset + count * incr2;
	re = raster[0] + length;
	for (rp = raster[0]; rp < re; rp++, rp3 += incr1)
	{
	    /* if( *rp )
	     *rp3 = screen_black;
	     else
	     *rp3 = screen_white; */
	    *rp3 = !!(*rp); 
	}
    }
    else if( visual_depth == 8 ) 
    {
        /* convert color rasters to 8 bit pixels in raster3 */
	rp3 = raster3 + offset + count * incr2;
	re = raster[0] + length;
	for (rp = raster[0]; rp < re; rp++, rp3 += incr1)
	{
	    *rp3 = map[*rp];
        }

    } 
    else
    {
        /* just keep an interpolated raster of the color indicies in raster3 */
	rp3 = raster3 + offset + count * incr2;
	re = raster[0] + length;
	for (rp = raster[0]; rp < re; rp++, rp3 += incr1)
	{
	    *rp3 = *rp;
	}
    }

    if (count != out_of - 1)
	return;

    if (xmono)
    {
	/* Convert raster3 to a bitmap AS 2/2/89 */
	raster2 = (unsigned char *) calloc ((widthpad / 2) * height / 8, 2);
	rp2 = raster2;
	rp3 = raster3;
	for (i = 0; i < height; i++)
	{
	    for (j = 0; j < widthpad / 8; j++, rp2++)
	    {
		for (k = 7; k >= 0; k--, rp3++)
		{
		    *rp2 |= (*rp3 << k);
		}
	    }
	}
    	depth = 1;
   	format = XYBitmap;
   	xoffset = 0;
    	bitmap_pad = 1;
    	bytes_per_line = widthpad/8;
    	pen_image = XCreateImage(pen_display,pen_visual,depth,format,xoffset,
				 (char *) raster2,width,height,bitmap_pad,bytes_per_line);

    } else if( visual_depth == 8 ){
	/* its easy to make 8 bit images, fast */
    	depth = 8;
   	format = ZPixmap;
   	xoffset = 0;
	xwidth = width;
	xheight = height;
    	bitmap_pad = 8;
    	bytes_per_line = widthpad;
    	pen_image = XCreateImage(pen_display,pen_visual,depth,format,xoffset,
				 (char *) raster3,xwidth,xheight,bitmap_pad,bytes_per_line);
	if (!pen_image) ERR(FATAL, name, "couldn't create XImage!");


    }else if( visual_depth == 16 ||
	      visual_depth==12 ||
	      visual_depth == 15 ){
	/* convert raster3 to a 12/15/16 bit image in raster2 */
	raster2 = (unsigned char *) calloc((size_t) (2*width*height),1);
	depth = visual_depth;
	format = ZPixmap;
	xoffset = 0;
	xwidth = width;
	xheight = height;
	bitmap_pad = 16;
	bytes_per_line = 0;
	pen_image = XCreateImage(pen_display,pen_visual,depth,format,xoffset,
				 (char *)raster2,xwidth,xheight,bitmap_pad,
				 bytes_per_line);
	if (!pen_image) ERR(FATAL, name, "couldn't create XImage!");
	 
	if (visual_depth == 12 && pen_image->bits_per_pixel != 16) {
	    ERR(FATAL, name,
		"Can't handle visual_depth==12 && bits_per_pixel != 16 ");
	}

	rp3 = raster3;
	
	if (pen_image->byte_order == MSBFirst) {
	   
	    for (i=0, rp2=raster2; i<height; i++) {
		rp3 = raster3 + i*widthpad;
		for (j=0, tip=rp2; j<width; j++, rp3++) {
		    xcol = map[*rp3];
		    *tip++ = xcol >>8 & 0xff;
		    *tip++ = xcol & 0xff;
		}
		rp2 += pen_image->bytes_per_line;
	    }
	   
	} else {  /* LSBFirst */
	   
	    for (i=0, rp2=raster2; i<height; i++) {
		rp3 = raster3 + i*widthpad;
		for (j=0, tip=rp2; j<width; j++, rp3++) {
		    xcol = map[*rp3];
		    *tip++ = xcol & 0xff;
		    *tip++ = xcol >>8 & 0xff;
		}
		rp2 += pen_image->bytes_per_line;
	    }
	}
       

    } else if( visual_depth == 24 || visual_depth == 32 ) {
        /* convert raster3 to a 24/32 bit image in raster2 */
	raster2 = (unsigned char *) calloc((size_t) (4*width*height),1);
    	depth = visual_depth;
   	format = ZPixmap;
   	xoffset = 0;
	xwidth = width;
	xheight = height;
    	bitmap_pad = 32;
    	bytes_per_line = 0;
    	pen_image = XCreateImage(pen_display,pen_visual,depth,format,xoffset,
				 (char *)raster2,xwidth,xheight,bitmap_pad,bytes_per_line);
        if (!pen_image) ERR(FATAL, name, "couldn't create XImage!");
 
        do32 = (pen_image->bits_per_pixel == 32);
 
 
        if (pen_image->byte_order == MSBFirst) {
	  
	    for (i=0, rp2=raster2; i<height; i++) {
		rp3 = raster3 + i*widthpad;
		for (j=0, tip=rp2; j<width; j++, rp3++) {
		    xcol = map[*rp3];
		    if (do32) *tip++ = 0;
		    *tip++ = (xcol>>16) & 0xff;
		    *tip++ = (xcol>>8) & 0xff;
		    *tip++ =  xcol & 0xff;
		}
		rp2 += pen_image->bytes_per_line;
            }

        } else {  /* LSBFirst */

            for (i=0, rp2=raster2; i<height; i++) {
		rp3 = raster3 + i*widthpad;
		for (j=0, tip=rp2; j<width; j++, rp3++) {
		    xcol = map[*rp3];
		    *tip++ =  xcol & 0xff;
		    *tip++ = (xcol>>8) & 0xff;
		    *tip++ = (xcol>>16) & 0xff;
		    if (do32) *tip++ = 0;
		}
		rp2 += pen_image->bytes_per_line;
            }
        }

    } else {
       
	/* this is a hard case, we will do it the slow way using XPutPixel */
    	depth = visual_depth;
   	format = ZPixmap;
   	xoffset = 0;
	xwidth = width;
	xheight = height;
    	/*bitmap_pad = 8;*/
    	bitmap_pad = 32;
    	bytes_per_line = 0;
	raster2 = (unsigned char *) calloc ((widthpad / 2) * 	
					    visual_depth * height / 8, 2); 
    	pen_image = XCreateImage(pen_display,pen_visual,depth,format,xoffset,
				 (char *)raster2,xwidth,xheight,bitmap_pad,bytes_per_line);
	for (i = 0; i < height; i++) {
	    rp = raster3 + widthpad*i;
	    for (j = 0; j < width; j++) {
		XPutPixel( pen_image , j, i, map[*rp++] );
	    }}
	
    }

    XPutImage(pen_display,pen_drawable,pen_gc,pen_image,
	      0,0,xloc,yloc,xwidth,xheight);
    XDestroyImage(pen_image);

    /* free up remaining work space */
    if (xmono || visual_depth != 8 ) {
	free(raster3 );
    }
    raster2 = raster3 = 0;
}
