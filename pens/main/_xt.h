/* include files */
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

#include <stdio.h>
#include <string.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>


#include "_vp.h"
#include "../include/params.h"
#include "../include/enum.h"
#include "../include/err.h"
#include "../include/extern.h"

/* define a structure to reference the input files */
typedef struct
{
	FILE	*stream;
	char	name[MAXFLEN+1];
	int	fnew;
} VplFileInfo;
 
/* global definitions */
#define	TBORDER 22
#define	BORDER 5

/* global variables */
extern Display *pen_display;
extern int screen_depth;
extern int visual_depth;
extern long screen_black, screen_white;
extern Visual *pen_visual;
extern Window pen_window;
extern GC pen_gc;
extern XtAppContext pen_context;
extern Widget xtpen, vpane, buttonbox, pen_picture, text_region;
extern unsigned long *pen_colors;
extern int xt_endframe;

/* picture scaling info */
extern unsigned long color;
extern int x_num_col;
extern int xmono;
extern int xtruecol; 
extern int skipit;
extern int xt_paused;
extern int xt_after_erase;

extern VplFileInfo *inFiles;
extern int mouse_clicked,cursx,cursy;

extern Drawable pen_drawable; /* drawable, could be the window or a pixmap */

extern int in_repaint;        /* Are we in the repaint proc ? */
extern int num_files;

/* Button functions */
extern void create_buttons();
extern void activate_buttons();
extern void inactivate_buttons();
extern int epause;   /* user sepcified pause in seconds */
extern float fpause; /* pause between frames in seconds */

/* "button press" values returned by xt_pause  */
enum {NEXT,PREV,RESTART,QUIT,RUN,STOP,CHANGED,TIME};  

/* "run_mode" values set by the mode function in xtcommands*/
extern int xt_run_mode;
enum {XT_FORWARD,XT_BACKWARD,XT_BOTH_WAYS};  

/*
 *
 *  source file:   ./filters/xtlib/xtpixmap.h
 *
 * Steve Cole (SEP), August 4 1991
 *	Wrote xtpen.  
 *	Inserted this sample edit history entry.
 */

extern Pixmap pen_pixmap;
extern int pen_width,pen_height;

extern int have_pixmap;

#if defined(__STDC__) || defined (__stdc__)
extern Pixmap MyCreatePixmap( Display *, Drawable,unsigned int,
		unsigned int,unsigned int );
#else
extern Pixmap MyCreatePixmap();
#endif

/* inline coordinate transforms;  X origin- upper left, pen origin- lower left */
#define	XCORD(x)	x
#define	YCORD(y)	(dev.ymax - y)
#define	IXCORD(x)	x
#define	IYCORD(y)	(dev.ymax - y)
 
extern Colormap pen_colormap;

/* pixel map for color manipulation */
#define SEP_MAP_SIZE 65536
extern unsigned long map[SEP_MAP_SIZE];
/* do we own the whole colormap ? */
extern int own_colormap;

/* plotting_started is set on first entry to xt_dovplot */
extern int plotting_started;

/* size of the message buffer */
#define TEXT_BUFFER_SIZE 4096

/* functions for dealing with the list of frames  and getting old images */

struct cmap_ {
unsigned char   red[65536] ;
unsigned char   green[65536];
unsigned char   blue[65536];
unsigned long  map[65536];
};

typedef struct xtf_ {
    int file_num;		/* file to read (in list of files) */
    long  file_position;        /* starting poition in the file */
    int  end_file;              /* file in which the frame ends */
    long  end_pos;              /* position in that file */
    long total_len;             /* total length in bytes */
    int   frame_num;            /* number of this frame */
    int has_image;              /* A stored image of this frame exists */
    int break_end;              /* This frame ended with a break (not erase)*/
    char filename[256];         /* File name */
    XImage* image;		/* An XIMage of this frame */
    Pixmap  pixmap;		/* A pixmap of this frame */
    struct cmap_ cmap;		/* the colormap for this frame */
    struct xtf_ *next;		/* The next frame */
    struct xtf_ *prev;		/* The previous frame */
  
    
} xtFrame ;


/* 
 * A list of frames that can be plotted, the parent list of which this list
 * is an orderd subset is the parent member.
 */

typedef struct xtfl_ {
    xtFrame* start;
    xtFrame* end;
    int num_frame;
    struct xtfl_ *parent;
} xtFrameList;

extern xtFrameList *frame_list; /* current list of frames to draw */
extern xtFrame *redraw_frame;  /* the frame to redraw not necessarily the same 
	     			 as the current frame  because pauses occur
				 at the  start of the next frame */

/* appplication data to be read from the resource database */
typedef struct {
        Boolean stretchy;
        Boolean pixmaps;
        Boolean images;
        Boolean buttons;
        Boolean labels;
        Boolean textpane;
        int num_col;
        int vis_depth;
        float pause;
} AppData, *AppDataPtr;

extern AppData app_data;



