#ifndef _vp_vplot_h
#define _vp_vplot_h

#include <stdbool.h>

enum {
    RPERIN=600,         /* vplot units per inch */
    HATCHPERIN=100,	/* Hatch units per inch */
    TXPERIN=33,	        /* Text units per inch */
    FATPERIN=200,	/* Fatness units per inch */
    MAX_GUN=255,	/* Maximum color gun strength */
    /*
     * This is the factor we scale our path and up vectors by before
     * running them through the local text coordinate transformation.
     * (The coordinate transformation, being in terms of device units,
     * gets done in integers. If we don't scale up we get severe roundoff
     * problems for small text sizes at odd angles. We can't make this
     * factor too big, though, or we risk very large text overflowing
     * the maximum possible integer.)
     */
    TEXTVECSCALE=10,
};

enum {
    VP_SETSTYLE         = 'S',
    VP_MOVE             = 'm',
    VP_DRAW	        = 'd',
    VP_PLINE	    	= 'L',
    VP_PMARK	   	= 'M',
    VP_TEXT		= 'T',
    VP_GTEXT		= 'G',
    VP_AREA		= 'A',
    VP_OLDAREA		= 'a',
    VP_BYTE_RASTER	= 'R',
    VP_BIT_RASTER	= 'r',
    
    VP_ERASE		= 'e',
    VP_BREAK		= 'b',
    VP_PURGE		= 'p',
    VP_NOOP		= 'n',
    
    VP_ORIGIN		= 'o',
    VP_WINDOW		= 'w',
    
    VP_FAT		= 'f',
    VP_SETDASH		= 's',
    VP_COLOR		= 'c',
    VP_SET_COLOR_TABLE	= 'C',
    VP_TXALIGN		= 'J',
    VP_TXFONTPREC	= 'F',
    VP_PATLOAD		= 'l',
    VP_OVERLAY		= 'v',
    
    VP_MESSAGE		= 'z',
    VP_BEGIN_GROUP	= '[',
    VP_END_GROUP	= ']',
    
    VP_OLDTEXT		= 't'
};

typedef enum {STANDARD=0,
	      ROTATED=1,
	      ABSOLUTE=3} vp_plotstyle;

/* for debuging purposes */
void vp_putint (int w);
int vp_getint (void);
void vp_putfloat (float w);
void vp_plot (float x, float y, bool down);
void vp_clip (float xmin, float ymin, float xmax, float ymax);
void vp_move (float x,float  y);
void vp_umove (float x,float  y);
void vp_draw (float x,float  y);
void vp_udraw (float x,float  y);
void vp_utext (float x, float y, int size, int orient, char *string);
void vp_text (float x, float y, int size, int orient, char *string);
void vp_orig (float x,float  y);
void vp_uorig (float x,float  y);
void vp_color (int col);
void vp_fat (int f);
void vp_scale (float xscale, float  yscale);
void vp_tjust (int xjust1, int yjust1);
void vp_arrow (float x1, float y1, float x, float y, float r);
void vp_uarrow (float x1, float y1, float x, float y, float r);

#endif
