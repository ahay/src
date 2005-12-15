#include "device.h"

#ifndef _vp_device_h

typedef struct Device *device;
/*^*/

typedef enum {BREAK_BREAK, BREAK_ERASE, BREAK_IGNORE} brake;
/*^*/

typedef enum {SET_COLOR=1,
	      SET_COLOR_TABLE,
	      SET_WINDOW,
	      NEW_DASH,		
	      NEW_PAT,		
	      NEW_FONT,		
	      NEW_OVERLAY,	
	      NEW_ALIGN,	
	      NEW_FAT,		
	      BEGIN_GROUP,	
	      END_GROUP} attr_command;
/*^*/


struct Device {
    float dashpos, dashsum, *dashes;
    float pixels_per_inch, aspect_ratio;
    float greyc, pixc;
    int invras;
    int afat;
    int xorigin, yorigin, xmin, ymin;
    float xscale, yscale, hscale, vscale;
    float mxx, mxy, myx, myy;
    float hshift, vshift;
    int no_stretch_text;
    int first_time;
    int mono;
    float redpow, greenpow, bluepow;
    float redmap[4], greenmap[4], bluemap[4];
    brake brk;
    /* control routines */
    int (*open)(void);
    int (*reset)(void);
    int (*message)(void);
    int (*erase)(void);
    int (*close)(void);    
    /* high level output */
    int (*vector)(int x1, int y1, int x2, int y2, int nfat, int dashon);
    int (*marker)(void);
    int (*text)(void);
    int (*area)(void);
    int (*raster)(void);
    int (*point)(void);
    int (*attributes)(device dev,attr_command,int,int,int,int);    
    /* input */
    int (*getpoint)(void);
    int (*interact)(void);  
    /* low level output */
    int (*plot)(void);
    int (*startpoly)(void);
    int (*midpoly)(void);
    int (*endpoly)(void);
};
/*^*/

#endif
