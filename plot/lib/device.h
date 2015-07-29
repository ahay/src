#ifndef _vp_device_h
#define _vp_device_h

#include <rsf.h>

enum {
    VP_DEFAULT_COLOR=7,
    VP_DEFAULT_LINESTYLE=0,
    VP_DEFAULT_FONT=0,
    VP_DEFAULT_HARDCOPY_FONT=3,
    VP_DEFAULT_PREC=0,
    VP_DEFAULT_HARDCOPY_PREC=2
};

typedef enum {SET_COLOR, SET_COLOR_TABLE, SET_WINDOW, NEW_DASH,
	      NEW_PAT, NEW_FONT, NEW_OVERLAY, NEW_ALIGN, NEW_FAT, 
	      BEGIN_GROUP, END_GROUP} vp_attribute;

typedef enum {CLOSE_DONE=-1,
	      CLOSE_NOTHING,
	      CLOSE_ERROR,
	      CLOSE_INTERRUPT,
	      CLOSE_NORMAL,
	      CLOSE_PAUSE,
	      CLOSE_FLUSH} vp_close;

typedef enum {INT_USER_PAUSE,
	      INT_GET_STRING,
	      INT_PAUSE,
	      INT_F_PAUSE} vp_interact;

typedef enum {ERASE_START,
	      ERASE_MIDDLE,
	      ERASE_END,
	      ERASE_BREAK} vp_eras;

struct vp_vertex {
    int x, y;
    struct vp_vertex *next, *last, *soft;
};

typedef struct vp_Device {
    int xwmax, xwmin, ywmax, ywmin; /* window */
} *vp_device;

/* device methods
   void (*open)(vp_device);
   void (*reset)(vp_device);
   void (*close)(vp_device,vp_close);
   void (*eras)(vp_device,vp_eras);
   void (*attributes)(vp_device,vp_attribute,int,int,int,int);
   void (*vector)(vp_device,int,int,int,int,int,bool);
   void (*marker)(vp_device,int,int,int,int*);
   void (*text) (vp_device,char*,float,float,float,float);
   void (*area) (vp_device,int,struct vp_vertex*);
   void (*raster) (vp_device,int,int,int,int,int,int,unsigned char*,int,int);
   void (*startpoly) (vp_device,int);
   void (*midpoly)(vp_device,int,int);
   void (*endpoly)(vp_device,bool);
*/

void vp_main(int argc, char* argv[], vp_device dev, 
	     void (*open)(vp_device),
	     void (*reset)(vp_device),
	     void (*close)(vp_device,vp_close),
	     void (*eras)(vp_device,vp_eras),
	     void (*attributes)(vp_device,
				vp_attribute,int,int,int,int),
	     void (*vector)(vp_device,int,int,int,int,int,bool),
	     void (*marker)(vp_device,int,int,int,int*),
	     void (*text) (vp_device,char*,float,float,float,float),
	     void (*area) (vp_device,int,struct vp_vertex*),
	     void (*raster) (vp_device,int,int,int,int,int,int,
			     unsigned char*,int,int));

#endif

