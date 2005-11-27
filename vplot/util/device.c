#include "device.h"

#ifndef _vp_device_h

typedef struct Device *device;
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
    int (*attributes)(void);    
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
