#ifndef _pen_device_h
#define _pen_device_h

/*
 * device routine table structure
 */

struct device{

    /* control routines */
    void (*open)(void);
    void (*reset)(void);
    void (*message)(int command, char *string);
    void (*erase)(int command);
    void (*close)(int status);
    
    /* high level output */
    void (*vector)(int x1, int y1, int x2, int y2, int nfat, int vpdashon);
    void (*marker)(int npts, int type, int size, int *pvec);
    void (*text)(char *string, float pathx, float pathy, float upx, float upy);
    void (*area)(int npts, struct vertex  *head);
    void (*raster)(int xpix, int ypix, int xmin, int ymin, int xmax, int ymax, 
		   unsigned char **raster_block, int orient, int dither_it);
    void (*point)(int x1, int y1);
    void (*attributes)(int command, int value, int v1, int v2, int v3);

    /* input */
    int (*getpoint)(int *x, int *y);
    int (*interact)(int what, FILE *controltty, char *string);

    /* low level output
     * - these are called by the generic high level output routines */
    void (*plot)(int x, int y, int draw);
    void (*startpoly)(int npts);
    void (*midpoly)(int x, int y);
    void (*endpoly)(int last);
};

#endif
