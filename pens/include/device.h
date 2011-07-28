#ifndef _pen_device_h
#define _pen_device_h

/*
 * device routine table structure
 */

struct device{
    int xmax, ymax, xmin, ymin;
    float pixels_per_inch;
    float aspect_ratio;
    int num_col;
    int lost;
    bool need_end_erase; /* Does device need erase at end? */
    bool smart_clip;     /* Can the device do its own clipping? (of vectors and polygons.) */
    bool smart_raster;   /* Can the device stretch AND clip its own raster? */
    bool smart_background; /* Can do smart background erase */
    bool cachepipe;
    /* Setting cachepipe = YES will copy any piped files to a temporary file,
     * this may get done in dev.open.
     * This is useful if the program may want to reverse seek. */
    int txfont,txprec,txovly;
    int xorigin, yorigin;	/* global "origin" */
    int brake;

    /* control routines */
    void (*open)(int argc, char* argv[]);
    void (*reset)(void);
    void (*message)(int command, const char *string);
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
    void (*reader)(int nn, FILE **inpltin, char** innames);
    int (*getpoint)(int *x, int *y);
    int (*interact)(int what, FILE *controltty, char *string);

    /* low level output
     * - these are called by the generic high level output routines */
    void (*plot)(int x, int y, int draw);
    void (*startpoly)(int npts);
    void (*midpoly)(int x, int y);
    void (*endpoly)(int last);
};

struct s_txalign {
    int hor;
    int ver;
};

#endif
