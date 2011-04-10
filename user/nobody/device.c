#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include <unistd.h>

#include <rsf.h>

#include "device.h"
#include "vplot.h"

#define NHATCH	20 	/* Maximum number of lines in a hatch pattern */
#define MAX_COL 511	/* Maximum color table number, maximum color number */
#define NPAT	512	/* Maximum internal pattern number, MAX_COL + 1 */
#define MAXDASH 10	/* Maximum number of dot-dashes */
#define NUM_PRIM  6     /* Number of Color table array elements */
#define MAX_GUN	255	/* Maximum color gun strength */

#define TXBUFLEN 250	/* Initial max length of a text string */
#define VXBUFLEN 250	/* Initial max number of vertices in a polygon */

/* BSD - MAXNAMELEN, Posix - NAME_MAX */
#ifndef NAME_MAX 
#define	NAME_MAX MAXNAMELEN
#endif

static bool need_devcolor = false;
static int xwmax_last, xwmin_last, ywmax_last, ywmin_last;
static int xmax = 0;
static int ymax = 0;
static int xmin = 0;
static int ymin = 0;
static float ppi = 0.;
static float aspect_ratio = 0.;
static int num_col = -1;

/*
static bool need_end_erase = false;
static bool buffer_output = true;
static bool buffer_input = true;
static bool allow_pipe = true;
static bool smart_clip = false;
*/

static bool smart_raster = false;

static float fatmult = 1.;
static float patternmult = 1.;

static char wstype[25];

/*
static char callname[25];
*/

static int xcenterflag, ycenterflag;

static bool ever_called = false;
static bool first_time = true;
static bool device_open = false;	/* used by ERR */

static int nplots = 0;	/* number of plots made */
static bool no_stretch_text = true;	/* Don't allow stretched text? */

static int xnew, ynew;	/* new pen location */
static int xold, yold;	/* old pen location */
static int xorigin = 0, yorigin = 0;	/* global "origin" */
static char *txbuffer;
static int txbuflen, vxbuflen;

struct vp_vertex *vxbuffer;

/*
static int xret, yret;
*/

static float mxx, mxy, myx, myy;

static int default_rotate = 0, default_hshift = 0, default_vshift = 0;
static float default_scale = 1., default_xscale = 1., default_yscale = 1.;

/*
static int linestyle;
*/

static int cur_color = VP_DEFAULT_COLOR;
static int pat_color = VP_DEFAULT_COLOR + 1;
static int next_color;

struct txalign {
    int hor, ver;
} txalign;

static int txfont = VP_DEFAULT_FONT;
static int txprec = VP_DEFAULT_PREC;
static int txovly = OVLY_NORMAL;
static int default_txfont, default_txprec, default_txovly;
static int fat = 0;
static int afat = 0;
static int ipat = 0;

struct pat /* area fill pattern */
{
    int xdim, ydim, xdim_orig, ydim_orig, *patbits;
} pat[NPAT + 1];

static float dashsum = 0.;
static bool dashon = false;
static float dashes[MAXDASH * 2];
static float dashpos = 0.;	/* Position in current dashing pattern */
static int color_set[MAX_COL + 1][NUM_PRIM];
static int num_col_8;

static bool mono;
/*
 * Invras determines the polarity of dithered raster.
 * invras=n means raster works the same way as vectors; what is normally
 * WHITE on most devices is BLACK on a hardcopy black and white device.
 * invras=y is the proper default to make dithered images not come out 
 * as negatives
 */
static bool invras = true;
static bool window = true;
static bool shade = true;

typedef enum {BREAK_BREAK, BREAK_ERASE, BREAK_IGNORE} breaks;

static breaks brake = BREAK_BREAK;
static bool framewindows = false;
static bool endpause = false;
static bool cachepipe = false;
/* Setting cachepipe = true will copy any piped files to a temporary file */

static bool wantras = true;
static bool colormask[5] = {true, true, true, true, true};

static vp_plotstyle style = VP_NO_STYLE_YET;
static vp_plotstyle default_style = VP_STANDARD;
static int rotate;

typedef enum {RELATIVE, ABSOLUTE} sizes;

static sizes size = RELATIVE;

enum {DO_NOTHING, DO_LITERALS=1, FORCE_INITIAL};

static int erase = FORCE_INITIAL | DO_LITERALS;

static int xcenter, ycenter;	/* Vplot of point to force as center */
static int fatbase = 0;
static int epause = 0;	/* time to pause before erasing screen */
static bool overlay = false;	/* 1=overlay 0=replace */
static bool default_overlay;

/*
 * 0 = none
 * 1 = random
 * 2 = Ordered
 * 3 = Floyd-Steinberg
 * 4 = Halftone
 */
static int dither = 1;	/* Dithering type */

static int xWmax, xWmin, yWmax, yWmin;
static float hshift, vshift;	/* Allow global translation of plot */
static float scale;	/* global scale */
static float xscale;	/* global x-scale */
static float yscale;	/* global y-scale */
static float hdevscale;	/* Vplot units to device units for x */
static float vdevscale;	/* Vplot units to device units for y */
static float txscale;/* global text scale */
static float mkscale;/* global marker scale */
static float dashscale;	/* global dashed line scale */
static char *interact; /* Coordinate file */
static float greyc = 1.;	/* Nonlinear correction */
static float pixc = 1.;	/* Pixel overlap correction */
static float redpow = 1., redmap[4] = {1., 0., 0., 0.};
static float greenpow = 1., greenmap[4] = {0., 1., 0., 0.};
static float bluepow = 1., bluemap[4] = {0., 0., 1., 0.};

/* filter options - resettable between plots */
static int user_rotate;
static float user_txscale;
static float user_mkscale;
static float user_dashscale;
static float user_scale;
static float user_xscale;
static float user_yscale;
static int user_size;
static float user_hshift;
static float user_vshift;
static int user_xwmax_flag;
static float user_xwmax;
static int user_ywmax_flag;
static float user_ywmax;
static int user_xwmin_flag;
static float user_xwmin;
static int user_ywmin_flag;
static float user_ywmin;

static int ifat = 0;
static float fatmult_orig;

/**************
static int     (*message) () = genmessage;
***************/

/*
static char group_name[NAME_MAX + 1];
static int group_number = 0;
*/

/*************
static int     (*genreader) () = gen_do_dovplot;
**************/

enum {STATUS, MAP, RED, GREEN, BLUE, GREY};

enum {UNSET, SET, MAPPED, MAP_SET};

static int greycorr (int colornum);
static void reset_parameters (vp_device dev, 
			      void (*attributes)(vp_device,vp_attribute,
						 int,int,int,int));
static void vptodevxy (int x, int y, int *outx, int *outy);
static void vptodevw (int x1, int y1, int x2, int y2, 
		      int *x1out, int *y1out, int *x2out, int *y2out);
static void devtovpxy (int x, int y, int *outx, int *outy);
static void devtovpw (int x1, int y1, int x2, int y2, 
		      int *x1out, int *y1out, int *x2out, int *y2out);
static void vptodevxy_text (int x, int y, int *outx, int *outy);
static void reset_windows (vp_device dev, 
			   void (*attributes)(vp_device,
					      vp_attribute,int,int,int,int));
static void init_colors (void);
static void outline_window (vp_device dev,
			    void (*attributes)(vp_device,
					       vp_attribute,int,int,int,int),
			    void (*vector)(vp_device,
					   int,int,int,int,int,bool));
static void setstyle (vp_plotstyle new_style, vp_device dev,
		      void (*attributes)(vp_device,
					 vp_attribute,int,int,int,int));
static void update_color (vp_device dev, 
			  void (*attributes)(vp_device,
					     vp_attribute,int,int,int,int));
static void getxy (int *x, int *y);
static void getxy_text (int *x, int *y);
static void getvpstring (void);
static int getpolygon (int npts);
static bool dupside (struct vp_vertex *base);
static void vecoutline (vp_device dev, struct vp_vertex *head, 
			void (*vector)(vp_device,int,int,int,int,int,bool));
static void reset_all (vp_device dev,
		       void (*attributes)(vp_device,
					  vp_attribute,int,int,int,int));
static void dithline (unsigned char *inpline, 
		      unsigned char *outline, 
		      int npixels, int linenum, int imethod);
static void drawpolygon (vp_device dev, int npts, int *x, int *y,
			 void (*attributes)(vp_device,vp_attribute,
					    int,int,int,int),
			 void (*area) (vp_device,int,struct vp_vertex*),
			 void (*vector)(vp_device,int,int,int,int,int,bool));
static void vplot_init (vp_device dev, void (*open)(vp_device));
static void vplot_do (vp_device dev, 
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

static void vplot_init (vp_device dev, void (*open)(vp_device))
{
    char *string;
    int i;
    float ftemp;

    txbuffer = sf_charalloc (TXBUFLEN);
    txbuflen = TXBUFLEN;

    vxbuflen = VXBUFLEN;
    vxbuffer = (struct vp_vertex *) sf_alloc (vxbuflen, sizeof (*vxbuffer));

    if (!sf_getbool("mono",&mono)) mono=false;

    if (NULL != (string = sf_getstring ("wstype")) ||
	NULL != (string = getenv ("WSTYPE"))) {
	strncpy (wstype, string, 24);
    } else {
	strncpy (wstype, "default",24);
    }

    for (i = 0; i <= NPAT; i++) {
	pat[i].patbits = NULL;
    }

    open (dev);
    device_open = true;

    if (!sf_getbool ("endpause",&endpause)) endpause=false;
    if (!sf_getbool ("cachepipe", &cachepipe)) cachepipe=false;
    if (!sf_getbool ("shade", &shade)) shade=true;
    if (!sf_getbool ("wantras", &wantras)) wantras=true;
    if (!sf_getbool ("window", &window)) window=true;
    if (!sf_getbool ("frame", &framewindows)) framewindows=false;
    if (!sf_getbool ("overlay", &overlay)) overlay=false;

    default_overlay = overlay;
    if (!sf_getbool ("invras", &invras)) invras=true;
    if (!sf_getbool ("txsquare", &no_stretch_text)) no_stretch_text=true;

    if (!sf_getbools ("colormask", colormask, 5)) {
	colormask[0] = 
	    colormask[1] = colormask[2] = colormask[3] = colormask[4] = true;
    } else if (!colormask[4]) {
	if (!colormask[0] && !colormask[1] && !colormask[2] && !colormask[3])
	    sf_error ("%s: Possibly invalid colormask=",__FILE__);
	if (colormask[3]) colormask[4] = true;
    }
    
    if (!sf_getfloat ("redpow",&redpow)) redpow=1.;
    if (!sf_getfloat ("greenpow",&greenpow)) greenpow=1.;
    if (!sf_getfloat ("bluepow",&bluepow)) bluepow=1.;

    if (redpow <= 0. || greenpow <= 0. || bluepow <= 0.)
	sf_error ("%s: Invalid (color)pow=",__FILE__);

    if (!sf_getfloats ("red", redmap,4)) {
	redmap[0]=1.; redmap[1]=redmap[2]=redmap[3]=0.;
    }
    if (!sf_getfloats ("green", greenmap,4)) {
	greenmap[1]=1.; greenmap[0]=greenmap[2]=greenmap[3]=0.;
    }
    if (!sf_getfloats ("blue", bluemap, 4)) {
	bluemap[2]=1.; bluemap[0]=bluemap[1]=bluemap[3]=0.;
    }

    if (!sf_getint ("dither",&dither)) dither=1;
    if (!sf_getfloat ("greyc",&greyc)) greyc=1.;
    if (!sf_getfloat ("pixc",&pixc)) pixc=1.;

    if (!sf_getint("txfont",&txfont)) txfont= VP_DEFAULT_FONT;
    if (!sf_getint ("txprec",&txprec)) txprec= VP_DEFAULT_PREC;
    if (!sf_getint ("txovly",&txovly)) txovly= OVLY_NORMAL;
    default_txfont = txfont;
    default_txprec = txprec;
    default_txovly = txovly;

    if (NULL != (string = sf_getstring("erase"))) {
	switch (string[0]) {
	    case 'n':
	    case 'N':
		erase = DO_NOTHING;
		break;
	    case 'o':
	    case 'O':
		erase = FORCE_INITIAL;
		break;
	    case 'l':
	    case 'L':
		erase = DO_LITERALS;
		break;
	    default:
		erase = FORCE_INITIAL | DO_LITERALS;
		break;
	}
    }

    if (NULL != sf_getstring ("break")) {
	switch (string[0]) {
	    case 'b':
		brake = BREAK_BREAK;
		break;
	    case 'e':
		brake = BREAK_ERASE;
		break;
	    default:
		brake = BREAK_IGNORE;
		break;
	}
    }

    xcenter = 0;
    ycenter = 0;
    xcenterflag = false;
    ycenterflag = false;
    if (sf_getfloat ("xcenter",&ftemp)) {
	xcenterflag = true;
	ftemp *= RPERIN;
	xcenter = (ftemp<0)? ftemp-0.5: ftemp+0.5;
    }
    if (sf_getfloat ("ycenter",&ftemp)) {
	ycenterflag = true;
	ftemp *= RPERIN;
	ycenter = (ftemp<0)? ftemp-0.5: ftemp+0.5;
    }

    if (!sf_getfloat ("patternmult",&patternmult)) {
	if (NULL != (string = getenv ("PATTERNMULT"))) {
	    sscanf (string, "%f", &patternmult);
	} else {
	    patternmult = 1.;
	}
    }


    if (!sf_getint ("pause",&epause)) epause = 0;

    if (NULL == (interact = sf_getstring("interact"))) {
	interact = "";
    } else if ('\0' != interact[0]) {
	epause = 0;		/* interact makes it own sort of pausing */
	endpause = false;
    }

    if (NULL != (string = sf_getstring("style")) ||
	NULL != (string = getenv ("PLOTSTYLE"))) {
	switch (string[0]) {
	    case 'r':
	    case 'R':
	    case 'm':
	    case 'M':
		default_style = VP_ROTATED;
		break;
	    case 'a':
	    case 'A':
		default_style = VP_ABSOLUTE;
		break;
	    default:
		default_style = VP_STANDARD;
		break;
	}
    }

    if (NULL == (string = sf_getstring("fatmult")) &&
	NULL == (string = getenv ("FATMULT"))) 
	fatmult=1.;
    fatmult_orig = fatmult;

    if (!sf_getint ("rotate",&user_rotate)) user_rotate=0;

    if (!sf_getfloat("txscale", &user_txscale)) user_txscale = 1.0;
    if (!sf_getfloat("mkscale", &user_mkscale)) user_mkscale = 1.0;
    if (!sf_getfloat("dashscale", &user_dashscale)) user_dashscale = 1.0;

    if (!sf_getfloat ("scale", &user_scale)) user_scale = 1.0;
    if (!sf_getfloat ("xscale", &user_xscale)) user_xscale = 1.0;
    if (!sf_getfloat ("yscale", &user_yscale)) user_yscale = 1.0;

    if (NULL != sf_getstring ("size")) {
	switch (string[0]) {
	    case 'a':
	    case 'A':
		user_size = ABSOLUTE;
		break;
	    default:
		user_size = RELATIVE;
		break;
	}
    } else {
	user_size = size;
    }

    if (!sf_getfloat ("hshift", &user_hshift)) user_hshift = 0.; 
    if (!sf_getfloat ("vshift", &user_vshift)) user_vshift = 0.;

    user_xwmax_flag = sf_getfloat ("xwmax", &user_xwmax);
    user_ywmax_flag = sf_getfloat ("ywmax", &user_ywmax);
    user_xwmin_flag = sf_getfloat ("xwmin", &user_xwmin);
    user_ywmin_flag = sf_getfloat ("ywmin", &user_ywmin);

    if (!sf_getint ("fat", &fatbase)) fatbase=0;
}

static void init_colors (void)
{
    int i;

    for (i = 0; i < 8; i++) {
	color_set[i][STATUS] = SET;
	color_set[i][RED] = MAX_GUN * ((i & 2) / 2);
	color_set[i][GREEN] = MAX_GUN * ((i & 4) / 4);
	color_set[i][BLUE] = MAX_GUN * ((i & 1) / 1);
	color_set[i][GREY] = greycorr ((MAX_GUN * i) / 7);
    }

    num_col_8 = (num_col > 8) ? num_col : 8;

    if (mono) {

	/* Monochrome devices are assumed to have no settable colors */
	num_col = 0;

	color_set[0][MAP] = 0;
	for (i = 1; i <= MAX_COL; i++) {
	    color_set[i][MAP] = 7;
	}
	for (i = num_col_8; i < 256; i++) {
	    color_set[i][GREY] = color_set[((i - 8) % 7) + 1][GREY];
	}
	/* Put a grey scale in the upper half of the color table */
	for (i = 256; i <= MAX_COL; i++) {
	    color_set[i][GREY] = greycorr (i - 256);
	}
    } else {
	/* Unmapped colors shouldn't be mapped; ie, they map to themselves */
	for (i = 0; i < num_col_8; i++) {
	    color_set[i][MAP] = i;
	}
	
	/* Colors outside the range of map cyclically back into 1 through 7 */ 
	for (i = num_col_8; i <= MAX_COL; i++) {
	    color_set[i][MAP] = ((i - 8) % 7) + 1;
	}
    }
}

static void setstyle (vp_plotstyle new_style, vp_device dev, 
		      void (*attributes)(vp_device,
					 vp_attribute,int,int,int,int))
{
    if (new_style != style) {
	style = new_style;
	reset_parameters (dev, attributes);
    }
}
    
static void reset_parameters (vp_device dev, 
			      void (*attributes)(vp_device,
						 vp_attribute,int,int,int,int))
{
    float inches, screenheight, screenwidth, f;
    int ix, iy;

    xorigin = 0;
    yorigin = 0;
    
    rotate = default_rotate;
    rotate += user_rotate;

    txscale = user_txscale;
    mkscale = user_mkscale;
    dashscale = user_dashscale;

    scale = default_scale * user_scale;
    xscale = default_xscale * user_xscale;
    yscale = default_yscale * user_yscale;

    fatmult = fatmult_orig;

    size = user_size;

    switch (style) {
	case VP_ROTATED:
	    rotate += 90;
	    inches = VP_ROTATED_HEIGHT;
	    break;
	case VP_ABSOLUTE:
	    size = ABSOLUTE;
	case VP_STANDARD:
	default:
	    inches = VP_STANDARD_HEIGHT;
	    break;
    }

    if (rotate >= 0) {
	rotate %= 360;
    } else {
	rotate = ((rotate % 360) + 360) % 360;
    }

    mxx = cosf (SF_PI * rotate / 180.);
    myy = cosf (SF_PI * rotate / 180.);
    mxy = sinf (SF_PI * rotate / 180.);
    myx = -sinf (SF_PI * rotate / 180.);

    if (ABSOLUTE == size) {
	vdevscale = ppi / (RPERIN * aspect_ratio);
	hdevscale = ppi / RPERIN;
    } else {
	/* Fit into a displayable box with aspect ratio SCREEN_RATIO */

	screenwidth  = (xmax - xmin) * ppi;
	screenheight = (ymax - ymin) * ppi * aspect_ratio;
	if ((screenheight / screenwidth) > VP_SCREEN_RATIO) {
	    vdevscale = VP_SCREEN_RATIO * (xmax - xmin)/aspect_ratio /
		(inches * RPERIN);
	    hdevscale = vdevscale * aspect_ratio;
	} else {
	    vdevscale = (ymax - ymin) / (inches * RPERIN);
	    hdevscale = vdevscale * aspect_ratio;
	}
    }

    hshift = default_hshift;
    vshift = default_vshift;

    if (VP_ROTATED == style) vshift += ymax - ymin;

    yscale *= scale;
    xscale *= scale;
    mkscale *= scale;

/* fatness multiplication factor */
    fatmult *= scale * hdevscale * RPERIN / FATPERIN;

    /* (xcenter,ycenter) in vplot coordinates to be centered in the screen. */
    if (xcenterflag || ycenterflag) {
	vptodevxy (xcenter, ycenter, &ix, &iy);
	if (xcenterflag) hshift += (xmax + xmin) / 2 - ix;
	if (ycenterflag) vshift += (ymax + ymin) / 2 - iy;
    }

    hshift += user_hshift * ppi;
    vshift += user_vshift * ppi / aspect_ratio;

    devtovpw (xmin, ymin, xmax, ymax, &xWmin, &yWmin, &xWmax, &yWmax);
    
    if (user_xwmax_flag) {
	f = user_xwmax * RPERIN; xWmax = (f < 0.)? f-0.5: f+0.5;
    }
    if (user_ywmax_flag) {
	f = user_ywmax * RPERIN; yWmax = (f < 0.)? f-0.5: f+0.5;
    }
    if (user_xwmin_flag) {
	f = user_xwmin * RPERIN; xWmin = (f < 0.)? f-0.5: f+0.5;
    }
    if (user_ywmin_flag) {
	f = user_ywmin * RPERIN; yWmin = (f < 0.)? f-0.5: f+0.5;
    }

    vptodevw (xWmin, yWmin, xWmax, yWmax, &xWmin, &yWmin, &xWmax, &yWmax);

    if (xWmin < xmin) xWmin = xmin;
    if (yWmin < ymin) yWmin = ymin;
    if (xWmax > xmax) xWmax = xmax;
    if (yWmax > ymax) yWmax = ymax;

    /* plot window parameters defaulted to maximum size */
    dev->xwmax = xWmax;		
    dev->xwmin = xWmin;	
    dev->ywmax = yWmax;
    dev->ywmin = yWmin;
    
    reset_windows (dev, attributes);
}

/* Utility to modify color tables for plotting grey rasters. */
static int greycorr (int colornum)
{
    float newval;

    newval = invras? 255-colornum: colornum;

/* correction to simulate nonlinearity of graphics displays */
    if (greyc != 1.) {
	newval /= 255.;
	newval *= ((-2. + 2. * greyc) * newval + 3. * (1. - greyc)) * newval 
	    + greyc;
	newval *= 255.;

	if      (newval < 0)    newval = 0.;
	else if (newval > 255.) newval = 255.;
    }

/* correction for pixel overlap on hardcopy devices */
    if (pixc != 1.) {
	if (newval < pixc * 128.) {
	    newval /= pixc;
	} else {
	    newval = 128. + (newval - pixc * 128.) / (2. - pixc);
	}

	if (newval < 0)         newval = 0.;
	else if (newval > 255.) newval = 255.;
    }

    return ((int) (newval+0.5));
}

/* convert vplot coordinates to device coordinates, or vice versa */

static void vptodevxy (int x, int y, int *outx, int *outy)
{
    float fx, fy, f;
    
    fx = (x - xorigin) * xscale;
    fy = (y - yorigin) * yscale;

    f =  mxx * fx + mxy * fy;
    fy = myx * fx + myy * fy;
    fx = f;

    fx = fx * hdevscale + xmin + hshift;
    fy = fy * vdevscale + ymin + vshift;

    *outx = (fx < 0.)? fx-0.5: fx+0.5;
    *outy = (fy < 0.)? fy-0.5: fy+0.5;
}

static void vptodevw (int x1, int y1, int x2, int y2, 
		      int *x1out, int *y1out, int *x2out, int *y2out)
{
    int x11, y11, x12, y12, x21, y21, x22, y22, a, b;

    vptodevxy (x1, y1, &x11, &y11);
    vptodevxy (x1, y2, &x12, &y12);
    vptodevxy (x2, y1, &x21, &y21);
    vptodevxy (x2, y2, &x22, &y22);

    a = (x11 > x12 ? x11 : x12);
    b = (x22 > x21 ? x22 : x21);
    *x2out = (a > b ? a : b);

    a = (y11 > y12 ? y11 : y12);
    b = (y22 > y21 ? y22 : y21);
    *y2out = (a > b ? a : b);

    a = (x11 < x12 ? x11 : x12);
    b = (x22 < x21 ? x22 : x21);
    *x1out = (a < b ? a : b);

    a = (y11 < y12 ? y11 : y12);
    b = (y22 < y21 ? y22 : y21);
    *y1out = (a < b ? a : b);
}

static void devtovpxy (int x, int y, int *outx, int *outy)
{
    float fx, fy, f;

    fx = (x - xmin - hshift) / hdevscale;
    fy = (y - ymin - vshift) / vdevscale;

    f = mxx * fx - mxy * fy;
    fy = -myx * fx + myy * fy;
    fx = f;

    fx = fx / xscale + xorigin;
    fy = fy / yscale + yorigin;

    *outx = (fx < 0.)? fx-0.5:fx+0.5;
    *outy = (fy < 0.)? fy-0.5:fy+0.5;
}

static void devtovpw (int x1, int y1, int x2, int y2, 
		      int *x1out, int *y1out, int *x2out, int *y2out)
{
    int x11, y11, x12, y12, x21, y21, x22, y22, a, b;

    devtovpxy (x1, y1, &x11, &y11);
    devtovpxy (x1, y2, &x12, &y12);
    devtovpxy (x2, y1, &x21, &y21);
    devtovpxy (x2, y2, &x22, &y22);

    a = (x11 > x12 ? x11 : x12);
    b = (x22 > x21 ? x22 : x21);
    *x2out = (a > b ? a : b);

    a = (y11 > y12 ? y11 : y12);
    b = (y22 > y21 ? y22 : y21);
    *y2out = (a > b ? a : b);

    a = (x11 < x12 ? x11 : x12);
    b = (x22 < x21 ? x22 : x21);
    *x1out = (a < b ? a : b);

    a = (y11 < y12 ? y11 : y12);
    b = (y22 < y21 ? y22 : y21);
    *y1out = (a < b ? a : b);
}

static void vptodevxy_text (int x, int y, int *outx, int *outy)
{
    float xscale_save=0., yscale_save=0.;

    if (no_stretch_text) {
	xscale_save = xscale;
	yscale_save = yscale;

	if (fabsf (xscale) < fabsf (yscale)) {
	    yscale = xscale;
	} else {
	    xscale = yscale;
	}
    }

    vptodevxy (x, y, outx, outy);

    if (no_stretch_text) {
	xscale = xscale_save;
	yscale = yscale_save;
    }
}

static void vplot_do (vp_device dev, 
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
				      unsigned char*,int,int))
{
    int i, j=0, k, c, key, size, npts, nmul, nx, ny, nx_orig, ny_orig;
    int orient, ras_orient, *ptr=NULL, nx_mult, ny_mult, nx_temp, ny_temp;
    int *tempbuf, new_style;
    int col_tab_no, red, green, blue, grey, dist, min_dist, best_col=0;
    int red_mask, green_mask, blue_mask;
    int red_0, green_0, blue_0;
    int red_7, green_7, blue_7;
    int hacol[NHATCH * 2], hafat[NHATCH * 2], haoff[NHATCH * 2];
    int	hasiz[NHATCH * 2], numhatch;
    int xpix, ypix, num_pat, num_byte, yrast, lastrast;
    int xvr_min, xvr_max, yvr_min, yvr_max;
    int xvru_min, xvru_max, yvru_min, yvru_max;
    int xr_min=0, xr_max=0, yr_min=0, yr_max=0;
    int xru_min=0, yru_max=0, xxx[4], yyy[4];
    int pos, ii, jj, kk, num_rep, ras_offset, dither_it;
    int xnewer, ynewer;
    int xtext0, xtext1, xtext2, ytext0, ytext1, ytext2, type;
    int *marker_vec, *mvec;
    int savefat;
    unsigned char *rasterline, *rasterline2, *outraster, *outraster2, ibyte;
    float red_f, green_f, blue_f, red_fm, green_fm, blue_fm;
    float angle, xrasmult=0., yrasmult=0., savefatmult, f;
    bool starterase, skip=false;
       
    starterase = false;
    while ((VP_ERASE == (c = getchar ())) || 
	   ((VP_BREAK == c) && (brake == BREAK_ERASE))) {
	starterase = true;
    }

    if (c == EOF) return;

    if (first_time) {
	reset (dev);

	if (xmax <= xmin ||
	    ymax <= ymin ||
	    ppi == 0. ||
	    aspect_ratio == 0. ||
	    num_col == -1)
	    sf_error ("%s: wrong values for critical device variables",
		      __FILE__);

	/* Set maximum clipping window */
	xwmax_last = xmax;
	xwmin_last = xmin;
	ywmax_last = ymax;
	ywmin_last = ymin;

        /* Set up color maps */
	init_colors ();

	ever_called = true;
    } else {
	close (dev,CLOSE_FLUSH);
	if (epause > 0) {
	    sleep ((unsigned int) epause);
	} else if (epause < 0) {
	    close (dev,CLOSE_FLUSH);
	    
/*
	    message (MESG_ERASE);
	    message (MESG_ON);
	    message (MESG_HOME);
	    message (MESG_READY);
	    message (MESG_HIGHLIGHT_ON);
	    message (MESG_TEXT, "Type Return to Continue...  ");
	    message (MESG_DONE);

	    ii = dev.interact (INT_PAUSE, controltty, string);
	    if (ii == DOVPLOT_EXIT)
	    {
		dev.close (CLOSE_FLUSH);
		return;
	    }
	    message (MESG_HIGHLIGHT_OFF);
	    if (!allowecho)
	    {
		message (MESG_READY);
		message (MESG_TEXT, CRLF);
	    }
	    message (MESG_DONE);
	    message (MESG_OFF);
	    message (MESG_ERASE);
	}

	if (interact[0] != '\0')
	{
	    getapoint ();
	}
*/
	}
    }

    /* Erase the screen to background color */
    if (((erase & FORCE_INITIAL) && (first_time || (erase & DO_LITERALS)))
	|| (starterase && (erase & DO_LITERALS))) {
	if (first_time) {
	    eras (dev,ERASE_START);
	} else {
	    eras (dev,ERASE_MIDDLE);
	    nplots++;
	}
	close (dev,CLOSE_FLUSH);
    }

    first_time = false;


/* Reset fatness, cur_color, etc. */
    new_style = default_style;
    setstyle (new_style,dev,attributes);
    reset_all (dev, attributes);

/* Make SURE the color is what it's supposed to be, just to be safe. */
    attributes (dev,SET_COLOR, cur_color, 0, 0, 0);
    need_devcolor = false;

/*    message (MESG_OFF); */
    close (dev,CLOSE_FLUSH);

/* Start a group that will contain this frame (always group 0). */

    /* ?????????????? 
    ii = ftell (pltin);
    sprintf (group_name, "%s.%d", pltname, nplots);
    dev.attributes (BEGIN_GROUP, group_number, ii, 0, 0);
    group_number++;
    ?????????????????? */

    if (framewindows) outline_window (dev, attributes, vector);

    while (EOF != c) {
	switch (c)  {
	    case VP_SETSTYLE: 
		c = getchar ();
		switch (c) {
		    case 'r': case 'R': case 'm': case 'M':
			new_style = VP_ROTATED;
			break;
		    case 'a': case 'A':
			new_style = VP_ABSOLUTE;
			break;
		    default:
			new_style = VP_STANDARD;
			break;
		}
		setstyle (new_style,dev,attributes);

		if (framewindows) outline_window (dev, attributes, vector);
		break;
	    case VP_MOVE:  
		/* Reset position in dash pattern. */
		dashpos = 0.;

		getxy(&xold,&yold);
		break;
	    case VP_DRAW:
		getxy (&xnew,&ynew);
		update_color (dev, attributes);

		while (VP_DRAW == (c = getchar ())) {
		    getxy (&xnewer,&ynewer);

		    /* Is it the same point? */
		    if (xnewer == xnew && ynewer == ynew) continue;
		
                    /* Is it colinear and horizontal or vertical? */

		    if ((ynewer == ynew && ynew == yold &&
			 ((xnewer > xnew) == (xnew > xold))) ||
			(xnewer == xnew && xnew == xold &&
			 ((ynewer > ynew) == (ynew > yold)))) {
			ynew = ynewer;
			xnew = xnewer;
			continue;
		    }

		    vector (dev, xold, yold, xnew, ynew, fat, dashon);
		    xold = xnew;
		    yold = ynew;
		    xnew = xnewer;
		    ynew = ynewer;
		}

		vector (dev, xold, yold, xnew, ynew, fat, dashon);
		xold = xnew;
		yold = ynew;

		skip = true;
		break;
	    case VP_PLINE:

		/* Reset position in dash pattern. */
		dashpos = 0.;

		npts = vp_getint ();
		if (0 == npts) break;

		getxy (&xold, &yold);
		npts--;
		if (0 == npts) break;

		getxy (&xnew, &ynew);
		npts--;

		update_color (dev,attributes);
		while (npts > 0) {
		    getxy (&xnewer, &ynewer);
		    npts--;

		    /* Is it the same point? */
		    if (xnewer == xnew && ynewer == ynew) continue;

		    /* Is it colinear and horizontal or vertical? */
		    if ((ynewer == ynew && ynew == yold &&
			 ((xnewer > xnew) == (xnew > xold))) ||
			(xnewer == xnew && xnew == xold &&
			 ((ynewer > ynew) == (ynew > yold)))) {
			ynew = ynewer;
			xnew = xnewer;
			continue;
		    }

		    vector (dev, xold, yold, xnew, ynew, fat, dashon);
		    xold = xnew;
		    yold = ynew;
		    xnew = xnewer;
		    ynew = ynewer;
		}

		vector (dev, xold, yold, xnew, ynew, fat, dashon);
		xold = xnew;
		yold = ynew;
		break;
	    case VP_PMARK:		/* polymarker */
		npts = vp_getint (); /* how many markers ? */
		type = vp_getint (); /* what symbol? (any positive integer) */
		size = vp_getint (); /* How big? */

		size *= mkscale * vdevscale * RPERIN / TXPERIN;

		if (0 == npts) break;

		/* allocate space for the points */
		marker_vec = sf_intalloc (npts * 2);
		mvec = marker_vec;

		/* read locations, and transform them to device coordinates */
		for (ii = 0; ii < 2*npts; ii += 2) {
		    getxy (mvec+ii,mvec+ii+1);
		}

		update_color (dev, attributes);

		marker (dev, npts, type, size, marker_vec);

		free (marker_vec);
		break;
	    case VP_ORIGIN:	
		xorigin = vp_getint();
		yorigin = vp_getint();
		break;
	    case VP_BEGIN_GROUP:
		sf_error("%s: VP_BEGIN_GROUP not implemented",__FILE__);
		/*
		  ii = ftell (pltin) - 1;
		  getvpstring ();
		  strncpy (group_name, txbuffer, MAXFLEN);
		  dev.attributes (BEGIN_GROUP, group_number, ii, 0, 0);
		  group_number++;
		*/
		break;
	    case VP_END_GROUP:
		sf_error("%s: VP_END_GROUP not implemented",__FILE__);
		/*
		  group_number--;
		  if (group_number < 1)
		  {
		  ERR (WARN, name,
		  "group invalidly nested");
		  group_number = 1;
		  }
		  else
		  {
		  ii = dev.attributes (END_GROUP, group_number, USER_EGROUP, 0, 0);
		  if (ii == DOVPLOT_EXIT)
		  {
		  dev.close (CLOSE_FLUSH);
		  return;
		  }
		  if (ii != DOVPLOT_CONT)
		  ERR (WARN, name, "dev.attributes(END_GROUP,...) returned junk.");
		  } 
		*/
		break;		
	    case VP_GTEXT:
	    case VP_TEXT:
	    case VP_OLDTEXT:
		xtext0 = 0;
		ytext0 = 0;

		if (VP_GTEXT == c) {
		    vptodevxy_text (xtext0, ytext0, &xtext0, &ytext0);
		    getxy_text (&xtext1, &ytext1);
		    getxy_text (&xtext2, &ytext2);
		} else { 
		    if (VP_TEXT == c) {
			size = vp_getint ();
			size *= RPERIN / (float) TXPERIN;
			
			orient = vp_getint();
		    } else {
			if ((key = vp_getint ()) < 0)
			    sf_error ("%s: invalid text key",__FILE__);
			size = (key & 037);
			size *= RPERIN / (float) TXPERIN;
			orient = (((key & 0140) >> 5) * 90);
		    }

		    /* Character path direction */
		    f = TEXTVECSCALE * size * cosf (orient * SF_PI / 180.);
		    xtext1 = (f < 0.)? f-0.5:f+0.5;
		    f = TEXTVECSCALE * size * sinf (orient * SF_PI / 180.);
		    ytext1 = (f < 0.)? f-0.5:f+0.5;

		    /* Character up vector direction */
		    orient += 90;
		    f = TEXTVECSCALE * size * cosf (orient * SF_PI / 180.);
		    xtext2 = (f < 0.)? f-0.5:f+0.5;
		    f = TEXTVECSCALE * size * sinf (orient * SF_PI / 180.);
		    ytext2 = (f < 0.)? f-0.5:f+0.5;

		    vptodevxy_text (xtext0, ytext0, &xtext0, &ytext0);
		    vptodevxy_text (xtext1, ytext1, &xtext1, &ytext1);
		    vptodevxy_text (xtext2, ytext2, &xtext2, &ytext2);
		}

		xtext1 -= xtext0;
		xtext2 -= xtext0;
		ytext1 -= ytext0;
		ytext2 -= ytext0;

		savefat = fat;
		savefatmult = fatmult;
		fatmult *= txscale;

		fat = (ifat >= 0)? fatmult * (ifat + fatbase): -1;

		update_color (dev, attributes);

		getvpstring ();

		text (dev, txbuffer,
		      txscale * xtext1 / TEXTVECSCALE,
		      txscale * ytext1 / TEXTVECSCALE,
		      txscale * xtext2 / TEXTVECSCALE,
		      txscale * ytext2 / TEXTVECSCALE);
		
		fat = savefat;
		fatmult = savefatmult;
		break;

	    case VP_OLDAREA:
		npts = vp_getint();
		afat = vp_getint ();

		if (afat >= 0) afat = fatmult * (fatbase + afat);
	    
		nx_temp = vp_getint ();
		ny_temp = vp_getint ();

		nx = nx_temp * patternmult;
		if (nx_temp == 1 || (nx_temp > 0 && nx == 0)) nx = 1;

		ny = ny_temp * patternmult;
		if (ny_temp == 1 || (ny_temp > 0 && ny == 0)) ny = 1;

		nx_orig = nx;
		ny_orig = ny;

                /* If a monochrome device, then fill with dot pattern.
		 * If a color device, fill with solid color. */

		if (!mono) {
		    if (nx_temp == 0 || ny_temp == 0) {
			nx = 0;
			ny = 0;
		    } else {
			nx = 1;
			ny = 1;
		    }
		}

		/* Create a temporary pattern */
		ipat = 0;
		if (nx * ny > 0) {
		    ptr = sf_intalloc (nx*ny);
		    for (i=0; i < nx*ny-1; i++) {
			ptr[i]=0;
		    }
		    ptr[nx*ny-1] = cur_color;
		    
		    pat[ipat].patbits = ptr;
		} else {
		    pat[ipat].patbits = NULL;
		}

		pat[ipat].xdim = nx;
		pat[ipat].ydim = ny;
		pat[ipat].xdim_orig = nx_orig;
		pat[ipat].ydim_orig = ny_orig;

		npts = getpolygon (npts);
		update_color (dev, attributes);

                /* Do the polygon */
		if (npts > 2 && shade && nx > 0 && ny > 0)
		    area (dev, npts, vxbuffer);

                /* And then its border */
		if (afat >= 0) vecoutline (dev,vxbuffer,vector);

		if (nx * ny > 0) free (ptr);
		break;
	    case VP_AREA:	
		npts = vp_getint();
		ipat = pat_color;

		if (NULL == pat[ipat].patbits && shade) {
		    /* Create default pattern (solid fill with this color) */ 
		    nx = 1;
		    ny = 1;

		    ptr = sf_intalloc (nx*ny);
		    for (i=0; i < nx*ny-1; i++) {
			ptr[i]=0;
		    }
		    ptr[nx*ny-1] = cur_color;

		    pat[ipat].patbits = ptr;
		    pat[ipat].xdim = nx;
		    pat[ipat].ydim = ny;
		    pat[ipat].xdim_orig = nx;
		    pat[ipat].ydim_orig = ny;
		}

		npts = getpolygon (npts);

		if (!shade) {
		    /* At least draw the boundary to show where it is */
		    afat = 0;
		    update_color (dev, attributes);
		    vecoutline (dev,vxbuffer,vector);
		} else {
		    /* See whether raster or hatch area */
		    if (pat[ipat].xdim >= 0) { /* raster */
			if (npts > 2 && 
			    pat[ipat].ydim > 0 && 
			    pat[ipat].xdim > 0) {
			    update_color (dev, attributes);
			    area (dev, npts, vxbuffer);
			}
		    } else { /* hatch */
			numhatch = -pat[ipat].xdim;
			angle = SF_PI * pat[ipat].ydim / 180.;
			ptr = pat[ipat].patbits;
			for (i = 0; i < numhatch * 2; i++) {
			    hafat[i] = (*ptr++);
			    hacol[i] = (*ptr++);
			    haoff[i] = (*ptr++);
			    hasiz[i] = (*ptr++);
			}
			if (npts > 2) {
			    /* Hatch patterns don't rotate. */
			    /* The polygon does, the fill pattern doesn't. 
			       genhatch (npts, numhatch, angle, hafat, 
			       hacol, haoff, hasiz, vxbuffer);
			    */
			}
		    }
		}
		break;
	    case VP_FAT:
		ifat = vp_getint ();

		fat = (ifat >= 0)? fatmult * (ifat + fatbase): -1;

		attributes (dev, NEW_FAT, fat, 0, 0, 0);
		break;

	    case VP_COLOR:
		pat_color = vp_getint();

		if (pat_color > MAX_COL || pat_color < 0) {
		    sf_warning("%s: bad color number %d (max %d, min 0)",
			       __FILE__,pat_color, MAX_COL);
		    pat_color = VP_DEFAULT_COLOR;
		}

		next_color = color_set[pat_color][MAP];
		if (next_color != cur_color) {
		    cur_color = next_color;
		    need_devcolor = true;
		}

		pat_color++;
		break;

	    case VP_SET_COLOR_TABLE:	
		/* The logic here is a bit tricky. 
		 * If the device actually has a settable color of that number,
		 * then reset the device's color as asked. Otherwise, take the
		 * closest color that we have and map to that. Color 0 is the
		 * background color, and is special. Merely "close" colors are
		 * not mapped to Color 0... it has to be exact. Devices with
		 * NO settable colors are still sent colors 0 through 7, and
		 * these colors are then always left set to their original
		 * defaults. In this way you can handle terminals like the GIGI
		 * which have colors but not SETTABLE colors. Remember that
		 * the device itself will have to permute colors 0 through 7
		 * to correspond with Vplot's definitions! */

		col_tab_no = vp_getint();
		if (col_tab_no > MAX_COL || col_tab_no < 0)
		    sf_error("%s: Bad color table value %d (%d max)",
			     __FILE__, col_tab_no, MAX_COL);

		red = vp_getint();
		green = vp_getint();
		blue = vp_getint();

		if (red > MAX_GUN || 
		    green > MAX_GUN || 
		    blue > MAX_GUN || 
		    red < 0 || 
		    green < 0 || 
		    blue < 0) 
		    sf_error("%s: Bad color level in color %d (%d,%d,%d),"
			     " %d max",__FILE__,col_tab_no, 
			     red, green, blue, MAX_GUN);

		red_f = red / (float) MAX_GUN;
		green_f = green / (float) MAX_GUN;
		blue_f = blue / (float) MAX_GUN;

		red_fm = redmap[3] + 
		    redmap[0] * red_f + 
		    redmap[1] * green_f + 
		    redmap[2] * blue_f;

		if      (red_fm < 0.) red_fm = 0.;
		else if (red_fm > 1.) red_fm = 1.;
		red_fm = powf (red_fm, redpow);
		red = red_fm * MAX_GUN + 0.5;

		green_fm = greenmap[3] + 
		    greenmap[0] * red_f + 
		    greenmap[1] * green_f + 
		    greenmap[2] * blue_f;

		if      (green_fm < 0.) green_fm = 0.;
		else if (green_fm > 1.) green_fm = 1.;
		green_fm = powf (green_fm, greenpow);
		green = green_fm * MAX_GUN + 0.5;

		blue_fm = bluemap[3] + 
		    bluemap[0] * red_f + 
		    bluemap[1] * green_f + 
		    bluemap[2] * blue_f;

		if      (blue_fm < 0.) blue_fm = 0.;
		else if (blue_fm > 1.) blue_fm = 1.;
		blue_fm = powf (blue_fm, bluepow);
		blue = blue_fm * MAX_GUN + 0.5;

		if (col_tab_no != 0) {
		    red_mask = red;
		    green_mask = green;
		    blue_mask = blue;

                    /* Background color */
		    red_0   = color_set[0][RED];
		    green_0 = color_set[0][GREEN];
		    blue_0  = color_set[0][BLUE];

                    /* Opposite of background color */
		    red_7   = MAX_GUN - color_set[0][RED];
		    green_7 = MAX_GUN - color_set[0][GREEN];
		    blue_7  = MAX_GUN - color_set[0][BLUE];

		    if (red == red_7 && green == green_7 && blue == blue_7) {
			if (!colormask[3]) {
			    red_mask = red_0;
			    green_mask = green_0;
			    blue_mask = blue_0;
			}
		    } else {
			if (colormask[4]) {
			    if (!colormask[0]) red_mask = red_0;
			    if (!colormask[1]) green_mask = green_0;
			    if (!colormask[2]) blue_mask = blue_0;
			} else {
			    if (colormask[0]) {
				green_mask = red_mask;
				blue_mask = red_mask;
			    }
			    if (colormask[1]) {
				red_mask = green_mask;
				blue_mask = green_mask;
			    }
			    if (colormask[2]) {
				red_mask = blue_mask;
				green_mask = blue_mask;
			    }
			}
		    }

		    red = red_mask;
		    green = green_mask;
		    blue = blue_mask;
		}

		/* Process the new color table value */
		if (col_tab_no < num_col) {
		    attributes (dev, SET_COLOR_TABLE, 
				col_tab_no, red, green, blue);
		    color_set[col_tab_no][STATUS] = SET;
		} else if (col_tab_no > 7) {
		    color_set[col_tab_no][STATUS] = MAPPED;
		} else if (col_tab_no > 0) { 
		    color_set[col_tab_no][STATUS] = MAP_SET;
		    /* Means that we need to map these but they are still set
		     * to the original colors and can't be changed. Color 0
		     * is always background, and any other color exactly the
		     * same color as color 0 "wants to be", even if color 0
		     * actually can't be and hasn't been reset, is mapped to
		     * color zero if it can't be set. */	
		}

		/* Save the color that this color table number wants to be */
		color_set[col_tab_no][RED] = red;
		color_set[col_tab_no][GREEN] = green;
		color_set[col_tab_no][BLUE] = blue;

		/* grey level is "Black and White TV style" mapped from color,
		 * with corrections added by the subroutine greycorr. */

		grey = (blue * 1 + red * 2 + green * 4 + 6) / 7;
		color_set[col_tab_no][GREY] = greycorr (grey);

		/* If the next command is also a set color table command 
		 * postpone doing this and kill 2 birds with one stone. */
	     
		c = getchar ();
		skip = true;

		if (EOF == c) break;

		if (VP_SET_COLOR_TABLE != c) {
		    if (mono) {
			for (ii = 1; ii <= MAX_COL; ii++) {
			    color_set[ii][MAP] = 
				(color_set[ii][RED] == color_set[0][RED] &&
				color_set[ii][GREEN] == color_set[0][GREEN] &&
				color_set[ii][BLUE] == color_set[0][BLUE])?
				0:7;
			}
		    } else {
			/* For all color table entries that aren't set, 
			 * find the color that IS set and use that instead. */
			for (ii = num_col; ii <= MAX_COL; ii++) {
			    if (color_set[ii][STATUS] & MAPPED) {
				min_dist = MAX_GUN * MAX_GUN * 8;
				for (i = num_col_8 - 1; i >= 0; i--) {
				    /* Colors 1 through 7 are guaranteed SET, 
				     * we always get an answer. Color zero is
				     * background and special. To map to it you
				     * have to hit its color exactly, and no
				     * other color matched exactly first. */
				    if (color_set[i][STATUS] & SET) {
					if (color_set[i][STATUS] == SET) {
					    k = color_set[i][RED] - 
						color_set[ii][RED];
					    dist = 2 * k * k;
					    k = color_set[i][GREEN] - 
						color_set[ii][GREEN];
					    dist += 4 * k * k;
					    k = color_set[i][BLUE] - 
						color_set[ii][BLUE];
					    dist += k * k;
					} else {
					    k = MAX_GUN * ((i & 2) / 2) - 
						color_set[ii][RED];
					    dist = 2 * k * k;
					    k = MAX_GUN * ((i & 4) / 4) - 
						color_set[ii][GREEN];
					    dist += 4 * k * k;
					    k = MAX_GUN * ((i & 1) / 1) - 
						color_set[ii][BLUE];
					    dist += k * k;
					}
					if (dist < min_dist && 
					    (i != 0 || dist == 0)) {
					    min_dist = dist;
					    best_col = i;
					    if (0 == dist) break;
					}
				    }
				}
				color_set[ii][MAP] = best_col;
			    }
			}
		    }
		}
		break;

	    case VP_PURGE:
		close (dev, CLOSE_FLUSH);
		break;
	    case VP_BREAK:		/* break */
		if (BREAK_IGNORE == brake) break;
		/* fall through */

	    case VP_ERASE: 
		close (dev, CLOSE_FLUSH);

/*
 * Erases and breaks can't occur inside groups.
 
 group_number--;
 if (group_number != 0)
 {
 ERR (WARN, name,
 "group contains erase or break");
		group_number = 0;
		}
		else
		{
		if (erase & DO_LITERALS)
		{
		if (c == VP_ERASE || brake)
		ii = dev.attributes (END_GROUP, group_number, ERASE_EGROUP, 0, 0);
		else
		ii = dev.attributes (END_GROUP, group_number, BREAK_EGROUP, 0, 0);
		}
		else
		ii = dev.attributes (END_GROUP, group_number, IGNORE_EGROUP, 0, 0);
		
		if (ii == DOVPLOT_EXIT)
		return;
		if (ii != DOVPLOT_CONT)
		    ERR (WARN, name, "dev.attributes(END_GROUP,...) returned junk.");
		    }
		    
		    if (epause < 0)
		    {
		    message (MESG_ERASE);
		    message (MESG_ON);
		    message (MESG_HOME);
		    message (MESG_READY);
		    message (MESG_HIGHLIGHT_ON);
		    message (MESG_TEXT, "Type Return to Continue...  ");
		    message (MESG_DONE);
		    ii = dev.interact (INT_PAUSE, controltty, string);
		if (ii == DOVPLOT_EXIT)
		return;
		message (MESG_HIGHLIGHT_OFF);
		if (!allowecho)
		{
		message (MESG_READY);
		message (MESG_TEXT, CRLF);
		}
		message (MESG_DONE);
		message (MESG_OFF);
		message (MESG_ERASE);
		}
		else
	    {
	    if (epause > 0)
	    sleep ((unsigned) epause);
	    }
	    
	     * Inquire point back from device
	     
	     if (interact[0] != '\0')
	     {
		getapoint ();
		}
		
*/

		if (erase & DO_LITERALS) {
		    if ((VP_ERASE == c) || brake) {
			eras (dev, ERASE_MIDDLE);
		    } else {
			eras (dev, ERASE_BREAK);
		    }
		    nplots++;
		}

		new_style = default_style;
		setstyle (new_style,dev,attributes);
		reset_all (dev,attributes);

	    /*
	     * Start a new group level 0 to contain the next frame (separated
	     * from the previous one by either an erase or break)
	     
	    ii = ftell (pltin);
	    sprintf (group_name, "%s.%d", pltname, nplots);
	    dev.attributes (BEGIN_GROUP, group_number, ii, 0, 0);
	    group_number++;
	    */

		if (framewindows) outline_window (dev, attributes, vector);
		break;

	    case VP_WINDOW:
		if (window) {
		    dev->xwmin = vp_getint();
		    dev->ywmin = vp_getint();
		    dev->xwmax = vp_getint();
		    dev->ywmax = vp_getint();

		    if (dev->xwmin > dev->xwmax || 
			dev->ywmin > dev->ywmax) {
			dev->xwmin = xWmax + 1;
			dev->xwmax = xWmin - 1;
			dev->ywmin = yWmax + 1;
			dev->ywmax = yWmin - 1;
		    } else {
			vptodevw (dev->xwmin, dev->ywmin, 
				  dev->xwmax, dev->ywmax, 
				  &(dev->xwmin), &(dev->ywmin), 
				  &(dev->xwmax), &(dev->ywmax));
			if (dev->xwmin < xWmin) dev->xwmin=xWmin;
			if (dev->xwmax > xWmax) dev->xwmax=xWmax;
			if (dev->ywmin < yWmin) dev->ywmin=yWmin;
			if (dev->ywmax > yWmax) dev->ywmax=yWmax;
		    }
		} else {
		    vp_getint();
		    vp_getint();
		    vp_getint();
		    vp_getint();
		}

		reset_windows (dev, attributes);
		if (framewindows) outline_window (dev, attributes, vector);
		break;

	    case VP_NOOP:
		break;

	    case VP_TXALIGN:
		txalign.hor = vp_getint();
		txalign.ver = vp_getint();
		attributes (dev, NEW_ALIGN, txalign.hor, txalign.ver, 0, 0);
		break;

	    case VP_TXFONTPREC:	
		ii = vp_getint();
		if (ii >= 0) {
		    txfont = ii;
		} else {
		    ii = -1;
		}

		jj = vp_getint();
		if (jj >= 0) {
		    txprec = jj;
		} else {
		    jj = -1;
		}

		kk = vp_getint();
		if (kk >= 0) {
		    txovly = kk;
		} else {
		    kk = -1;
		}

		attributes (dev, NEW_FONT, ii, jj, kk, 0);
		break;

	    case VP_OVERLAY:
		overlay = vp_getint();
		attributes (dev, NEW_OVERLAY, overlay, 0, 0, 0);
		break;

	    case VP_PATLOAD:	
		nmul = vp_getint();
		ny = vp_getint();

		/* See whether Raster or Hatch */
		if (ny >= 0) { /* Raster */
		    ny_mult = ny * (ppi / nmul) * patternmult / aspect_ratio;
		    if (ny_mult == 0 && ny > 0) ny_mult = 1;

		    nx = vp_getint();
		    nx_mult = nx * (ppi / nmul) * patternmult;
		    if (nx_mult == 0 && nx > 0) nx_mult = 1;

		    ipat = vp_getint() + 1;
		    if (ipat > NPAT || ipat < 1)
			sf_error ("%s: Bad pattern number %d (max %d, min 1)", 
				  __FILE__,ipat, NPAT);

		    if (NULL != pat[ipat].patbits)
			free (pat[ipat].patbits);

		    if (nx_mult * ny_mult > 0) {
			pat[ipat].patbits = sf_intalloc(nx_mult * ny_mult);
		    } else {
			pat[ipat].patbits = sf_intalloc(1);
			pat[ipat].patbits[0] = 0;
		    }

		    tempbuf = (nx * ny > 0)? sf_intalloc (nx * ny): NULL;

		    /* read in pattern */
		    for (j = 0; j < nx * ny; j++) {
			k = vp_getint();
			if (k > MAX_COL || k < 0) {
			    sf_warning("%s: bad color number in pattern "
				       "%d (max %d, min 0)",
				       __FILE__,k,MAX_COL);
			    k = VP_DEFAULT_COLOR;
			}
			tempbuf[j] = color_set[k][MAP];
		    }

		    /* copy into patbits, with "stretching" */
		    for (i = 0; i < ny_mult; i++) {
			ny_temp = i * ny / ny_mult;
			for (j = 0; j < nx_mult; j++) {
			    nx_temp = j * nx / nx_mult;
			    pat[ipat].patbits[j + nx_mult * i] = 
				tempbuf[nx_temp + nx * ny_temp];
			}
		    }

		    /* set dimensions of pattern */
		    pat[ipat].xdim = nx_mult;
		    pat[ipat].ydim = ny_mult;
		    pat[ipat].xdim_orig = nx_mult;
		    pat[ipat].ydim_orig = ny_mult;

		    if (NULL != tempbuf) free (tempbuf);

		    attributes (dev, NEW_PAT, ipat, 0, 0, 0);
		} else { /* Hatch Pattern */
		    nx = vp_getint();
		    if (nx <= 0 || nx * 2 > NHATCH)
			sf_error ("%s: Bad numhatch %d (max %d/2, min 1)", 
				  __FILE__,nx, NHATCH);

		    ipat = vp_getint() + 1;
		    if (ipat > NPAT || ipat < 1)
			sf_error("%s: Bad pattern number %d (max %d, min 1)",
				 __FILE__,ipat, NPAT);

		    if (NULL != pat[ipat].patbits) free (pat[ipat].patbits);

		    pat[ipat].patbits = sf_intalloc(nx*8);

		    for (i = 0; i < nx * 2; i++) {
			hafat[i] = vp_getint();
			if (hafat[i] >= 0)
			    hafat[i] = fatmult * (fatbase + hafat[i]);

			k = vp_getint();
			if (k > MAX_COL || k < 0) {
			    sf_warning("%s: Bad color number in "
				       "hatch %d (max %d, min 0)",
				       __FILE__,k, MAX_COL);
			    k = VP_DEFAULT_COLOR;
			}

			hacol[i] = color_set[k][MAP];
			haoff[i] = vp_getint() * patternmult * ppi / RPERIN;
			hasiz[i] = vp_getint();
		    }

		    /* Find smallest hatch interval. 1/2 that, and then force
		     * everything to be a multiple of that. If this were not
		     * done, then hatch intervals that originally differed by
		     * some simple fraction might end up slightly off, causing
		     * very unsightly beats. */

		    /* Upper quantization limit of 1/10 inch */
		    k = (RPERIN / 10) * 2;
		    for (i = 0; i < nx * 2; i++) {
			if (hasiz[i] > 0 && hasiz[i] < k)
			    k = hasiz[i];
		    }
		    j = k * (patternmult * ppi / RPERIN) / 2.;
		    if (j < 1) j = 1;

		    for (i = 0; i < nx * 2; i++) {
			hasiz[i] = ((hasiz[i] * 2) / k) * j;
		    }

		    ptr = pat[ipat].patbits;
		    for (i = 0; i < nx * 2; i++) {
			if (haoff[i] >= hasiz[i]) haoff[i] = 0;
			*ptr++ = hafat[i];
			*ptr++ = hacol[i];
			*ptr++ = haoff[i];
			*ptr++ = hasiz[i];
		    }

		    pat[ipat].xdim = -nx;
		    pat[ipat].ydim = nmul;	/* angle */
		}
		break;

	    case VP_BIT_RASTER:	  
	    case VP_BYTE_RASTER:  
		ras_orient = vp_getint();
		if (rotate % 90 != 0) {
		    if (wantras) {
			sf_warning ("%s: Raster only possible in "
				    "4 principal orientations",__FILE__);
			wantras = false;
		    }
		} else {
		    ras_orient += rotate / 90;
		    if (ras_orient >= 0) {
			ras_orient %= 4;
		    } else {
			ras_orient = ((ras_orient % 4) + 4) % 4;
		    }
		}

		ras_offset = vp_getint();

		if (ras_offset > MAX_COL || ras_offset + 255 < 0)
		    sf_error("%s: Absurd raster offset %d",
			     __FILE__,ras_offset);
		
		xvr_min = vp_getint();
		yvr_min = vp_getint();
		xvr_max = vp_getint();
		yvr_max = vp_getint();
		vptodevw (xvr_min, yvr_min, xvr_max, yvr_max,
			  &xvr_min, &yvr_min, &xvr_max, &yvr_max);
		xvru_min = xvr_min;
		yvru_min = yvr_min;
		xvru_max = xvr_max;
		yvru_max = yvr_max;

		xpix = vp_getint();
		ypix = vp_getint();

	    switch (ras_orient) {
		case 0:
		    xrasmult = xpix / (float) (xvr_max - xvr_min);
		    yrasmult = ypix / (float) (yvr_max - yvr_min);
		    xvr_max--;
		    yvr_max--;
		    break;
		case 1:
		    yrasmult = ypix / (float) (xvr_max - xvr_min);
		    xrasmult = xpix / (float) (yvr_max - yvr_min);
		    xvr_max--;
		    yvr_min++;
		    break;
		case 2:
		    xrasmult = xpix / (float) (xvr_max - xvr_min);
		    yrasmult = ypix / (float) (yvr_max - yvr_min);
		    xvr_min++;
		    yvr_min++;
		    break;
		case 3:
		    yrasmult = ypix / (float) (xvr_max - xvr_min);
		    xrasmult = xpix / (float) (yvr_max - yvr_min);
		    xvr_min++;
		    yvr_max--;
		    break;
	    }

	    if (wantras && smart_raster) {
/*
		rasterline = (unsigned char *) malloc ((unsigned) ((xpix + 7 + 2) * sizeof (unsigned char)));
		rasterline2 = (unsigned char *) malloc ((unsigned) ((xpix + 7 + 2) * sizeof (unsigned char)));
		if (rasterline == NULL || rasterline2 == NULL)
		    ERR (FATAL, name, "cannot alloc memory to load raster line");

		outraster = (unsigned char *) malloc ((unsigned) (xpix * ypix) * sizeof (unsigned char));
		if (outraster == NULL)
		    ERR (FATAL, name, "cannot alloc memory to load raster image");


 * See whether we want to dither or not.
 * Dither if monochrome device and dithering has been asked for.
 * It's up to the device to decide whether to actually do this.

		dither_it = dither && mono;


 * If the device is color, but the raster happens to all be shades
 * of grey, let the device know so it can take advantage of the
 * fact. This is available to the routine via the external variable
 * "ras_allgrey".	joe & david 10/03/94
 
              ras_allgrey = YES;
              if (!mono)
              {
                  if (c == VP_BIT_RASTER)
                  {
                      for (ii=0; ii<2; ii++)
                      {
                          if (
                          (color_set [ras_offset * ii][_RED] !=
                          color_set [ras_offset * ii][_GREEN]) ||
                          (color_set [ras_offset * ii][_GREEN] !=
                          color_set [ras_offset * ii][_BLUE])
                          )
                          {
                              ras_allgrey = NO;
                              break;
                          }
                      }
                  }
                  else
                  {
                      for (ii=0; ii<256; ii++)
                      {
                          if (
                          (color_set [ras_offset + ii][_RED] !=
                          color_set [ras_offset + ii][_GREEN]) ||
                          (color_set [ras_offset + ii][_GREEN] !=
                          color_set [ras_offset + ii][_BLUE])
                          )
                          {
                              ras_allgrey = NO;
                              break;
                          }
                      }
                  }
              }

		
		 * Read in the Raster data for "Smart" devices, ie, those
		 * which can stretch (and dither) their own raster.
		 
		num_rep = 0;
		for (yrast = 0; yrast < ypix; yrast++)
		{
		    
		     * Read in the next raster line, if we have a new one
		     
		    if (num_rep <= 0)
		    {
			num_rep = vp_getint();
			if (num_rep <= 0)
			    ERR (FATAL, name, "Bad Raster line multiplier");
			pos = 0;
		new_pat:num_pat = vp_getint();
			num_byte = vp_getint();
			if (num_pat <= 0 || num_byte <= 0 ||
			    pos + num_pat * num_byte > xpix)
			    ERR (FATAL, name, "Raster line not length promised");

			if (num_pat > 1)
			{
			    if (dither_it)
			    {
				READ_RASTER (
					     rasterline2[j] = GREY_MAP (ras_offset + (int) fgetc (pltin)),
					     rasterline2[j + jj] = GREY_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
				 );
			    }
			    else
			    {
				READ_RASTER (
					     rasterline2[j] = COLOR_MAP (ras_offset + (int) fgetc (pltin)),
					     rasterline2[j + jj] = COLOR_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
				 );
			    }
			    for (j = 0; j < num_pat; j++)
			    {
				for (k = 0; k < num_byte; k++)
				{
				    rasterline[pos] = rasterline2[k];
				    pos++;
				}
			    }
			}
			else
			{
			    if (dither_it)
			    {
				READ_RASTER (
					     rasterline[pos + j] = GREY_MAP (ras_offset + (int) fgetc (pltin)),
					     rasterline[pos + j + jj] = GREY_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
				 );
			    }
			    else
			    {
				READ_RASTER (
					     rasterline[pos + j] = COLOR_MAP (ras_offset + (int) fgetc (pltin)),
					     rasterline[pos + j + jj] = COLOR_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
				 );
			    }
			    pos += num_byte;
			}

			if (pos < xpix)
			    goto new_pat;
			if (pos != xpix)
			    ERR (FATAL, name, "Raster line not length promised");
		    }
		    num_rep--;
		    for (ii = 0; ii < xpix; ii++)
		    {
			outraster[ii + yrast * xpix] = rasterline[ii];
		    }
		}

		Smart form 
		dev.raster (xpix, ypix, xvr_min, yvr_min, xvr_max, yvr_max,
			    outraster, ras_orient, dither_it);

		free ((char *) rasterline);
		free ((char *) rasterline2);
		free ((char *) outraster);
*/
	    } else {
		if (xvr_min < dev->xwmin) xvr_min=dev->xwmin;
		if (xvr_max > dev->xwmax) xvr_max=dev->xwmax;
		if (yvr_min < dev->ywmin) yvr_min=dev->ywmin;
		if (yvr_max > dev->ywmax) yvr_max=dev->ywmax;
		
		switch (ras_orient) {
		    case 0:
			xvr_max++;
			yvr_max++;
			xr_max = xvr_max - xvru_min;
			xr_min = xvr_min - xvru_min;
			yr_max = yvr_max - yvru_min;
			yr_min = yvr_min - yvru_min;
			yru_max = yvru_max - yvru_min;
			break;
		    case 1:
			xvr_max++;
			yvr_min--;
			xr_max = yvru_max - yvr_min;
			xr_min = yvru_max - yvr_max;
			yr_max = xvr_max - xvru_min;
			yr_min = xvr_min - xvru_min;
			yru_max = xvru_max - xvru_min;
			break;
		    case 2:
			xvr_min--;
			yvr_min--;
			xr_max = xvru_max - xvr_min;
			xr_min = xvru_max - xvr_max;
			yr_max = yvru_max - yvr_min;
			yr_min = yvru_max - yvr_max;
			yru_max = yvru_max - yvru_min;
			break;
		    case 3:
			xvr_min--;
			yvr_max++;
			xr_max = yvr_max - yvru_min;
			xr_min = yvr_min - yvru_min;
			yr_max = xvru_max - xvr_min;
			yr_min = xvru_max - xvr_max;
			yru_max = xvru_max - xvru_min;
			break;
		}
		xru_min = 0;
		
		if (yr_max < yr_min || 
		    xr_max < xr_min || !wantras) {
		    yr_max = yr_min;
		    xr_max = xr_min;
		}

		rasterline = (unsigned char *) 
		    sf_alloc (xpix + 7 + 2,sizeof (unsigned char));
		rasterline2 = (unsigned char *) 
		    sf_alloc (xpix + 7 + 2,sizeof (unsigned char));

		dither_it = dither && mono;

		if (xr_max > xr_min) {
		    outraster2 = (unsigned char *) 
			sf_alloc (xr_max - xr_min, sizeof (unsigned char));
		    outraster = dither_it? (unsigned char *) 
			sf_alloc (xr_max - xr_min, sizeof (unsigned char)):
			outraster2;
		} else {
		    outraster2 = NULL;
		    outraster = NULL;
		}

		/* Read in the Raster data */
		lastrast = -1;
		num_rep = 0;
		for (i = yr_max - 1; i >= yr_min - 1; i--) {
		    yrast = (i == yr_min - 1)? 
			ypix-1: (yru_max - 1 - i) * yrasmult;

		    for (ii = 0; ii < (yrast - lastrast); ii++) {
			if (num_rep <= 0) {
			    num_rep = vp_getint();
			    if (num_rep <= 0)
				sf_error ("%s: Bad Raster line multiplier",
					  __FILE__);
			    pos = 0;
		    
			    do {
				num_pat = vp_getint();
				num_byte = vp_getint();
				if (num_pat <= 0 || num_byte <= 0 ||
				    pos + num_pat * num_byte > xpix)
				    sf_error ("%s: Raster not length promised",
					      __FILE__);
				
				if (ii + num_rep >= yrast - lastrast
				    && xr_max > xr_min && i >= yr_min) {
				    if (num_pat > 1) {
					if (dither_it) {
					    if (c == VP_BYTE_RASTER) {
						for (j=0; j < num_byte; j++) {
						    k = ras_offset + getchar();
						    rasterline2[j] = 
							color_set[k][GREY];
						}
					    } else {
						for (jj=0;jj<num_byte;jj+=8) {
						    ibyte = (unsigned char) 
							getchar();
						    k = ras_offset * 
							((ibyte & 
							  (001 << (7 - j))) 
							 != 0);
						    rasterline2[j + jj] = 
							color_set[k][GREY];
						}
					    }
					} else {
					    if (c == VP_BYTE_RASTER) {
						for (j=0; j < num_byte; j++) {
						    k = ras_offset + getchar();
						    rasterline2[j] = 
							color_set[k][MAP];
						}
					    } else {
						for (jj=0;jj<num_byte;jj+=8) {
						    ibyte = (unsigned char) 
							getchar();
						    k = ras_offset * 
							((ibyte & 
							  (001 << (7 - j)))
							 != 0);
						    rasterline2[j + jj] = 
							color_set[k][MAP];
						}
					    }
					}
					for (j = 0; j < num_pat; j++) {
					    for (k = 0; k < num_byte; k++) {
						rasterline[pos] = 
						    rasterline2[k];
						pos++;
					    }
					}
				    } else {
					if (dither_it) {
					    if (c == VP_BYTE_RASTER) {
						for (j=0; j < num_byte; j++) {
						    k = ras_offset + getchar();
						    rasterline2[pos + j] = 
							color_set[k][GREY];
						}
					    } else {
						for (jj=0;jj<num_byte;jj+=8) {
						    ibyte = (unsigned char) 
							getchar();
						    k = ras_offset * 
							((ibyte & 
							  (001 << (7 - j)))
							 != 0);
						    rasterline2[pos + j + jj] =
							color_set[k][GREY];
						}
					    }
					} else {
					    if (c == VP_BYTE_RASTER) {
						for (j=0; j < num_byte; j++) {
						    k = ras_offset + getchar();
						    rasterline2[pos + j] = 
							color_set[k][MAP];
						}
					    } else {
						for (jj=0;jj<num_byte;jj+=8) {
						    ibyte = (unsigned char) 
							getchar();
						    k = ras_offset * 
							((ibyte & 
							  (001 << (7 - j)))
							 != 0);
						    rasterline2[pos + j + jj] =
							color_set[k][MAP];
						}
					    }
					}
					pos += num_byte;
				    }
				} else {
				    if (c == VP_BYTE_RASTER) {
					for (j = 0; j < num_byte; j++) {
					    getchar ();
					}
				    } else {
					for (j = 0; j < num_byte; j += 8) {
					    getchar ();
					}
				    }
				    pos += num_pat * num_byte;
				}
			    } while (pos < xpix);

			    if (pos != xpix)
				sf_error("%s: Raster not length promised",
					 __FILE__);

			    if (wantras && xr_max > xr_min && i >= yr_min) {
				for (j = xr_min; j < xr_max; j++) {
				    outraster2[j - xr_min] =
					rasterline[(int)((j - xru_min) * 
							 xrasmult)];
				}
			    }
			}
			num_rep--;
		    }
		    lastrast = yrast;

		    if (xr_max > xr_min && i >= yr_min) {
			if (dither_it)
			    dithline (outraster2, outraster, 
				      xr_max - xr_min, yr_max - 1 - i, dither);

			switch (ras_orient) {
			    case 0:
				raster (dev,yr_max - 1 - i, yr_max - yr_min,
					xvr_min, i + yvru_min,
					xr_max - xr_min, ras_orient, 
					outraster, 0, 0);
			    break;
			    case 1:
				raster (dev,yr_max - 1 - i, yr_max - yr_min,
					xvru_min + i, yvr_max,
					xr_max - xr_min, ras_orient, 
					outraster, 0, 0);
				break;
			    case 2:
				raster (dev,yr_max - 1 - i, yr_max - yr_min,
					xvr_max, yvru_max - i,
					xr_max - xr_min, ras_orient, 
					outraster, 0, 0);
				break;
			    case 3:
				raster (dev,yr_max - 1 - i, yr_max - yr_min,
					xvru_max - i, yvr_min,
					xr_max - xr_min, ras_orient, 
					outraster, 0, 0);
				break;
			}
		    }
		}
		free (rasterline);
		free (rasterline2);

		if (NULL != outraster2) {
		    free (outraster2);
		    if (dither_it) free (outraster);
		}

		if (!wantras) {
		    xxx[0] = xvru_min;
		    yyy[0] = yvru_min;
		    xxx[1] = xvru_min;
		    yyy[1] = yvru_max;
		    xxx[2] = xvru_max;
		    yyy[2] = yvru_max;
		    xxx[3] = xvru_max;
		    yyy[3] = yvru_min;
		    drawpolygon (dev,4, xxx, yyy,attributes,area,vector);
		}
	    }	    
	    break;
		
	    case VP_MESSAGE:
		getvpstring ();
/*
		message (MESG_READY);
		message (MESG_MESSAGE);
		message (MESG_TEXT, txbuffer);
		message (MESG_TEXT, CRLF);
		message (MESG_DONE);
*/
		break;
	    case VP_SETDASH:
		npts = vp_getint();
		if (npts > MAXDASH) 
		    sf_error("%s: Too complicated dash line pattern",__FILE__);
		dashon = npts;
		k = 0;
		dashsum = 0.;
		for (ii = 0; ii < npts * 2; ii++) {
		    dashes[ii] = dashscale * vp_getint() / RPERIN;
		    if (dashes[ii] < 0.)
			sf_error("%s: Negative dash distance.",__FILE__);

		    if (dashes[ii] != 0.) k = 1;

		    dashsum += dashes[ii];
		}
		if (!k) dashon = false;
		attributes (dev,NEW_DASH, dashon, 0, 0, 0);
		break;
	    default:	
		sf_error("%s: invalid VPLOT command decimal %d character %c",
			 __FILE__, c, (char) c);
		break;
	}
	if (!skip) c = getchar();
	skip = false;
    }

    close (dev,CLOSE_FLUSH);

/* End the group for this frame 
    group_number--;
    if (group_number != 0)
    {
	ERR (WARN, name,
	     "group left unclosed at end of file");
	group_number = 0;
    }
    else
    {
	if ((erase & FORCE_INITIAL) && (erase & DO_LITERALS))
	    ii = dev.attributes (END_GROUP, group_number, EOF_ERASE_EGROUP, 0, 0);
	else
	    ii = dev.attributes (END_GROUP, group_number, EOF_IGNORE_EGROUP, 0, 0);

	if (ii == DOVPLOT_EXIT)
	    return;
	if (ii != DOVPLOT_CONT)
	    ERR (WARN, name, "dev.attributes(END_GROUP,...) returned junk.");
    }

*/
}

static void reset_windows (vp_device dev, 
			   void (*attributes)(vp_device,
					      vp_attribute,int,int,int,int))
{
    if (dev->xwmax != xwmax_last || 
	dev->ywmax != ywmax_last || 
	dev->xwmin != xwmin_last || 
	dev->ywmin != ywmin_last)
	attributes (dev,SET_WINDOW, 
		    dev->xwmin, dev->ywmin, dev->xwmax, dev->ywmax);

    xwmin_last = dev->xwmin;
    ywmin_last = dev->ywmin;
    xwmax_last = dev->xwmax;
    ywmax_last = dev->ywmax;
}

static void outline_window (vp_device dev, 
			    void (*attributes)(vp_device,
					       vp_attribute,int,int,int,int),
			    void (*vector)(vp_device,int,int,int,int,int,bool))
{
    int color;

    color = color_set[VP_DEFAULT_COLOR][MAP];

    if (need_devcolor || cur_color != color) {
	attributes (dev, SET_COLOR, color, 0, 0, 0);
	need_devcolor = false;
    }

    vector (dev, dev->xwmin, dev->ywmin, dev->xwmax, dev->ywmin, 0, false);
    vector (dev, dev->xwmax, dev->ywmin, dev->xwmax, dev->ywmax, 0, false);
    vector (dev, dev->xwmax, dev->ywmax, dev->xwmin, dev->ywmax, 0, false);
    vector (dev, dev->xwmin, dev->ywmax, dev->xwmin, dev->ywmin, 0, false);

    if (cur_color != color) 
	attributes (dev, SET_COLOR, cur_color, 0, 0, 0);
}

static void update_color (vp_device dev, 
			  void (*attributes)(vp_device,
					     vp_attribute,int,int,int,int))
{
    if (need_devcolor) {
	attributes (dev, SET_COLOR, cur_color, 0, 0, 0);
	need_devcolor = false;
    }
}

static void getxy (int *x, int *y)
{
    int x0, y0;

    x0 = vp_getint();
    y0 = vp_getint();

    vptodevxy(x0,y0,x,y);
}

static void getxy_text (int *x, int *y)
{
    int x0, y0;

    x0 = vp_getint();
    y0 = vp_getint();

    vptodevxy_text(x0,y0,x,y);
}

static void getvpstring (void)
{
    char *txptr;
    int i;

    txptr = txbuffer;

    while (1) {
	i = getchar ();

	if (EOF == i) 
	    sf_error("%s: Unexpected EOF encountered in text string",__FILE__);

	*txptr = i;
	txptr++;

	if (0 == i) break;

	if ((txptr - txbuffer) == txbuflen) {	    
	    txbuffer = (char *) 
		sf_realloc (txbuffer, (txbuflen + TXBUFLEN), sizeof(char));
	    txptr = txbuffer + txbuflen;
	    txbuflen += TXBUFLEN;
	}
    }
}

static int getpolygon (int npts)
{
    int i;
    struct vp_vertex *vertex;

    if (npts > vxbuflen - 1) {
	vxbuflen = npts+1;
	vxbuffer = (struct vp_vertex *) 
	    sf_realloc (vxbuffer,vxbuflen,sizeof(*vxbuffer));
    }
    vertex = vxbuffer;

    getxy (&xnew, &ynew);
    vertex->x = xnew;
    vertex->y = ynew;
    vertex->next = vertex + 1;
    vertex++;
    i = npts - 1;
    while (i--) {
	getxy (&(vertex->x), &(vertex->y));
	if ((vertex->x != xnew) || 
	    (vertex->y != ynew)) {
	    xnew = vertex->x;
	    ynew = vertex->y;
	    vertex->next = vertex + 1;
	    vertex->last = vertex - 1;
	    vertex++;
	    continue;
	}
	npts--;
    }
    vertex--;
    vxbuffer->last = vertex;
    vertex->next = vxbuffer;
    
    return npts;
}

/* Draw the outline of the polygon pointed to by 'head'.  If any side of
 * the polygon is identical to any other (that is they have identical
 * endpoints), then neither line is drawn.  This allows among things a
 * doughnut to be defined by a single polygon without the line that
 * connects the inner and outer being plotted. */

static void vecoutline (vp_device dev, struct vp_vertex *head, 
			void (*vector)(vp_device,int,int,int,int,int,bool))
{ 
    int xlast, ylast;
    struct vp_vertex *v;

    xlast = head->last->x;
    ylast = head->last->y;

    for (v = head; v->next != head; v = v->next) {
	if (!dupside (v))
	    vector (dev,v->x, v->y, xlast, ylast, afat, false);
	xlast = v->x;
	ylast = v->y;
    } 
}

/* Determine if other sides in the polygon are
 * identical to the side specified by the
 * vertices v and v->b. */
static bool dupside (struct vp_vertex *base)
{
    struct vp_vertex *v;
    int x1, x2, y1, y2;

    x1 = base->x;
    x2 = base->last->x;
    y1 = base->y;
    y2 = base->last->y;
    for (v = base->next; v != base; v = v->next) {
	if (x1 == v->x && 
	    y1 == v->y && 
	    x2 == v->last->x && 
	    y2 == v->last->y)
	    return true;

	if (x2 == v->x && 
	    y2 == v->y && 
	    x1 == v->last->x && 
	    y1 == v->last->y)
	    return true;
    }

    return false;
}

/* reset variables that can be affected by vplot commands when
 * processing multiple plots, and don't stay set across pages */
static void reset_all (vp_device dev, 
		       void (*attributes)(vp_device,
					  vp_attribute,int,int,int,int))
{
    int i, j, k, color;

    dev->xwmin = xWmin;		/* plot window parameters defaulted */
    dev->xwmax = xWmax;		/* to maximum size	 */
    dev->ywmin = yWmin;
    dev->ywmax = yWmax;

    reset_windows (dev, attributes);

    fat = fatmult * fatbase;
    attributes (dev, NEW_FAT, fat, 0, 0, 0);

    color =  color_set[VP_DEFAULT_COLOR][MAP];
    if (cur_color != color) {
	need_devcolor = true;
	cur_color = color;
    }

    txalign.hor = TH_NORMAL;
    txalign.ver = TV_NORMAL;
    attributes (dev, NEW_ALIGN, txalign.hor, txalign.ver, 0, 0);

    i = -1;
    j = -1;
    k = -1;
    if (txfont != default_txfont) {
	txfont = default_txfont;
	i = txfont;
    }
    if (txprec != default_txprec) {
	txprec = default_txprec;
	j = txprec;
    }
    if (txovly != default_txovly) {
	txovly = default_txovly;
	k = txovly;
    }
    attributes (dev, NEW_FONT, i, j, k, 0);

    dashon = false;
    attributes (dev, NEW_DASH, dashon, 0, 0, 0);

    overlay = default_overlay;
    attributes (dev, NEW_OVERLAY, overlay, 0, 0, 0);
}

/* 1 random threshold
 * 2 256 element ordered dither (oriented at 0 degrees)
 * 3 Floyd-Steinberg minimized average error method
 * 4 32 element halftone (oriented at 45 degrees) */
static void dithline (unsigned char *inpline, 
		      unsigned char *outline, 
		      int npixels, int linenum, int imethod)
{
    int greydata, i1, ipoint, jpoint, irand;
    float pixel, pixerr, nexterr;
    const int pix_on = 0, pix_off = 7;
    static float *errline = NULL;
    static int ialloc = 0;
    static float    alpha = 0.4375;
    const float beta = 0.1875, gamma = 0.3125, delta = 0.0625;
    const int dith256[256] = {
	1,128,32,160,8,136,40,168,2,130,34,162,10,138,42,170,
	192,64,224,96,200,72,232,104,194,66,226,98,202,74,234,106,
	48,176,16,144,56,184,24,152,50,178,18,146,58,186,26,154,
	240,112,208,80,248,120,216,88,242,114,210,82,250,122,218,90,
	12,140,44,172,4,132,36,164,14,142,46,174,6,134,38,166,
	204,76,236,108,196,68,228,100,206,78,238,110,198,70,230,102,
	60,188,28,156,52,180,20,148,62,190,30,158,54,182,22,150,
	252,124,220,92,244,116,212,84,254,126,222,94,246,118,214,86,
	3,131,35,163,11,139,43,171,1,129,33,161,9,137,41,169,
	195,67,227,99,203,75,235,107,193,65,225,97,201,73,233,105,
	51,179,19,147,59,187,27,155,49,177,17,145,57,185,25,153,
	243,115,211,83,251,123,219,91,241,113,209,81,249,121,217,89,
	15,143,47,175,7,135,39,167,13,141,45,173,5,133,37,165,
	207,79,239,111,199,71,231,103,205,77,237,109,197,69,229,101,
	63,191,31,159,55,183,23,151,61,189,29,157,53,181,21,149,
	254,127,223,95,247,119,215,87,253,125,221,93,245,117,213,85};
    const int halftone32[64] = {
	92,100,124,148,164,156,132,108,
	28,20,76,220,228,236,180,36,
	4,12,84,212,252,244,172,44,
	52,60,116,188,204,196,140,68,
	164,156,132,108,92,100,124,148,
	228,236,180,36,28,20,76,220,
	252,244,172,44,4,12,84,212,
	204,196,140,68,52,60,116,188
    };

    switch (imethod) {
	case 1: /* Random Dither */
	    for (i1 = 0; i1 < npixels; i1++) {
		greydata = inpline[i1];
		irand = (random () & 255);

		outline[i1] = (greydata > irand)? pix_off: pix_on;
	    }
	    break;
	case 2: /* Ordered Dither */
	    for (i1 = 0; i1 < npixels; i1++) {
		greydata = inpline[i1];
		ipoint = i1 % 16;
		jpoint = linenum % 16;
		ipoint = ipoint * 16 + jpoint;

		outline[i1] = (greydata > dith256[ipoint])? pix_off:pix_on;
	    }
	    break;
	case 3: /* Floyd-Steinberg */
	    if (ialloc < npixels) {
		if (ialloc > 0) {
		    free (errline);
		    ialloc = 0;
		}
		errline = sf_floatalloc (npixels);
		ialloc = npixels;
		for (i1 = 0; i1 < npixels; i1++) {
		    errline[i1] = 0.;
		}
	    }

	    nexterr = errline[0];
	    for (i1 = 0; i1 < npixels; i1++) {
		pixel = inpline[i1];
		pixel += nexterr;
		if (pixel < 128) {
		    outline[i1] = pix_on;
		    pixerr = pixel;
		} else {
		    outline[i1] = pix_off;
		    pixerr = pixel - 255;
		}
		if (i1 < npixels - 1) {
		    nexterr = errline[i1 + 1] + pixerr * alpha;
		    errline[i1 + 1] = pixerr * delta;
		}
		if (i1 > 0) {
		    errline[i1 - 1] += pixerr * beta;
		}
		if (i1 == 0) {
		    errline[i1] = pixerr * gamma;
		} else {
		    errline[i1] += pixerr * gamma;
		}
	    }
	    break;
	case 4: /* 32 element halftone at 45 degrees */
	default:
	    for (i1 = 0; i1 < npixels; i1++) {
		greydata = inpline[i1];
		ipoint = i1 % 8;
		jpoint = linenum % 8;
		ipoint = ipoint * 8 + jpoint;

		outline[i1] = (greydata > halftone32[ipoint])? pix_off: pix_on;
	    }
	    break;
    }
}

static void drawpolygon (vp_device dev, int npts, int *x, int *y,
			 void (*attributes)(vp_device,vp_attribute,
					    int,int,int,int),
			 void (*area) (vp_device,int,struct vp_vertex*),
			 void (*vector)(vp_device,int,int,int,int,int,bool))
{
    int i, j;
    static int point;
    struct vp_vertex *vertex;

    j = 0;
    if (npts > vxbuflen - 1)
	vxbuffer = (struct vp_vertex *) 
	    sf_realloc (vxbuffer,npts+1,sizeof(*vxbuffer));
    vertex = vxbuffer;

    xnew = x[j];
    ynew = y[j];
    j++;
    vertex->x = xnew;
    vertex->y = ynew;
    vertex->next = vertex + 1;
    vertex++;

    i = npts - 1;
    while (i--) {
	vertex->x = x[j];
	vertex->y = y[j];
	j++;
	if ((vertex->x != xnew) || 
	    (vertex->y != ynew)) {
	    xnew = vertex->x;
	    ynew = vertex->y;
	    vertex->next = vertex + 1;
	    vertex->last = vertex - 1;
	    vertex++;
	    continue;
	}
	npts--;
    }
    vertex--;
    vxbuffer->last = vertex;
    vertex->next = vxbuffer;

    ipat = 0;
    point = cur_color;
    pat[ipat].patbits = &point;
    pat[ipat].xdim = 1;
    pat[ipat].ydim = 1;
    pat[ipat].xdim_orig = 1;
    pat[ipat].ydim_orig = 1;
    update_color (dev, attributes);
    if (npts > 2) {
	if (shade) {
	    area (dev, npts, vxbuffer);
	} else {
	    vecoutline (dev,vxbuffer,vector);
	}
    }
}

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
			     unsigned char*,int,int))
{
    int i;
    char *arg;

    vplot_init (dev, open);
    vplot_do(dev,reset,close,eras,attributes,vector,marker,text,area,raster);

    for (i=1; i < argc; i++) {
	arg = argv[i];
	if (NULL == strchr(arg,'=')) { /* not a parameter */
	    if (NULL == freopen(arg,"rb",stdin))
		sf_error("%s: Cannot open file %s:",__FILE__,arg);
	    vplot_do(dev,reset,close,eras,attributes,
		     vector,marker,text,area,raster);
	}
    }
}
