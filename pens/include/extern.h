/*
 * defines structures referred to herein
 */
#include <rsf.h>

#include "../include/vertex.h"
#include "../include/device.h"
/*
 * external variables (quite a few of them)
 */

/*
 * these must be DECLARED in dev.conf
 */
extern const char *documentation[];
extern int doclength;
extern char name[];

/*
 * these MUST be initialized in either dev.open or dev.reset
 * (Not setting them will cause a fatal error.)
 */
/*
extern int dev_xmax, dev_ymax, dev_xmin, dev_ymin;
extern float pixels_per_inch;
extern float aspect_ratio;
extern int num_col;
*/

/*
 * options and variables that may also need to be reset in dev.open
 * They can't be reset in dev.reset, because the user may override
 * the values set by the device. (Except that mono will not be reset
 * to NO if the device says it is YES.)
 */
extern bool mono, invras;
extern float fatmult;
extern float patternmult;
extern bool shade;
extern bool wantras;
extern int rotate;
extern float  hshift, vshift;
extern int dither;
extern bool endpause;
extern bool honor_background;
/* extern int txfont,txprec,txovly; */
extern float pixc, greyc;

/*
 * these can also be set in dev.open or dev.reset if dovplot gets them wrong,
 * but can usually be left at their default values.

extern bool need_end_erase;
extern bool smart_clip;
extern bool smart_raster;
extern bool cachepipe;
*/

/*
 * These variables may be useful for referring to in dev.open,
 * so that similar devices can be merged under one pen filter.
 */
extern char callname[];

/*
 * Usual place to read from and place to send the plot to.
 * The device can use these if they are appropriate, or reject
 * these and handle things on its own instead.
 */
extern FILE *pltin;
extern FILE *pltout;

extern struct device dev;
extern struct s_txalign txalign;
extern void (*message)(int command, const char* string);

/*
 * options
 */
extern int epause;
extern int size;
extern int echo;
extern int  xorig, yorig;
extern bool overlay;
extern float  xscale, yscale, txscale, mkscale, dashscale;
extern float  hdevscale, vdevscale;
extern bool serifs_OK;

/*
 * variables
 */
extern int xold,yold;
extern int xnew,ynew;
extern int xwmin,xwmax,ywmin,ywmax;
extern int fat,fatbase,dashon;
extern int afat;
extern int ipat;
extern float mxx, mxy, myx, myy;
extern float dashes[];
extern float dashpos, dashsum;
extern char pltname[];
extern char group_name[];
extern int group_number;
extern int ras_allgrey;		/* david,  10/03/94	*/

extern double   path_orient_dx, path_orient_dy;
extern double   up_orient_dx, up_orient_dy;
