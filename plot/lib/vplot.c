#include <math.h>
#include <stdio.h>

#include <rsf.h>

#include "vplot.h"
#include "coltab.h"

static float fx=0.0, fy=0.0;             /* origin in inches */
static float ufx=0.0, ufy=0.0;           /* origin in user units */
static float xold=0.0, yold=0.0;         /* old pen position (in inches) */
static float dashpos=0.0;                /* position in dashes (in inches) */
static float ddef[4]={0.0,0.0,0.0,0.0};  /* definition of dashes (in inches) */
static float xscl=1.0, yscl=1.0;         /* scaling from user units to inches */
static bool pendown=false;               /* is pen down (or up) */
static bool dashon=false;                /* true if dashed line */
static vp_plotstyle style=STANDARD;          
static int font=-1, prec=-1, ovly=-1;    /* font, precision, and overlay */
static int xjust=0, yjust=0;             /* x and y text justification */

static void pout (float xp, float  yp, int  down);

int vp_getint (void)
{
    unsigned short w;

    w = getchar ();
    w += (getchar () << 8);

    return ((int) ((short) w));
}

void vp_putint (int w)
{
    char j;
    
    j = w & 255;        /* low order byte of halfword value */
    putchar (j);
    
    j = (w >> 8) & 255;	/* high order byte of halfword value */
    putchar (j);
}

void vp_putfloat (float w)
{
    w *= RPERIN;
    vp_putint((int) (w < 0.0)? w-0.5 : w+0.5);
}

void vp_putfloat0 (float w)
{
    vp_putint((int) (w < 0.0)? w-0.5 : w+0.5);
}

void vp_egroup (void)
{
    putchar (VP_END_GROUP);
}

void vp_erase (void)
{
    putchar (VP_ERASE);
}

void vp_fat (int f)
{
    putchar (VP_FAT);
    vp_putint (f);
}

void vp_fill (const float *xp, const float *yp, int  np)
{
    int i;

    putchar (VP_AREA);
    vp_putint (np);

    for (i = 0; i < np; i++) {
	vp_putfloat (xp[i]);
	vp_putfloat (yp[i]);
    }
}

void vp_area (const float *xp, const float *yp, int np, 
	      int fat, int xmask, int ymask)
{
    int i;

    putchar (VP_OLDAREA);
    vp_putint (np);
    vp_putint (fat);
    vp_putint (xmask);
    vp_putint (ymask);
    for (i = 0; i < np; i++) {
	vp_putfloat (xp[i]);
	vp_putfloat (yp[i]);
    }
}

void vp_uarea (const float *xp, const float *yp, int np, 
	       int fat, int xmask, int ymask)
{
    int i;
    float x, y;

    putchar (VP_OLDAREA);
    vp_putint (np);
    vp_putint (fat);
    vp_putint (xmask);
    vp_putint (ymask);
    for (i = 0; i < np; i++) {
	x = fx + (xp[i]-ufx) * xscl;
	y = fy + (yp[i]-ufy) * yscl;
	vp_putfloat (x);
	vp_putfloat (y);
    }
}

void vp_coltab (int color, float r, float g, float b)
{
    r *= MAX_GUN; r += (r < 0.0)? -0.5: +0.5;
    g *= MAX_GUN; g += (g < 0.0)? -0.5: +0.5;
    b *= MAX_GUN; b += (b < 0.0)? -0.5: +0.5;

    putchar (VP_SET_COLOR_TABLE);
    vp_putint (color);
    vp_putint ((int) r);
    vp_putint ((int) g);
    vp_putint ((int) b);
}


void vp_gtext (float x, float y, 
	       float xpath, float ypath, 
	       float xup, float yup, const char *string)
{
    pout (x, y, 0);
    putchar (VP_GTEXT);
    vp_putfloat (TEXTVECSCALE * xpath);
    vp_putfloat (TEXTVECSCALE * ypath);
    vp_putfloat (TEXTVECSCALE * xup);
    vp_putfloat (TEXTVECSCALE * yup);

    do {
	putchar (*string);
    } while (*string++);
}

void vp_hatchload (int angle,int nhatch, int ihatch, int *hatch)
{
    int c, i;

    putchar (VP_PATLOAD);
    vp_putint (angle);
    vp_putint (-1);
    vp_putint (nhatch);
    vp_putint (ihatch);

    for (i = 0; i < 2 * 4 * nhatch; i++) {
	c = hatch[i];
	if (i % 4 > 1) c *= RPERIN / HATCHPERIN;
	vp_putint (c);
    }
}

void vp_message (const char *string)
{
    putchar (VP_MESSAGE);

    do {
	putchar (*string);
    } while (*string++);
}

void vp_move (float x,float  y)
{
    vp_plot (x, y, false);
}

void vp_umove (float x,float  y)
{
    x = fx + (x - ufx) * xscl;
    y = fy + (y - ufy) * yscl;
    vp_plot (x, y, false);
}

void vp_orig (float x,float  y)
{
    fx = x;
    fy = y;
}

void vp_uorig (float x,float  y)
{
    ufx = x;
    ufy = y;
}

void vp_patload (int ppi, int  nx, int ny, int ipat, int *col)
{
    int c, i;

    putchar (VP_PATLOAD);
    c = ppi;
    vp_putint (c);
    c = nx;
    vp_putint (c);
    c = ny;
    vp_putint (c);
    c = ipat;
    vp_putint (c);

    for (i = 0; i < nx * ny; i++) {
	c = col[i];
	vp_putint (c);
    }
}

/* go to location (x,y) and then put the pen down */
void vp_pendn (float x, float y)
{
    vp_plot (x, y, pendown);
    pendown = true;
}

void vp_upendn (float x, float y)
{
    vp_uplot (x, y, pendown);
    pendown = true;
}

void vp_penup (void)
{
    pendown = false;
}

void vp_pline (const float *xp, const float *yp, int np)
{
    int i;

    putchar (VP_PLINE);
    vp_putint (np);
    for (i = 0; i < np; i++) {
	vp_putfloat (xp[i]);
	vp_putfloat (yp[i]);
    }
}

void vp_upline (const float *xp, const float *yp, int np)
{
    int i;
    float x, y;

    putchar (VP_PLINE);
    vp_putint (np);
    for (i = 0; i < np; i++) {
	x = fx + (xp[i]-ufx) * xscl;
	y = fy + (yp[i]-ufy) * yscl;
	vp_putfloat (x);
	vp_putfloat (y);
    }
}

void vp_plot (float x, float y, bool  down)
{
    float dx, dy, dist, dpos, xp, yp, tonext, cosine, sine;
    int i;

/* if(vp_pc._pltout==0){ vp_plot_init();} */

    if (!down || !dashon) { /* if move or no dashes */
	pout (x, y, down);	/* output a move or draw */
	xold = x;
	yold = y;		/* save old x and y */
	dashpos = 0.0;		/* reset position in dashes */
	return;
    }
    
    dx = x - xold;
    dy = y - yold;	/* change in x and y */
    dist = hypotf (dx,dy);	/* distance */
    if (dist == 0.) return;	/* return if no change */
    
    cosine = dx / dist;
    sine = dy / dist;
    dpos = dashpos;		/* current position in dashes */
    dashpos = fmodf (dpos + dist, ddef[3]);	/* next position in
							 * dashes */
    for (i = 0; i < 4 && ddef[i] <= dpos; i++);	/* index to dash def */
    xp = xold;
    yp = yold;		/* initialize xp and yp */
    while (dist > 0.0) {
	tonext = ddef[i] - dpos;	/* dist to next gap or dash */
	if (tonext > dist) tonext = dist;
	xp += tonext * cosine;
	yp += tonext * sine;
	pout (xp, yp, !(i % 2));
	dpos = ddef[i];	/* new position */
	i = (i + 1) % 4;	/* i = 0,1,2, or 3 */
	if (i == 0) dpos = 0.0;		/* back to start of dashes */
	dist -= tonext;
    }
    xold = xp;
    yold = yp;
}

static void pout (float xp, float  yp, int  down)
{
    
    if      (xp >  VP_MAX) xp =  VP_MAX;
    else if (xp < -VP_MAX) xp = -VP_MAX;
    if      (yp >  VP_MAX) yp =  VP_MAX;
    else if (yp < -VP_MAX) yp = -VP_MAX;
    putchar ((down ? VP_DRAW : VP_MOVE));
    vp_putfloat (xp);
    vp_putfloat (yp);
}

void vp_draw (float x,float  y)
{
    vp_plot (x, y, true);
}

void vp_udraw (float x,float  y)
{
    x = fx + (x - ufx) * xscl;
    y = fy + (y - ufy) * yscl;
    vp_plot (x, y, true);
}

void vp_pmark (int npts, int mtype, int msize, const float *xp, const float *yp)
{
    int i;

    putchar (VP_PMARK);
    vp_putint (npts);
    vp_putint (mtype);
    vp_putint (msize);
    
    for (i = 0; i < npts; i++) {
	vp_putfloat (xp[i]);
	vp_putfloat (yp[i]);
    }
}

void vp_purge ()
{
    putchar (VP_PURGE);
    fflush (stdout);
}

/*
 * This routine sets up a "raster color table", to be used with vp_raster
 * with an offset of 256.
 * It sets this color table up in colors 256 through 511.
 * Colors 0 through nreserve-1 are left untouched; you can either let
 * these default or set them yourself.
 * "colname" is a string which defines what raster color table you
 * want to use. For now, it is the same as the "color" option of "Movie",
 * with the extensions that appending a "C" means to flag clipped colors
 * in red, and a string longer than 2 characters is interpreted as a
 * "colfile" name.
 */

void vp_rascoltab (int nreserve, char *colname)
{
    int i, j, k, incr, smap[256];
/*
 * float arrays with numbers in the range 0. through 1.,
 * that define our "raster" color table.
 */
    float red[256], green[256], blue[256];

/*
 * This call makes a color table given its name. We ask for
 * a color table 256 long.
 */
    name2coltab (colname, 256, red, green, blue);

    /*
     * In colors 256 through 512, we define our "raster" color table,
     * unscrambled and with a full range of values. When we call "vp_raster",
     * we'll give an offset of 256 to shift 0-255 into this range. 
     */
    for (i = 0; i < 256; i++)
    {
	vp_coltab (i+256, red[i], green[i], blue[i]);
    }


    /*
     * set the lower "scrambled" scale, which will give us the optimal raster
     * color scale possible for any device, no matter how many color table
     * entries it has. 
     */
    smap[0] = 0;
    smap[1] = 255;

    /*
     * This obscure code does the scrambling. A uniform grey scale, for
     * example, that goes "0, 1, 2, 3, ..., 255" gets reordered to "0, 255,
     * 128, 64, 196, 32, 96, 160, ... 1, 3, ... 251, 253". This way the
     * "most important" colors are closer to the front. This way the first
     * colors lost due to a less-than-256-color device or less than 256
     * colors set aside for raster are the least important ones. 
     */
    for (incr = 2; incr < 512; incr *= 2) {
	for (k = 1, i=2; k < incr && i < 256; k += 2, i++) {
	    j = k * 256 / incr;
	    smap[i] = j;
	}
    }

    for (i = 0; i < 256 - nreserve; i++) {
	j = smap[i];
	vp_coltab (i + nreserve, red[j],green[j],blue[j]);
    }
}

#define ARRAY(A,B)	array[((invert==0)?(A):(ypix-1-(A)))*xpix+(B)]

#define MAXREP	8
#define PUTHSIZE sizeof(int)

#define RASOUT(NUMPAT,NUMBYTE,BITFLAG,BYTES)	\
{\
    vp_putint ((int) (NUMPAT));\
    vp_putint ((int) (NUMBYTE));\
\
    if (BITFLAG)\
    {\
	for (n = 0; n + 7 < (NUMBYTE); n += 8)\
	{\
	    obyte = 0x00;\
\
	    obyte =\
		((*((BYTES) + n + 0) != 0) << 7) |\
		((*((BYTES) + n + 1) != 0) << 6) |\
		((*((BYTES) + n + 2) != 0) << 5) |\
		((*((BYTES) + n + 3) != 0) << 4) |\
		((*((BYTES) + n + 4) != 0) << 3) |\
		((*((BYTES) + n + 5) != 0) << 2) |\
		((*((BYTES) + n + 6) != 0) << 1) |\
		((*((BYTES) + n + 7) != 0) << 0);\
	    putchar ((char) obyte);\
	}\
	if (n < (NUMBYTE))\
	{\
	    obyte = 0x00;\
	    for (nn = 7; n < (NUMBYTE); n++, nn--)\
	    {\
		obyte |= ((*((BYTES) + n) != 0) << nn);\
	    }\
	    putchar ((char) obyte);\
	}\
    }\
    else\
    {\
	for (n = 0; n < (NUMBYTE); n++)\
	{\
	    putchar ((char) * ((BYTES) + n));\
	}\
    }\
}

void vp_raster (unsigned char *array, int blast, int bit, int offset, 
		int xpix, int ypix, 
		float xll, float yll, float ppi,
		float *xur,float *yur, int orient, int invert)
{
    int m, n, nn, l, k, j, i, count, ucount, bitbyte;
    unsigned char obyte;

    if (bit) {
	bitbyte = 8;
	putchar (VP_BIT_RASTER);
    } else {
	bitbyte = 1;
	putchar (VP_BYTE_RASTER);
    }

    if (orient >= 0) {
	orient %= 4;
    } else {
	orient = ((orient % 4) + 4) % 4;
    }

    vp_putint (orient);
    vp_putint (offset);
    vp_putfloat (xll);
    vp_putfloat (yll);
    if (ppi > 0) {
	xll += xpix / ppi;
	yll += ypix / ppi;
	*xur = xll;
	*yur = yll;
    } else {
	xll = *xur;
	yll = *yur;
    }
    vp_putfloat (xll);
    vp_putfloat (yll);

    vp_putint (xpix);
    vp_putint (ypix);

    if (blast) {
	for (i = 0; i < ypix; i++) {
	    vp_putint ((int) 1);
	    RASOUT (1, xpix, bit, &ARRAY (i, 0));
	}
    } else { /* Try to compact it */
	count = 1;
	for (i = 0; i < ypix; i++) {
	    for (j = 0; j < xpix; j++) {
		if (i == ypix - 1 || ARRAY (i, j) != ARRAY (i + 1, j)) {
/* Write this row out, it's not just a repeat of the next one */
/* Write out how many times this line is to be repeated */
		    vp_putint (count);
/*
 * This entire section tries to efficiently represent ONE RASTER LINE
 */
		    
/*
 * Keep track of "unaccounted for" bytes between
 * bytes containted in a pattern
 */
		    ucount = 0;
/* Loop through a raster line */
		    for (k = 0; k < xpix; k++) {
/* Try different size patterns */
			for (l = 1; l <= MAXREP; l++) {
/* See how far this pattern size works */
			    for (m = 0; m <= xpix - (k + l); m++) {
/* Pattern stopped working */
				if (m == xpix - (k + l) ||
				    ARRAY (i, k + m) != ARRAY (i, k + m + l)) {
/* See how many repetitions we got */
				    m = l * (1 + m / l);
/* Is it worth it to use this? */
/* (PUTHSIZE is the number of bytes needed for a puth */
				    if (PUTHSIZE * 2 + l <= 
					1 + (m - 1) / bitbyte) {
/* First take care of bytes not in a pattern */
					if (ucount > 0) {
					    RASOUT (1, ucount, bit, 
						    &ARRAY (i, k - ucount));
					    ucount = 0;
					}
/* Now take care of the bytes in the pattern */
					RASOUT (m / l, l, bit, &ARRAY (i, k));
					k += m - 1;
/* Stop looking, we already found what we wanted */
					goto compact_found;
				    }
/* This pattern didn't work. Stop looking at this size and try a longer one */
				    break;
				}
			    }
			}
/* No pattern including this byte, add it to the unaccounted for list */
			ucount++;
		    compact_found: 
			continue;
		    }
/* Take care of any bytes left hanging on the end */
		    if (ucount > 0) {
			RASOUT (1, ucount, bit, &ARRAY (i, k - ucount));
		    }
/*
 * END of section that does ONE RASTER LINE
 */


/* Reset counter */
		    count = 1;
/* Exit the loop */
		    break;
		}
	    }
/* If we didn't write it out, it will just be a repeat of the next one */
	    if (j == xpix) count++;
	}
    }
}

void vp_scale (float xscale, float  yscale)
{
    xscl = xscale;
    yscl = yscale;
}

void vp_setdash (const float *dash, const float *gapp, int np)
{
    int i;

    putchar (VP_SETDASH);
    vp_putint (np);
    for (i = 0; i < np; i++) {
	vp_putfloat (dash[i]);
	vp_putfloat (gapp[i]);
    }
}

void vp_stretch (float xmin, float ymin, float xmax, float ymax) {
    vp_uorig (xmin, ymin);
    vp_orig (0., 0.);

    if (ROTATED == style) {
	vp_scale (VP_ROTATED_HEIGHT / (xmax - xmin), 
		  (VP_ROTATED_HEIGHT / VP_SCREEN_RATIO) / (ymax - ymin));
    } else {
	vp_scale ((VP_STANDARD_HEIGHT / VP_SCREEN_RATIO) / (xmax - xmin), 
		  VP_STANDARD_HEIGHT / (ymax - ymin));
    }
}

void vp_style (vp_plotstyle st) {
    style = st;
    putchar (VP_SETSTYLE);
    switch (style) {
	case ROTATED:
	    putchar ('r');
	    break;
	case ABSOLUTE:
	    putchar ('a');
	    break;
	case STANDARD:
	default:
	    putchar ('s');
	    break;
    }
}

void vp_text (float x, float y, int size, int orient, const char *string)
{
    if (0 == size) return;

    pout (x, y, 0);
    putchar (VP_TEXT);
    vp_putint (size);
    vp_putint (orient);

    do {
	putchar (*string);
    }
    while (*string++);
}

void vp_tfont (int font1, int prec1, int ovly1)
{
    font = font1;
    prec = prec1;
    ovly = ovly1;

    putchar (VP_TXFONTPREC);
    vp_putint (font);
    vp_putint (prec);
    vp_putint (ovly);
}

void vp_tjust (int xjust1, int yjust1)
{
    xjust = xjust1;
    yjust = yjust1;

    putchar (VP_TXALIGN);
    vp_putint (xjust);
    vp_putint (yjust);
}

void vp_clip (float xmin, float ymin, float xmax, float ymax)
{
    putchar (VP_WINDOW);
    vp_putfloat (xmin);
    vp_putfloat (ymin);
    vp_putfloat (xmax);
    vp_putfloat (ymax);
}


void vp_uclip (float xmin, float ymin, float xmax, float ymax)
{
    xmin = fx + (xmin - ufx) * xscl;
    ymin = fy + (ymin - ufy) * yscl;
    xmax = fx + (xmax - ufx) * xscl;
    ymax = fy + (ymax - ufy) * yscl;
    vp_clip (xmin, ymin, xmax, ymax);
}

void vp_ugtext (float x, float y, 
		float xpath, float ypath,
		float xup, float yup, const char *string)
{
    x = fx + (x - ufx) * xscl;
    y = fy + (y - ufy) * yscl;
    xpath *= xscl;
    ypath *= yscl;
    xup *= xscl;
    yup *= yscl;

    vp_ugtext (x, y, xpath, ypath, xup, yup, string);
}

void vp_uplot (float x, float y, bool down)
{
    x = fx + (x - ufx) * xscl;
    y = fy + (y - ufy) * yscl;
    vp_plot (x, y, down);
}

void vp_uraster (unsigned char *array, int blast, int bit, int offset,
		 int xpix, int ypix, 
		 float xll, float yll, float ppi,
		 float *xur,float *yur, int orient, int invert)
{
    float x1, y1, x2, y2;
    
    x1 = fx + (xll - ufx) * xscl;
    y1 = fy + (yll - ufy) * yscl;

    if (0. == ppi) {
	x2 = fx + (*xur - ufx) * xscl;
	y2 = fy + (*yur - ufy) * yscl;
    }

    vp_raster (array, blast, bit, offset, 
	       xpix, ypix, 
	       x1, y1, ppi, &x2, &y2, orient, invert);

    if (ppi != 0. && xscl != 0. && yscl != 0.) {
	*xur = (x2 - fx)/xscl + ufx;
	*yur = (y2 - fy)/yscl + ufy;
    }
}

void vp_utext (float x, float y, int size, int orient, const char *string)
{
    x = fx + (x - ufx) * xscl;
    y = fy + (y - ufy) * yscl;
    vp_text (x, y, size, orient, string);
}

void vp_where (float *x, float *y)
{
    *x = xold;
    *y = yold;
}

void vp_color (int col)
{
    putchar (VP_COLOR);
    vp_putint (col);
}

/* plot an arrow from (x1,y1) to (x,y) with arrow-size r */
void vp_arrow (float x1, float y1, float x, float y, float r)
{
    float beta, alpha, xp[4], yp[4];
    const float pio4=0.785398;
    bool flag;

    flag = (r < 0.);
    if (flag) r = -r;

    if (x == x1 && y == y1) {
	xp[0] = x - r / 3.;
	yp[0] = y - r / 3.;
	xp[1] = x - r / 3.;
	yp[1] = y + r / 3.;
	xp[2] = x + r / 3.;
	yp[2] = y + r / 3.;
	xp[3] = x + r / 3.;
	yp[3] = y - r / 3.;
	if (flag) {
	    vp_area (xp, yp, 4, 0, 0, 0);
	} else {
	    vp_area (xp, yp, 4, 0, 1, 1);
	}
    } else {
	vp_move (x1, y1);
	vp_draw (x, y);

	beta = atan2f (y - y1, x - x1);

	xp[0] = x;
	yp[0] = y;
	alpha = pio4 + beta;
	xp[1] = x - r * cosf (alpha);
	yp[1] = y - r * sinf (alpha);
	alpha = pio4 - beta;
	xp[2] = x - r * cosf (alpha);
	yp[2] = y + r * sinf (alpha);
	if (flag) {
	    vp_area (xp, yp, 3, 0, 0, 0);
	} else {
	    vp_area (xp, yp, 3, 0, 1, 1);
	}
    }
}


/* plot an arrow from (x1,y1) to (x,y) with arrow-size r */
void vp_uarrow (float x1, float y1, float x, float y, float r)
{
    float beta, alpha, xp[4], yp[4];
    const float pio4=0.785398;
    bool flag;

    flag = (r < 0.);
    if (flag) r = -r;

    if (x == x1 && y == y1) {
	xp[0] = x - r / 3.;
	yp[0] = y - r / 3.;
	xp[1] = x - r / 3.;
	yp[1] = y + r / 3.;
	xp[2] = x + r / 3.;
	yp[2] = y + r / 3.;
	xp[3] = x + r / 3.;
	yp[3] = y - r / 3.;
	if (flag) {
	    vp_uarea (xp, yp, 4, 0, 0, 0);
	} else {
	    vp_uarea (xp, yp, 4, 0, 1, 1);
	}
    } else {
	vp_umove (x1, y1);
	vp_udraw (x, y);

	beta = atan2f (y - y1, x - x1);

	xp[0] = x;
	yp[0] = y;
	alpha = pio4 + beta;
	xp[1] = x - r * cosf (alpha);
	yp[1] = y - r * sinf (alpha);
	alpha = pio4 - beta;
	xp[2] = x - r * cosf (alpha);
	yp[2] = y + r * sinf (alpha);
	if (flag) {
	    vp_uarea (xp, yp, 3, 0, 0, 0);
	} else {
	    vp_uarea (xp, yp, 3, 0, 1, 1);
	}
    }
}

void vp_dash (float dash1, float gap1, float dash2, float gap2)
{
    if (dash1 < 0. || gap1 < 0. || dash2 < 0. || gap2 < 0.) {
	dashon = false;
	return;
    }

    ddef[0] = dash1;
    ddef[1] = dash1 + gap1;
    ddef[2] = ddef[1] + dash2;
    ddef[3] = ddef[2] + gap2;

    if (ddef[3] <= 0.) {
	dashon = false;
    } else {
	dashon = true;
	dashpos = 0.0;
    }
}
