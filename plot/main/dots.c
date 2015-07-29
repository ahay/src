/* Plot signal with lollipops.

Takes: > plot.vpl

The axis is displayed only if label1= is present in the input
file or the command line.  
*/
/*
  Copyright (C) 1987 The Board of Trustees of Stanford University
  Copyright (C) 2004 University of Texas at Austin

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
/* Modified from the original version by Jon F. Claerbout */

#include <math.h>

#include <rsf.h>
#include <rsfplot.h>

static float screenwide, screenhigh, screenratio;
static bool transp;

#define XOF(Y) (screenwide *       (Y) / screenhigh)
#define YOF(X) (screenhigh * (1. - (X) / screenwide))

static void move(float x, float y);
static void draw(float x,  float y);
static void text(float x, float y, 
		 int size, int orient, char *label);
static void circle(int corners, 
		   float *vx, float *vy, float x, float y, float radius);

int main (int argc, char* argv[])
{
    int i, i1,n1, i2,n2, i3, n3, ir;
    int labelsz, connect, corners, dots, newsize, font;
    size_t len;
    float **data, xxscale, yyscale, clip, f, vx[5], vy[5];
    float epsilon, dd1, dd2, axis, hi=0., lo=0., av, maxab, range;
    float marginl, marginr, margint, marginb, x, y, radius;
    float tracehigh, overlap, zerosignal, o1,d1;
    float d, full, abs, sgn, signal, *cx, *cy;
    char *label, *label1, *unit1, *title, **labels;
    bool seemean, strings, silk, gaineach, yreverse, constsep, seedead;
    sf_file in=NULL;

    sf_init(argc,argv);
    in = sf_input ("in");
    vp_init();

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if(!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if(!sf_histint(in,"n2",&n2)) n2=1;
    n3 = sf_leftsize(in,2);
 
    if(!sf_getfloat("o1",&o1) && !sf_histfloat(in,"o1",&o1)) o1=0.;
    if(!sf_getfloat("d1",&d1) && !sf_histfloat(in,"d1",&d1)) d1=1.;

    if (NULL == (label1 = sf_getstring("label1"))) 
	/*( label1 label for the axis )*/
	label1 = sf_histstring(in,"label1");
    if (label1 != NULL &&
	(*label1 == '\0' || (*label1 == ' ' && *(label1+1) == '\0'))) {
	free (label1);
	label1 = NULL;
    }

    if (NULL == (unit1 = sf_getstring("unit1"))) 
	/*( unit1 unit for the axis )*/
	unit1 = sf_histstring(in,"unit1");
    if (unit1 != NULL &&
	(*unit1 == '\0' || (*unit1 == ' ' && *(unit1+1) == '\0'))) {
	free (unit1);
	unit1 = NULL;
    }

    if (NULL != label1 && NULL != unit1) {
	len = strlen(label1)+strlen(unit1)+4;
	label = sf_charalloc(len);
	snprintf(label,len,"%s (%s)",label1,unit1);
	free(label1);
	free(unit1);
    } else {
	label = label1;
    }
 
    if (NULL == (title = sf_getstring("title")))
	/*( title plot title )*/ 
	title = sf_histstring(in,"title");
    if (title != NULL &&
	(*title == '\0' || (*title == ' ' && *(title+1) == '\0'))) {
	free (title);
	title = NULL;
    }

    data = sf_floatalloc2 (n1,n2);
 
    if (!sf_getint("dots",&dots)) dots = (n1 <= 130)? 1: 0;
    /* type of dots: 1 - baloon, 0 - no dots, 2 - only for non-zero data */
    if (!sf_getbool("seemean",&seemean)) seemean = (bool) (n2 <= 30);
    /* if y, draw axis lines */
    if (!sf_getbool("strings",&strings)) strings = (bool) (n1 <= 400);
    /* if y, draw strings */
    if (!sf_getint("connect",&connect)) connect = 1; 
    /* connection type: 1 - diagonal, 2 - bar, 4 - only for non-zero data */
    if (!sf_getint("corners",&corners)) {
	/* number of polygon corners (default is 6) */
	corners=7;
    } else {
	corners++;
    }

    cx = sf_floatalloc(corners);
    cy = sf_floatalloc(corners);

    if (!sf_getbool("silk",&silk)) silk=false; 
    /* if y, silky plot */
    if (!sf_getbool("gaineach",&gaineach)) gaineach=true;
    /* if y, gain each trace independently */
    if (!sf_getint("labelsz",&labelsz)) labelsz=8;
    /* label size */
    if (!sf_getbool("yreverse",&yreverse)) yreverse=false;
    /* if y, reverse y axis */
    if (!sf_getbool("constsep",&constsep)) constsep=false;
    /* if y, use constant trace separation */
    if (!sf_getbool("seedead",&seedead)) seedead=false;
    /* if y, show zero traces */
    if (!sf_getbool("transp",&transp)) transp=false;
    /*  if y, transpose the axis */

    labels = (char**) sf_alloc(n2,sizeof(char*));
    if (!sf_getstrings("labels",labels,n2)) labels[0] = NULL;
    /* trace labels */

    if (!sf_getfloat("xxscale",&xxscale)) xxscale=1.;
    /* x scaling */
    if (!sf_getfloat("yyscale",&yyscale)) yyscale=1.;
    /* y scaling */
    if (!sf_getfloat("clip",&clip)) clip=-1.; 
    /* data clip */
    if (!sf_getfloat("overlap",&overlap)) overlap=0.9; 
    /* trace overlap */

    /* get screen size */
    if (!sf_getfloat ("screenratio",&screenratio))
	screenratio = VP_SCREEN_RATIO;
    /* screen aspect ratio */
    if (!sf_getfloat ("screenht",&screenhigh))
	screenhigh = VP_STANDARD_HEIGHT;
    /* screen height */
    if (!sf_getfloat ("screenwd",&screenwide))
	screenwide = screenhigh / screenratio;
    /* screen width */

    screenwide *= xxscale;
    screenhigh *= yyscale;  
    epsilon = .0002 * screenhigh;

    marginl = screenwide * ((NULL == labels[0])? 0.03: 0.15);
    margint = screenhigh * ((NULL == title)? 0.03: 0.07);
    marginb = screenhigh * ((NULL == label)? 0.03: 0.15);
    marginr = screenwide * 0.03;

    if (NULL == title && NULL != labels[0]) margint += 0.03*labelsz;

    dd1 = (screenwide - marginl - marginr) / ( n1             );
    dd2 = (screenhigh - marginb - margint) / ((n2-1) + overlap);

    if(!sf_getfloat("radius",&radius)) radius = dd1/3;
    /* dot radius */
    if(radius > dd2/15) radius = dd2/15;

    tracehigh = overlap * (dots? dd2 - 3*radius: dd2); 

    if (!sf_getint("font",&font)) font=-1; /* font to use in text */

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(data[0],n1*n2,in);
	
	vp_erase();
	if (-1 != font) vp_tfont (font,VP_NO_CHANGE,VP_NO_CHANGE);

	if(!gaineach) {
	    hi = lo = data[0][0];
	    for(i2=0; i2<n2; i2++) {
		for(i1=0; i1<n1; i1++) {
		    f = data[i2][i1];
		    if(f > hi) {
			hi = f;
		    } else if (f < lo) {
			lo = f;
		    }
		}
	    }
	}
	
	for(i2=0; i2<n2; i2++) {
	    ir = yreverse? n2-1-i2: i2;
	    axis = marginb + dd2*(ir + 0.5*overlap);

	    if (clip > 0.) {
		hi =  clip;
		lo = -clip;
	    } else {	    
		if(gaineach) {
		    hi = lo = data[i2][0];
		    for(i1=0; i1<n1; i1++) {
			f = data[i2][i1];
			if(f > hi) {
			    hi = f;
			} else if (f < lo) {
			    lo = f;
			}
		    }
		}
		
		if(constsep || silk) {
		    maxab = ( -lo > hi)? -lo: hi;
		    hi =  maxab;
		    lo = -maxab;
		}
	    }

	    if (hi > 0.0 && lo > 0.0) {
		lo=-lo;
	    } else if (hi < 0.0 && lo < 0.0) {
		hi=-hi;
	    } else if (hi==0.0) {
		hi=1.0;
		lo=-1.0;
	    }

	    av = (hi + lo) / 2.;
	    range = 1./(hi-lo);

	    zerosignal = axis +  tracehigh * (0.-av)*range;

	    vp_bgroup("labels");
	    if(NULL != labels[0] && NULL != labels[i2]) {
		vp_color(VP_CYAN);
		text(0.03*screenwide, zerosignal+.2*dd2, 
		     labelsz, 0, labels[i2]);
	    }
	    vp_egroup();

	    if (silk) {
		for (i1=0; i1<n1; i1++) {
		    d = 2. * data[i2][i1] * range;
		    abs = (d > 0.) ? d  : -d ;
		    sgn = (d > 0.) ? 1. : -1.;
		    full = .4 + .45 * sgn * sqrt(abs);
		    full = (full > 0) ? full : 0.;
		    y = axis;
		    x = marginl + dd1/2 + i1*dd1;
		    i = 0;
		    vx[i] = x + dd1*full; vy[i] = y           ; i++;
		    vx[i] = x           ; vy[i] = y + dd2*full; i++;
		    vx[i] = x - dd1*full; vy[i] = y           ; i++;
		    vx[i] = x           ; vy[i] = y - dd2*full; i++;
		    vx[i] = x + dd1*full; vy[i] = y           ; i++;
		    vp_area(vx,vy,5,1,1,1);
		}
	    } else {
		if(!seedead) {
		    for(i1=0; i1<n1; i1++) {
			if (data[i2][i1] != 0.) break;
		    }
		    if (i1==n1) continue;
		}

		if(seemean) {
		    vp_bgroup("axis");
		    vp_color(VP_BLUE);
		    move(marginl, zerosignal);
		    draw(-marginr+screenwide, zerosignal);
		    vp_egroup();
		}

		if (connect) {
		    vp_color(VP_WHITE);
		    vp_bgroup("connectors");
		    for(i1=0; i1<n1; i1++) {
			y = axis +  tracehigh * (data[i2][i1]-av)*range;
			switch (connect) {
			    case 4: /* no segment */
				x = marginl + dd1/2 + i1*dd1;
				if (i1 == 0) {
				    move (x, y);
				} else if(data[i2][i1]   != 0. &&  
					  data[i2][i1-1] != 0.) {
				    draw (x, y); 
				} else { 
				    move (x, y);
				}
				break;
			    case 1: /* diagonal segments */
				x = marginl + dd1/2 + i1*dd1;
				if (i1 == 0) {
				    move (x, y);
				} else {
				    draw (x, y);
				}
				break;
			    case 2: /* bar graph */
				x = marginl + i1*dd1;
				if (i1 == 0) {
				    move (x, y);
				} else {
				    draw (x, y);
				}
				draw (x + dd1, y); /* bar */
				break;
			    case 3: /* horizontal segments */
				x = marginl + i1*dd1;
				move (x, y);
				draw (x + dd1, y);
				break;
			}
		    }
		    vp_egroup();
		}

		vp_bgroup("data");
		for (i1=0; i1<n1; i1++) {
		    y = axis +  tracehigh * (data[i2][i1]-av)*range;
		    x = marginl + dd1/2 + i1*dd1;
		    vp_color(VP_YELLOW);
		    if (strings) {
			signal = y - zerosignal;
			if(fabsf(signal)  > epsilon ) {
			    move (x, y);
			    draw (x, zerosignal);
			}
		    }

		    if (1==dots || (2==dots && data[i2][i1] != 0.))
			circle (corners, cx, cy, x, y, radius);
		}	
		vp_egroup();
	    }
	} /* i2 */

	if(NULL != label) 
	    vp_simple_axis(marginl+dd1/2,        marginb*0.8,  
			   marginl+(n1-0.5)*dd1, marginb*0.8,
			   o1, o1+(n1-1)*d1, 0., 0., 
			   .25, label, 0.03*labelsz);

	if(NULL != title) {
	    newsize = 1.2*labelsz;
	    x = marginl + .5*( screenwide-marginl-marginr);
	    y = screenhigh - margint;

	    vp_bgroup("title");
	    vp_color(VP_CYAN);
	    vp_tjust(TH_CENTER,TV_NORMAL);
	    text (x, y, newsize, 0, title);
	    vp_tjust(TH_NORMAL,TV_NORMAL);
	    vp_egroup();
	}

    } /* i3 */


    exit(0);
}

static void move(float x, float y)
{
    if (transp) {
	vp_move( XOF(y), YOF(x));
    } else {
	vp_move(     x,      y);
    }
}

static void draw(float x,  float y)
{
    if (transp) { 
	vp_draw( XOF(y), YOF(x));
    } else {
	vp_draw(     x,      y);
    }
}

static void text(float x, float y, 
		 int size, int orient, char* label)
{
    if (transp) {
	vp_text( XOF(y), YOF(x), size, orient, label);
    } else {
	vp_text(     x,      y , size, orient, label);
    }
}

static void circle(int corners, 
		   float *vx, float *vy, float x, float y, float radius)
{
    float arg;
    int i;

    for(i=0; i<corners; i++) {
	arg = (2.*SF_PI*i)/(corners-1);
	vx[i] = x + radius * sinf(arg);
	vy[i] = y + radius * cosf(arg);
    }

    vp_color(VP_PURPLE);
    vp_area(vx,vy,corners,1,1,1);
}
