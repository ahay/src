/* Plot signal with lollipops.

Takes: < data.rsf > plot.vpl
*/

#include <math.h>

#include <rsf.h>
#include <rsfplot.h>

static float screenwide, screenhigh;
static bool transp;
static const float eps=1.e-20;

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
    int i, i1,n1, i2,n2, i3, n3, ir, labelsz, connect, corners, dots, newsize;
    float **data, xxscale, yyscale, clip, f, vx[5], vy[5];
    float epsilon, dd1, dd2, axis, hi=0., lo=0., av, maxab, range;
    float marginl, marginr, margint, marginb, x, y, radius;
    float tracehigh, overlap, zerosignal, o1,d1;
    float d, full, abs, sgn, signal, *cx, *cy;
    char *label1, *title, **labels;
    bool seemean, strings, silk, gaineach, yreverse, constsep, seedead;
    sf_file in;

    sf_init(argc,argv);
    in = sf_input ("in");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if(!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if(!sf_histint(in,"n2",&n2)) n2=1;
    n3 = sf_leftsize(in,2);
 
    if(!sf_getfloat("o1",&o1) && !sf_histfloat(in,"o1",&o1)) o1=0.;
    if(!sf_getfloat("d1",&d1) && !sf_histfloat(in,"d1",&d1)) d1=1.;

    if (NULL == (label1 = sf_getstring("label1"))) 
	/* label for the axis */
	label1 = sf_histstring(in,"label1");
    if (label1 != NULL &&
	(*label1 == '\0' || (*label1 == ' ' && *(label1+1) == '\0'))) {
	free (label1);
	label1 = NULL;
    }
 
    if (NULL == (title = sf_getstring("title")))
	/* plot title */ 
	title = sf_histstring(in,"title");
    if (title != NULL &&
	(*title == '\0' || (*title == ' ' && *(title+1) == '\0'))) {
	free (title);
	title = NULL;
    }

    data = sf_floatalloc2 (n1,n2);
 
    if (!sf_getint("dots",&dots)) dots = (n1 <= 130)? 1: 0;
    /* type of dots: 1 - baloon, 0 - no dots, 2 - only for non-zero data */
    if (!sf_getbool("seemean",&seemean)) seemean = (n2 <= 30);
    /* if y, draw axis lines */
    if (!sf_getbool("strings",&strings)) strings = (n1 <= 400);
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

    if (!sf_getfloat("xxscale",&xxscale)) xxscale=1.;
    /* x scaling */
    if (!sf_getfloat("yyscale",&yyscale)) yyscale=1.;
    /* y scaling */
    if (!sf_getfloat("clip",&clip)) clip=-1.; 
    /* data clip */
    if (!sf_getfloat("overlap",&overlap)) overlap=0.9; 
    /* trace overlap */

    screenwide = VP_STANDARD_HEIGHT/VP_SCREEN_RATIO * xxscale;
    screenhigh = VP_STANDARD_HEIGHT  * yyscale;  
    epsilon = .0002 * screenhigh;

    marginl = screenwide * ((NULL == labels[0])? 0.03: 0.15);
    margint = screenhigh * ((NULL == title)? 0.03: 0.07);
    marginb = screenhigh * ((NULL == label1)? 0.03: 0.15);
    marginr = screenwide * 0.03;

    dd1 = (screenwide - marginl - marginr) / ( n1             );
    dd2 = (screenhigh - marginb - margint) / ((n2-1) + overlap);

    if(!sf_getfloat("radius",&radius)) radius = dd1/3;
    /* dot radius */
    if(radius > dd2/15) radius = dd2/15;

    tracehigh = overlap * (dots? dd2 - 3*radius: dd2);

    for (i3=0; i3 < n3; i3++) {
	sf_read(data[0],sizeof(float),n1*n2,in);
	
	vp_erase();
    
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
	    
	    av = (hi + lo) / 2.;
	    range = (hi>lo)? 1./(hi-lo): 1./eps;
	    zerosignal = axis +  tracehigh * (0.-av)*range;
	    
	    if(NULL != labels[0] && NULL != labels[i2]) {
		vp_color(5);
		text(0.03*screenwide, zerosignal+.2*dd2, 
		     labelsz, 0, labels[i2]);
	    }
	    
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
		    vp_color(1);
		    move(marginl, zerosignal);
		    draw(-marginr+screenwide, zerosignal);
		}

		if (connect) {
		    vp_color(7);
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
		}

		for (i1=0; i1<n1; i1++) {
		    y = axis +  tracehigh * (data[i2][i1]-av)*range;
		    x = marginl + dd1/2 + i1*dd1;
		    vp_color( 6);
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
	    }
	} /* i2 */

	if(NULL != label1) 
	    vp_simple_axis(marginl, marginb*0.8,  
			   screenwide-marginr, marginb*0.8,
			   o1, o1+(n1-1)*d1, 0., 0., 
			   .25, label1, 0.03*labelsz);

	if(NULL != title) {
	    newsize = 1.2*labelsz;
	    x = marginl + .5*( screenwide-marginl-marginr);
	    y = screenhigh - margint;

	    vp_color( 5);
	    vp_tjust(TH_CENTER,TV_NORMAL);
	    text (x, y, newsize, 0, title);
	    vp_tjust(TH_NORMAL,TV_NORMAL);
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

    vp_color(3);
    vp_area(vx,vy,corners,1,1,1);
}

/* 	$Id: dots.c,v 1.5 2003/10/01 23:41:18 fomels Exp $	 */

