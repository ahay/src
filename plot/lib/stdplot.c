/* Setting up frames for a generic plot. */
/*
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

#include <math.h>
#include <stdio.h>
#include <string.h>

#include <rsf.h>
/*^*/

#include "stdplot.h"
#include "axis.h"
#include "vplot.h"
#include "plot.h"

static float min1,min2, max1,max2, mid1,mid2, inch1,inch2, orig1,orig2, inch3, tickscale, tickscale1, tickscale2, tickscale3, tickscale4;
static float labelsz, barlabelsz, barmin, barmax, bar0, dbar, sinth, costh, griddash;
static float d1, d2, d3, frame1, frame2, frame3, l1min, l1max, l2min, l2max, l3min, l3max;
static int framecol, gridcol, gridfat=1;
static int cubelinecol=VP_WHITE;
static int framelabelcol=VP_YELLOW;
static bool labelrot, transp, barmove, wheretics, scalebar, vertbar, wherebartics;
static bool cube=false, movie=false, flat, xreverse, yreverse;
static char blank[]=" ";
static const float aspect=0.8, ticksep=1.75;

static struct Label {
    char where;
    int fat;
    float x, y, xpath, ypath, xup, yup, min, max;
    char* text;
}   /*@null@*/ *label1=NULL, 
    /*@null@*/ *label2=NULL, 
    /*@null@*/ *label3=NULL, 
    /*@null@*/ *label4=NULL,
    /*@null@*/ *title=NULL, 
    /*@null@*/ *barlabel=NULL;

static struct Axis {
    bool parallel;
    int ntic, maxstrlen;
    float orig, dnum, num0;
    const char *format;
}   /*@null@*/ *axis1=NULL, 
    /*@null@*/ *axis2=NULL, 
    /*@null@*/ *axis3=NULL, 
    /*@null@*/ *axis4=NULL, 
    /*@null@*/ *baraxis=NULL, 
    /*@null@*/ *grid1=NULL, 
    /*@null@*/ *grid2=NULL, 
    /*@null@*/ *grid3=NULL;

static void make_title (sf_file in, char wheret);
static void make_barlabel (sf_file in);
static void make_labels (sf_file in, char where1, char where2);
static void make_baraxis (float min, float max);
static void make_axes (void);
static void make_grid (bool grid);
static void swap(float*a,float*b);

static void swap(float*a,float*b)
{
    float f;

    f = *a; *a = *b; *b=f;
}

void vp_stdplot_init (float umin1, float umax1 /* user's frame for axis 1 */, 
		      float umin2, float umax2 /* user's frame for axis 2 */,
		      bool transp1             /* default transpose flag */, 
		      bool xreverse1           /* default x reverse flag */, 
		      bool yreverse1           /* default y reverse flag */, 
		      bool pad1                /* default padding */)
/*< Initializing standard plot >*/
{
    bool pad, set;
    int font;
    float mid, off, crowd, barwd;
    float xll, xur, yll, yur, screenratio, screenht, screenwd, marg;
    char* bartype;

    transp = transp1;
    if (!sf_getbool ("xreverse",&xreverse)) xreverse = xreverse1;
    /* reverse horizontal axis */
    if (!sf_getbool ("yreverse",&yreverse)) yreverse = yreverse1;
    /* reverse vertical axis */
    if (!sf_getbool ("pad",&pad)) pad = pad1;
    /* pad plotting area */
    if (!sf_getbool ("wantscalebar",&scalebar) && 
	!sf_getbool ("scalebar",&scalebar)) scalebar = false;
    /* plot a scalebar */

    if (!sf_getbool ("barmove", &barmove)) barmove = false;
    /* adjust scalebar position, if bartype=h */

    /* plot a scalebar */
    if (pad) { /* 4% stretch */
	mid = 0.5*(umin1+umax1);
	off = 1.04*0.5*(umax1-umin1);
	umin1 = mid-off;
	umax1 = mid+off;
	
	mid = 0.5*(umin2+umax2);
	off = 1.04*0.5*(umax2-umin2);
	umin2 = mid-off;
	umax2 = mid+off;
    }

    if (!sf_getfloat ("tickscale",&tickscale)) tickscale=0.5;
    /* ticks scaling */
    if (!sf_getfloat ("tickscale1",&tickscale1)) tickscale1=tickscale;
    /* ticks scaling on first axis */
    if (!sf_getfloat ("tickscale2",&tickscale2)) tickscale2=tickscale;
    /* ticks scaling on second axis */    
    if (!sf_getfloat ("tickscale3",&tickscale3)) tickscale3=tickscale;
    /* ticks scaling on third axis */
    if (!sf_getfloat ("tickscale4",&tickscale4)) tickscale4=tickscale;
    /* ticks scaling on fourth axis */    

    /* get max and min */
    if (!sf_getfloat ("min1",&min1)) min1=umin1;
    /* minimum on the first axis */
    if (!sf_getfloat ("min2",&min2)) min2=umin2;
    /* minimum on the second axis */
    if (!sf_getfloat ("max1",&max1)) max1=umax1;
    /* maximum on the first axis */
    if (!sf_getfloat ("max2",&max2)) max2=umax2;
    /* maximum on the second axis */

    if (transp) {
	swap(&min1,&min2);
	swap(&max1,&max2);
    }

    if (xreverse) swap(&min1,&max1);
    if (yreverse) swap(&min2,&max2);
    
    vp_style (VP_STANDARD);

    if (!sf_getint("font",&font)) font=-1; /* font to use in text */
    if (-1 != font) vp_tfont (font,VP_NO_CHANGE,VP_NO_CHANGE);

    /* get screen size */
    if (!sf_getfloat ("screenratio",&screenratio))
	screenratio = VP_SCREEN_RATIO;
    /* ratio of screen height to screen width */
    if (!sf_getfloat ("screenht",&screenht))
	screenht = VP_STANDARD_HEIGHT;
    /* screen height (in "inches") */
    if (!sf_getfloat ("screenwd",&screenwd))
	screenwd = screenht / screenratio;
    /* screen width (in "inches") */

    /* get inches */
    if (!sf_getfloat ("crowd",&crowd)) crowd  = 0.75;
    /* window size as ratio of the total area */
    if (!sf_getfloat ("xinch",&inch1)) {
	/* window width as ratio of the screen width */
	if (!sf_getfloat ("crowd1",&inch1)) inch1 = crowd;
	/* window width as ratio of the screen width */
	inch1 *= screenwd;
    }
    if (!sf_getfloat ("yinch",&inch2)) {
	/* window height as ratio of the screen height */
	if (!sf_getfloat ("crowd2",&inch2)) inch2 = crowd;
	/* window height as ratio of the screen height */
	inch2 *= screenht;
    }

    /* get corners */
    set = sf_getfloat ("xll",&xll);
    /* lower left horizontal coordinate */
    if (!sf_getfloat ("xur",&xur)) {
	/* upper right horizontal coordinate */
	if (set) {
	    xur = xll + inch1;
	} else {
	    marg =  screenwd - inch1;
	    xll = marg * 2./3.;
	    xur = screenwd - marg * 1./3.;
	}
    } else if (!set) {
	xll = xur - inch1;
    }
	
    set = sf_getfloat ("yll",&yll);
    /* lower left vertical coordinate */
    if (!sf_getfloat ("yur",&yur)) {
	/* upper right vertical coordinate */
	if (set) {
	    yur = yll + inch2;
	} else {
	    marg =  screenht - inch2;
	    yll = marg * 1./2.;
	    yur = screenht - marg * 1./2.;
	}
    } else if (!set) {
	yll = yur - inch2;
    }

    /* make frame smaller to accomodate scale bar */
    if (scalebar) {
	bartype = sf_getstring("bartype");
	/* [v,h] vertical or horizontal bar (default is v) */
	if (NULL == bartype) {
	    vertbar = true;
	} else {
	    vertbar = (bool) (bartype[0] == 'v');
	    free (bartype);
	}
	if (!sf_getfloat("barwidth",&barwd)) barwd = 0.36;
	/* scale bar size */
	if (vertbar) {
	    xur -= (0.07*screenwd + barwd);
	    xll -= 0.5*barwd;
	    barmax = 0.88*screenwd;
	    barmin = barmax-barwd;
	} else {
	    if (cube) {
                if (!barmove) {
		    yll -= 0.4*barwd;
   		    yur -= (0.04*screenht + barwd);
	   	    barmax = 0.98*screenht;
		    barmin = barmax-barwd;
                } else {
                    barmin = yll - 0.05*screenht;
                    barmax = barmin+barwd;
                    yur += 0.4*barwd;
                    yll += (0.05*screenht + barwd);
                }
	    } else {
		barmin = barmove ? yll : 0.12*screenht;
		barmax = barmin+barwd;
		yur += 0.4*barwd;
		yll += (0.04*screenht + barwd);
	    }
	}
    }

    /* set origin */
    if (xll != xur) {
	inch1 = xur - xll;
	orig1 = xll + (xur - xll)*0.5;
    } else {
	orig1 = screenwd*0.5;
    }

    if (yll != yur) {
	inch2 = yur - yll;
	orig2 = yll + (yur - yll)*0.5;
    } else {
	orig2 = screenht*0.5;
    }
    
    if (min1 == 0. && max1 == 0.) {
	min1 = -1.;
	max1 = 1.;
    } else if (min1 == max1) {
        max1 *= 1.04;
        min1 *= 0.96;
    }
    
    if (min2 == 0. && max2 == 0.) {
	min2 = -1.;
	max2 = 1.;
    } else if (min2 == max2) {
        max2 *= 1.04;
        min2 *= 0.96;
    }

    if (!sf_getint ("axiscol",&framecol)) framecol=VP_WHITE;
    /* axes color */

    if (!sf_getint ("framelabelcol",&framelabelcol)) framelabelcol=VP_YELLOW;
    /* frame labels color */

    if (!sf_getint ("cubelinecol",&cubelinecol)) cubelinecol=framelabelcol;
    /* cube lines color */

    vp_coordinates();
}

void vp_cubecoord (int side /* 1=top, 2=side, 3=front */,
		   float umin1, float umax1, float umin2, float umax2)
/*< setup the coordinate system for the cube front >*/
{
    float scale1, scale2;

    /* find scales and user origin */
    switch (side) {
	case 3: /* front */
	    scale1 = inch1 * (mid1 - min1) / (max1 - min1) / (umax1 - umin1);
	    scale2 = inch2 * (mid2 - min2) / (max2 - min2) / (umax2 - umin2);
	    vp_orig (orig1-0.5*inch1,orig2-0.5*inch2);
	    break;
	case 2: /* side */
	    scale1 = inch1 * (max1 - mid1) / (max1 - min1) / (umax1 - umin1);
	    scale2 = inch2 * (mid2 - min2) / (max2 - min2) / (umax2 - umin2);
	    vp_orig (orig1+inch1*((mid1 - min1)/(max1 - min1)-0.5),
		     orig2-0.5*inch2);
	    break;
	case 1:
	default:
	    scale1 = inch1 * (mid1 - min1) / (max1 - min1) / (umax1 - umin1);
	    scale2 = inch2 * (max2 - mid2) / (max2 - min2) / (umax2 - umin2);
	    vp_orig (orig1-0.5*inch1,
		     orig2+inch2*((mid2 - min2)/(max2 - min2)-0.5));
    }
    
    /* setup the coordinate system */
    vp_scale (scale1, scale2);
    vp_uorig (umin1, umin2);
    vp_uclip (umin1, umin2, umax1, umax2);
}


void vp_coordinates (void)
/*< setup the coordinate system >*/
{
    float scale1, scale2, uorig1, uorig2;

    /* find scales and user origin */
    scale1 = inch1 / (max1 - min1);
    scale2 = inch2 / (max2 - min2);
    uorig2 = (min2 + max2)*0.5;
    uorig1 = (min1 + max1)*0.5;
    
    /* setup the coordinate system */
    vp_scale (scale1, scale2);
    vp_orig (orig1, orig2);
    vp_uorig (uorig1, uorig2);
}


void vp_cubeplot_init (int n1pix, int n2pix,      /* total pixels */ 
		       int n1front, int n2front,  /* front face pixels */
		       bool flat1                 /* flat flag */,
                       bool movie1                /* movie flag */) 
/*< Initializing 3-D cube plot. >*/
{
    min1 = -0.5;
    max1 = n1pix-0.5;
    mid1 = n1front-0.5;

    min2 = -0.5;
    max2 = n2pix-0.5;
    mid2 = n2front-0.5;

    cube = true;
    flat = flat1;
    movie = movie1;

    vp_stdplot_init (min1,max1,min2,max2, true, false, false, false);
    swap(&mid1,&mid2); /* transp=true */
}

static void make_labels (sf_file in, char where1, char where2)
{
    int n1, n2, n3;
    float vs, xc, yc, x1, y1, x2, y2, o1, o2, o3;
    bool want;
    char *where, *labl, *unit;
    size_t len;
    struct Label *label;

    if (cube) {
	if (!sf_histint(in,"n1",&n1)) n1=1;
	if (!sf_histfloat(in,"o1",&o1)) o1=0.;
	if (!sf_histfloat(in,"d1",&d1)) d1=1.;
	l2max = o1-0.5*d1;
	l2min = o1+(n1-0.5)*d1;

	if (!sf_histint(in,"n2",&n2)) n2=1;
	if (!sf_histfloat(in,"o2",&o2)) o2=0.;
	if (!sf_histfloat(in,"d2",&d2)) d2=1.;
	l1min = o2-0.5*d2;
	l1max = o2+(n2-0.5)*d2;

	if (!movie || !sf_histint(in,"n3",&n3)) n3 = sf_leftsize(in,2);
	if (!sf_histfloat(in,"o3",&o3)) o3=0.;
	if (!sf_histfloat(in,"d3",&d3)) d3=1.;
	l3min = o3-0.5*d3;
	l3max = o3+(n3-0.5)*d3;
    }

    if (sf_getbool ("wantaxis", &want) && !want) {
	/*( wantaxis if draw axes )*/
	label1 = NULL;
	label2 = NULL;
	if (cube) label3 = NULL;
	return;
    }

    if (sf_getbool ("wantaxis1",&want) && !want) {
	/*( wantaxis1 if draw first axis )*/
	label1 = NULL;
    } else if (NULL == label1) { 
	label1 = (struct Label*) sf_alloc(1,sizeof(struct Label));
    }

    if (sf_getbool ("wantaxis2",&want) && !want) {
	/*( wantaxis2 if draw second axis )*/
	label2 = NULL;
	if (!cube && NULL == label1) return;
    } else if (NULL == label2) {
	label2 = (struct Label*) sf_alloc(1,sizeof(struct Label));
    }

    if (cube) {
	if (sf_getbool ("wantaxis3",&want) && !want) {
	    /*( wantaxis3 if draw third axis in cube plots )*/
	    label3 = NULL;
	    if (NULL == label1 && NULL == label2) return;
	} else if (NULL == label3) {
	    label3 = (struct Label*) sf_alloc(1,sizeof(struct Label));
	}
    }

    if (transp) {
	label = label1;
	label1 = label2;
	label2 = label;
    }

    if (cube && flat && NULL != label3 && NULL != label2) {
	label4 = (struct Label*) sf_alloc(1,sizeof(struct Label));
    } else {
	label4 = NULL;
    }

    if (!sf_getfloat ("labelsz",&labelsz)) labelsz=8.;
    /* label size */
    if (cube) {
	labelsz *= 0.03; /* slightly smaller */
    } else {
	labelsz /= 33.;
    }
    vs = 2.5*labelsz;

    if (!sf_getbool ("labelrot",&labelrot)) labelrot = false;
    /* if rotate vertical label */
    
    if (NULL != label1) { /* horizontal axis */
	if (!cube && NULL != (where = sf_getstring("wherexlabel"))) {
	    /* where to put horizontal axis (top,bottom) */
	    label1->where = *where;
	} else {
	    label1->where = where1;
	}
	if (!sf_getint ("labelfat",&(label1->fat))) label1->fat=0;
	/*( labelfat label fatness )*/
	if ((NULL == (labl=sf_getstring(transp? "label2":"label1"))) &&
	    (NULL == (labl=sf_histstring(in,
					 transp? "label2":"label1")))) {
	    label1->text = blank;
	    /*( label1 label on the first axis )*/
	} else if (((NULL == (unit=sf_getstring(transp? "unit2":"unit1"))) &&
		    (NULL == (unit=sf_histstring(in,
						 transp? "unit2":"unit1")))) ||
		   *unit == '\0' || (*unit == ' ' && *(unit+1) == '\0')) {
	    label1->text = labl;
	    /*( unit1 unit on the first axis )*/
	} else {
	    len = strlen(labl)+strlen(unit)+4;
	    label1->text = sf_charalloc(len);
	    snprintf(label1->text,len,"%s (%s)",labl,unit);
	    free(labl);
	    free(unit);
	}
	
	label1->xpath = labelsz;
	label1->xup = 0.;
	label1->ypath = 0.;
	label1->yup = labelsz;

	xc = cube? 0.5*(mid1 + min1) : 0.5*(max1 + min1);
	yc = (label1->where == 't') ? max2: min2;

	vp_umove (xc, yc);
	vp_where (&xc, &yc);

	label1->x = xc;
	label1->y = (label1->where == 't') ? yc+vs: yc-vs;

	if (cube) {
	    label1->min = l1min;
	    label1->max = l1max;
	} else {
	    label1->min = min1;
	    label1->max = max1;
	}
    }

    if (NULL != label3) { /* out of plane axis */
	if (!sf_getint ("labelfat",&(label3->fat))) label3->fat=0;
	/* label fatness */
	if ((NULL == (labl=sf_getstring("label3"))) &&
	    (NULL == (labl=sf_histstring(in,"label3")))) {
	    label3->text = blank;
	    /*( label3 label on the third axis )*/
	} else if (((NULL == (unit=sf_getstring("unit3"))) &&
		    (NULL == (unit=sf_histstring(in,"unit3")))) ||
		   *unit == '\0' || (*unit == ' ' && *(unit+1) == '\0')) {
	    label3->text = labl;
	    /*( unit3 unit on the third axis )*/
	} else {
	    len = strlen(labl)+strlen(unit)+4;
	    label3->text = sf_charalloc(len);
	    snprintf(label3->text,len,"%s (%s)",labl,unit);
	    free(labl);
	    free(unit);
	}

	if (flat) {
	    vp_umove (mid1,min2);
	    vp_where (&x1, &y1);
	    vp_umove (max1,min2);
	    vp_where (&x2, &y1);
	    inch3 = fabsf(x2-x1);

	    label3->xpath = labelsz;
	    label3->xup = 0.;
	    label3->ypath = 0.;
	    label3->yup = labelsz;

	    xc = 0.5*(x2+x1);
	    yc = y1;

	    label3->y = yc-vs;	
	    label3->x = xc;
	} else {
	    vp_umove (mid1,min2);
	    vp_where (&x1, &y1);
	    vp_umove (max1,max2-mid2);
	    vp_where (&x2, &y2);
	    inch3 = hypotf(x2-x1,y2-y1);
	    costh = (x2-x1)/inch3;
	    sinth = (y2-y1)/inch3;

	    label3->xpath = labelsz*costh;
	    label3->xup = -labelsz*sinth;
	    label3->ypath = labelsz*sinth;
	    label3->yup = labelsz*costh;

	    xc = 0.5*(x1+x2);
	    yc = 0.5*(y1+y2);

	    label3->y = yc-vs*costh;	
	    label3->x = xc+vs*sinth;
	}

	label3->min = l3min;
	label3->max = l3max;
    }

    if (NULL != label2) { /* vertical axis */
	if (!cube && NULL != (where = sf_getstring("whereylabel"))) {
	    /* where to put vertical label (left,right) */
	    label2->where = *where;
	} else {
	    label2->where = where2;
	}
	if (!sf_getint ("labelfat",&(label2->fat))) label2->fat=0;
	/* label fatness */
	if ((NULL == (labl=sf_getstring(transp? "label1":"label2"))) &&
	    (NULL == (labl=sf_histstring(in,
					 transp? "label1":"label2")))) {
	    label2->text = blank;
	    /*( label2 label on the second axis )*/
	} else if (((NULL == (unit=sf_getstring(transp? "unit1":"unit2"))) &&
		    (NULL == (unit=sf_histstring(in,
						 transp? "unit1":"unit2")))) ||
		   *unit == '\0' || (*unit == ' ' && *(unit+1) == '\0')) { 
	    label2->text = labl;
	    /*( unit2 unit on the second axis )*/	    
	} else {
	    len = strlen(labl)+strlen(unit)+4;
	    label2->text = sf_charalloc(len);
	    snprintf(label2->text,len,"%s (%s)",labl,unit);
	    free(labl);
	    free(unit);
	}

	label2->ypath = labelrot? -labelsz: labelsz;
	label2->yup = 0.;
	label2->xpath = 0.;
	label2->xup = labelrot? labelsz: -labelsz;

	xc  = (label2->where == 'l')? min1: max1;
	yc = cube? 0.5*(min2 + mid2) : 0.5*(min2 + max2);

	vp_umove (xc, yc);
	vp_where (&xc, &yc);

	label2->y = yc;	
	label2->x = (label2->where == 'l')? xc-vs: xc+vs;
	
	if (cube) {
	    label2->max = l2max;
	    label2->min = l2min;
	} else {
	    label2->min = min2;
	    label2->max = max2;
	}
    } 

    if (NULL != label4) { /* vertical axis for the out of plane axis */
	label4->where = label2->where;
	label4->fat   = label2->fat;
	label4->text  = label3->text;

	label4->ypath = label2->ypath;
	label4->yup   = label2->yup;
	label4->xpath = label2->xpath;
	label4->xup   = label2->xup;

	xc  = min1;
	yc  = 0.5*(max2 + mid2);

	vp_umove (xc, yc);
	vp_where (&xc, &yc);

	label4->y = yc;	
	label4->x = label2->x;
      
	label4->max = l3max;
	label4->min = l3min;
    } 
}

static void make_baraxis (float min, float max)
{
    char* where;
    bool modify;

    if (barlabel == NULL) return;

    if (NULL == baraxis) 
	baraxis = (struct Axis*) sf_alloc(1,sizeof(struct Axis));

    if (!sf_getfloat ("axisorbar",&(baraxis->orig))) {
	/* bar origin */
	if (vertbar) {
	    baraxis->orig = (barlabel->where == 'l'? barmin: barmax);
	} else {
	    baraxis->orig = (barlabel->where == 'b'? barmin: barmax);
	}
    }

    if (!sf_getbool ("parallelbar",&(baraxis->parallel)) &&
	(( vertbar && !sf_getbool ("parallel2",&(baraxis->parallel))) ||
	 (!vertbar && !sf_getbool ("parallel1",&(baraxis->parallel)))) &&
	!sf_getbool ("parallel", &(baraxis->parallel))) baraxis->parallel=true;
	/* tickmark labels parallel to the axis */

    if (!sf_getint ("nbartic",&(baraxis->ntic))) 
	baraxis->ntic = 0.5 + (vertbar? inch2: inch1)/(aspect*labelsz);
    /* nbartic number of scalebar ticmarks */

    if (NULL == (baraxis->format = sf_getstring("formatbar"))) {
	/* format for ticmark labels in the scalebar */
	if ((vertbar  && NULL == (baraxis->format = sf_getstring("format2"))) ||
	    (!vertbar && NULL == (baraxis->format = sf_getstring("format1"))))
	    baraxis->format="%1.5g";
    }

    modify = (bool) (!sf_getfloat ("dbarnum", &(baraxis->dnum)) ||
		     /*( dbarnum scalebar tic increment )*/
		     !sf_getfloat ("obarnum", &(baraxis->num0)));
    /*( obarnum scalebar tic origin )*/

    baraxis->ntic = vp_optimal_scale(baraxis->ntic, 
				     modify,
				     baraxis->parallel,
				     baraxis->format,
				     min, max, 
				     &(baraxis->num0), 
				     &(baraxis->dnum),
				     &(baraxis->maxstrlen));

    if (fabsf(max-min) < SF_EPS) {
        baraxis->ntic = 1;
        bar0 = 0.0;
        dbar = SF_EPS;
    } else {
        bar0 = (baraxis->num0-0.5*(min+max))/(max-min);
        dbar = (baraxis->dnum)/(max-min);
    }

    if (vertbar) {
	bar0 = orig2 + inch2*bar0;
	dbar *= inch2;
    } else {
	bar0 = orig1 + inch1*bar0;
	dbar *= inch1;
    }

    wherebartics = (bool) ((NULL != (where = sf_getstring ("wherebartics"))) &&
			   ('a' == *where));
    /*( wherebartics where to put scalebar ticmarks )*/
}	

static void make_axes (void)
/* put tickmarks */
{
    char *where;
    bool modify;

    if (label1 != NULL) { /* horizontal axis */ 
	if (NULL == axis1) 
	    axis1 = (struct Axis*) sf_alloc(1,sizeof(struct Axis));

	if (!sf_getfloat ("axisor1",&(axis1->orig)))
	    axis1->orig = (label1->where == 'b'? min2: max2);
	/* first axis origin */

	if (!sf_getbool ("parallel1",&(axis1->parallel)) &&
	    !sf_getbool ("parallel", &(axis1->parallel))) axis1->parallel=true;
	/* tickmark labels parallel to the axis */

	if (NULL == (axis1->format = sf_getstring("format1")))
	    axis1->format = "%1.5g"; /* tick mark format */

	if (!sf_getint ("n1tic",&(axis1->ntic))) 
	    /*( n1tic axis1 number of ticmarks )*/
	    axis1->ntic = cube? 
		0.5 + inch1*(mid1-min1)/((max1-min1)*aspect*labelsz):
		0.5 + inch1/(aspect*labelsz);

	modify = (bool) (!sf_getfloat ("d1num", &(axis1->dnum)) ||
			 /*( d1num axis1 tic increment )*/
			 !sf_getfloat ("o1num", &(axis1->num0)));
	/*( o1num axis1 tic origin )*/
 
	axis1->ntic = vp_optimal_scale(axis1->ntic,
				       modify,
				       axis1->parallel,
				       axis1->format,
				       label1->min, label1->max, 
				       &(axis1->num0), 
				       &(axis1->dnum),
				       &(axis1->maxstrlen));
    }	
    
    if (label2 != NULL) { /* vertical axis */
	if (NULL == axis2)
	    axis2 = (struct Axis*) sf_alloc(1,sizeof(struct Axis));

	if (!sf_getfloat ("axisor2",&(axis2->orig)))
	    axis2->orig = (label2->where == 'l'? min1: max1);
	/* second axis origin */

	if (!sf_getbool ("parallel2",&(axis2->parallel)) &&
	    !sf_getbool ("parallel", &(axis2->parallel))) axis2->parallel=true;
	/* tickmark labels parallel to the axis */

	if (NULL == (axis2->format = sf_getstring("format2")))
	    axis2->format = "%1.5g"; /* tickmark format */

	if (!sf_getint ("n2tic",&(axis2->ntic)))
	    /*( n2tic axis2 number of ticmarks )*/
	    axis2->ntic = cube?
		0.5 + inch2*(mid2-min2)/((max2-min2)*aspect*labelsz):
		0.5 + inch2/(aspect*labelsz);

	modify = (bool) (!sf_getfloat ("d2num", &(axis2->dnum)) ||
			 /*( d2num axis2 tic increment )*/
			 !sf_getfloat ("o2num", &(axis2->num0)));
	/*( o2num axis2 tic origin )*/

	axis2->ntic = vp_optimal_scale(axis2->ntic,
				       modify,
				       axis2->parallel,
				       axis2->format,
				       label2->min, label2->max, 
				       &(axis2->num0), 
				       &(axis2->dnum),
				       &(axis2->maxstrlen));
    }

    if (label3 != NULL) {
	if (NULL == axis3)
	    axis3 = (struct Axis*) sf_alloc(1,sizeof(struct Axis));

	if (!sf_getfloat ("axisor3",&(axis3->orig)))
	    axis3->orig = min1;
	/* third axis origin */

	if (!sf_getbool ("parallel3",&(axis3->parallel)) &&
	    !sf_getbool ("parallel", &(axis3->parallel))) axis3->parallel=true;
	/* tickmark labels parallel to the axis */

	if (NULL == (axis3->format = sf_getstring("format3")))
	    axis3->format = "%1.5g"; /* tickmark format */
	
	if (!sf_getint ("n3tic",&(axis3->ntic)))
	    /*( n3tic axis3 number of ticmarks )*/
	    axis3->ntic = 0.5 + inch3/(aspect*labelsz);

	modify = (bool) (!sf_getfloat ("d3num", &(axis3->dnum)) ||
			 /*( d3num axis3 tic increment )*/
			 !sf_getfloat ("o3num", &(axis3->num0)));
	/*( o3num axis3 tic origin )*/
 
	axis3->ntic = vp_optimal_scale(axis3->ntic,
				       modify,
				       axis3->parallel,
				       axis3->format,
				       label3->min, label3->max, 
				       &(axis3->num0), 
				       &(axis3->dnum),
				       &(axis3->maxstrlen));
    }

    if (label4 != NULL) { /* vertical axis for out of plane */
	if (NULL == axis4)
	    axis4 = (struct Axis*) sf_alloc(1,sizeof(struct Axis));

	axis4->orig       = axis2->orig;
	axis4->parallel = axis2->parallel;
	axis4->format   = axis2->format;
	
	if (!sf_getint ("n4tic",&(axis4->ntic)))
	    /*( n4tic axis4 number of ticmarks )*/
	    axis4->ntic = 0.5+inch2*(max2-mid2)/((max2-min2)*aspect*labelsz);

	modify = (bool) (!sf_getfloat ("d4num", &(axis4->dnum)) ||
			 /*( d4num axis4 tic increment )*/
			 !sf_getfloat ("o4num", &(axis4->num0)));
	/*( o4num axis4 tic origin )*/

	axis4->ntic = vp_optimal_scale(axis4->ntic,
				       modify,
				       axis4->parallel,
				       axis4->format,
				       label4->min, label4->max, 
				       &(axis4->num0), 
				       &(axis4->dnum),
				       &(axis4->maxstrlen));
    }

    wheretics = (bool) ((NULL != (where = sf_getstring ("wheretics"))) &&
			('a' == *where));
    /*( wheretics where to put ticmarks )*/
}

static void make_grid (bool grid)
{
    bool need;
    float num;

    if (axis1 != NULL) { 
	if (!sf_getbool("grid",&need) && !sf_getbool("grid1",&need))
	    need = transp? false: grid;
	/*( grid1 to draw grid on first axis )*/

	if (need) {
	    grid1 = (struct Axis*) sf_alloc(1,sizeof(struct Axis));

	    if (!sf_getfloat ("g1num0",&(grid1->num0))) 
		grid1->num0 = axis1->num0;
	    /*( g1num0 grid mark origin on first axis )*/
	    grid1->orig = axis1->orig;

	    if (!sf_getfloat ("g1num",&(grid1->dnum))) {
		/*( g1num grid mark sampling on first axis )*/
		grid1->dnum = axis1->dnum;
		grid1->ntic = axis1->ntic;
	    } else {
		grid1->ntic=0; 
		if (xreverse) {
		    grid1->dnum = -grid1->dnum;
		    for (num=grid1->num0; num >= max1; num += grid1->dnum) {
			grid1->ntic++;
		    }
		} else {
		    for (num=grid1->num0; num <= max1; num += grid1->dnum) {
			grid1->ntic++;
		    }
		}
	    }
	}
    }

    if (axis2 != NULL) {
	if (!sf_getbool("grid",&need) && !sf_getbool("grid2",&need))
	    need = transp? grid: false;
	/*( grid2 to draw grid on second axis )*/

	if (need) {
	    grid2 = (struct Axis*) sf_alloc(1,sizeof(struct Axis));
	    
	    if (!sf_getfloat ("g2num0",&(grid2->num0))) 
		grid2->num0 = axis2->num0;
	    /*( g2num0 grid mark origin on second axis )*/
	    grid2->orig = axis2->orig;
	    
	    if (!sf_getfloat ("g2num",&(grid2->dnum))) {
		/*( g2num grid mark sampling on second axis )*/
		grid2->dnum = axis2->dnum;
		grid2->ntic = axis2->ntic;
	    } else {
		grid2->ntic=0; 
		if (yreverse) {
		    grid2->dnum = -grid2->dnum;
		    for (num=grid2->num0; num >= max2; num += grid2->dnum) {
			grid2->ntic++;
		    }
		} else {
		    for (num=grid2->num0; num <= max2; num += grid2->dnum) {
			grid2->ntic++;
		    }
		}
	    }
			      
	    
	} 
    }
    
    if ((axis3 != NULL) && 
	((sf_getbool("grid",&need) && need) || 
	 (sf_getbool("grid3",&need) && need))) {
	grid3 = (struct Axis*) sf_alloc(1,sizeof(struct Axis));
	/* to draw grid */

	grid3->num0 = axis3->num0;
	grid3->orig = axis3->orig;
	grid3->dnum = axis3->dnum;
	grid3->ntic = axis3->ntic;
    }

    if (NULL != grid1 || NULL != grid2 || NULL != grid3) {
	if (!sf_getint("gridcol",&gridcol)) gridcol=grid? VP_RED: framecol;
	/* grid color */
	if (!sf_getint("gridfat",&gridfat)) gridfat=1;
	/* grid fatness */
	if (!sf_getfloat("griddash",&griddash)) griddash=0.0f;
	/* grid dash pattern */
    } 
}

static void make_title (sf_file in, char wheret)
{
    bool want;
    char* where;
    float titlesz, vs, xc, yc;
    
    if (sf_getbool ("wanttitle",&want) && !want) {
	/* if include title */
	title = NULL;
	return;
    }

    if (NULL == title)
	title = (struct Label*) sf_alloc(1,sizeof(struct Label));

    if (!cube && NULL != (where = sf_getstring("wheretitle"))) {
	/* where to put title (top,bottom,left,right) */
	title->where = *where;
    } else {
	title->where = wheret;
    }

    if (!sf_getint ("titlefat",&(title->fat))) title->fat=0;
    /* title fatness */

    if (!sf_getfloat ("titlesz",&titlesz)) titlesz=10.;
    /* title font size */

    if (NULL == (title->text = sf_getstring("title")) &&
	NULL == (title->text = sf_histstring(in,"title")) &&
	NULL == (title->text = sf_histstring(in,"in"))) 
	title->text = blank;
    /*( title plot title )*/

    titlesz /= 33.;
    vs = titlesz * 0.6;

    if (title->where == 'l' || title->where == 'r') {	
	if (NULL != label2 && title->where == label2->where) {
	    if (label1->text != blank)
		vs += 3.25*labelsz; /* !!! fix that - use text justification */
	    else
                vs += 2.0*labelsz;
        }

	title->ypath = labelrot? -titlesz: titlesz;
	title->yup = 0.;
	title->xpath = 0.;
	title->xup = labelrot? titlesz: -titlesz;
	
    	yc = 0.5*(max2+min2);
	xc = (title->where == 'l')? min1: max1;
	vp_umove (xc, yc);
	vp_where (&xc, &yc);

	title->y = yc;
	title->x = (title->where == 'l')? xc-vs: xc+vs;
    } else {
	if (NULL != label1 && title->where == label1->where) 
	    vs += 3.25*labelsz; /* !!! fix that - use text justification */
	if (cube) vs += 2.00*labelsz;

	title->ypath = 0.;
	title->yup = titlesz;
	title->xpath = titlesz;
	title->xup = 0.;
	
	xc = 0.5*(max1+min1);
	yc = (title->where == 'b')? min2: max2;
	vp_umove (xc, yc);
	vp_where (&xc, &yc);
	
	title->x = xc;
	title->y = (title->where == 'b')? yc-vs: yc+vs;
    }	   
}    

static void make_barlabel (sf_file in)
{
    float vs;
    bool want;
    char* where, *labl, *unit;
    size_t len;

    if (sf_getbool ("wantbaraxis", &want) && !want) {
	/* if include scale bar axis */
	barlabel = NULL;
	return;
    }

    barlabel = (struct Label*) sf_alloc(1,sizeof(struct Label));

    if (!sf_getfloat ("barlabelsz",&barlabelsz)) {
	/* bar label font size */
	barlabelsz=labelsz;
    } else {
	barlabelsz /= 33.;
    }
    vs = 2.5*barlabelsz;
    
    if (!sf_getint ("barlabelfat",&(barlabel->fat)) &&
	!sf_getint ("labelfat",&(barlabel->fat))) barlabel->fat=0;
    /*( barlabelfat bar label fatness )*/
    if (NULL == (labl=sf_getstring("barlabel")) && 
	NULL == (labl=sf_histstring(in,"label"))) { 
	/*( barlabel bar label )*/
	barlabel->text = blank;
    } else if ((NULL == (unit=sf_getstring("barunit")) && 
		NULL == (unit=sf_histstring(in,"unit"))) ||
	       *unit == '\0' || (*unit == ' ' && *(unit+1) == '\0')) {
	/*( barunit bar unit )*/
	barlabel->text = labl;
    } else {
	len = strlen(labl)+strlen(unit)+4;
	barlabel->text = sf_charalloc(len);
	snprintf(barlabel->text,len,"%s (%s)",labl,unit);
	free(labl);
	free(unit);
    }
  
    if (NULL != (where = sf_getstring("wherebarlabel"))) {
	/* where to put bar label (top,bottom,left,right) */ 
	barlabel->where = *where;
	if (vertbar) {
	    if ('l' == barlabel->where) {
		barmin += 1.2*vs;
		barmax += 1.2*vs;
	    } else if ('r' != barlabel->where) {
		sf_error("For bartype=v, wherebarlabel must be 'l' or 'r'");
	    }
	} else {
	    if ('t' == barlabel->where) {
		barmin -= 1.2*vs;
		barmax -= 1.2*vs;
	    } else if ('b' != barlabel->where) {
		sf_error("For bartype=h, wherebarlabel must be 't' or 'b'");
	    }
	}
    } else {
	barlabel->where = vertbar? 'r':'b';
    }

    if (vertbar) {
	barlabel->ypath = labelrot? -barlabelsz: barlabelsz;
	barlabel->yup = 0.;
	barlabel->xpath = 0.;
	barlabel->xup = labelrot? barlabelsz: -barlabelsz;
	barlabel->y = orig2;	
	barlabel->x = (barlabel->where == 'l')? barmin-vs: barmax+vs;
    } else {
	barlabel->xpath = barlabelsz;
	barlabel->xup = 0.;
	barlabel->ypath = 0.;
	barlabel->yup = barlabelsz;
	barlabel->x = orig1;
	barlabel->y = (barlabel->where == 'b')? barmin-vs: barmax+vs;
    }
}

void vp_frame_init (sf_file in        /* input file */, 
		    const char* where /* three chars: label1, label2, title */,
		    bool grid         /* grid flag */)
/*< Initializing generic frame >*/
{
    make_labels(in,where[0],where[1]);
    make_axes();
    make_grid(grid);
    make_title(in,where[2]);
}

void vp_barframe_init (sf_file in, float min, float max)
/*< Initializing bar label frame >*/
{
    make_barlabel(in);
    make_baraxis(min,max);
}

void vp_plot_unset (void)
/*< return to defaults >*/
{
    vp_fat (gridfat);
    vp_color (gridcol);
    vp_set_dash (griddash);
    vp_clip(orig1-inch1,orig2-inch2,orig1+inch1,orig2+inch2);
}

void vp_simpleframe(void)
/*< Drawing simple frame >*/
{
    int i;
    float xc, yc, num;

    if (NULL != grid1 || NULL != grid2 || NULL != grid3) {
	vp_bgroup("grid");
	vp_fat (gridfat);
	vp_color (gridcol);
	vp_set_dash (griddash);

	if (NULL != grid1) {
	    for (i=0; i < grid1->ntic; i++) {
		num = grid1->num0 + i*(grid1->dnum);
		if (fabsf(grid1->dnum) > SF_EPS && 
		    fabsf(num) < SF_EPS) num=0.;
		
		if (cube) {
		    xc = (num-label1->min)*(mid1-min1)/
			(label1->max-label1->min);
		} else {
		    xc = num;
		}
		
		vp_umove (xc, min2);
		vp_udraw (xc, max2);
	    }
	}

	if (NULL != grid2) {
	    for (i=0; i < grid2->ntic; i++) {
		num = grid2->num0 + i*(grid2->dnum);
		if (fabsf(grid2->dnum) > SF_EPS && 
		    fabsf(num) < SF_EPS) num=0.;

		if (cube) {
		    yc = mid2+(num-label2->max)*(mid2-min2)/
			(label2->max-label2->min);
		} else {
		    yc = num;
		}	    

		vp_umove (min1, yc);
		vp_udraw (max1, yc);
	    }
	}

	vp_set_dash (0);
	vp_egroup();
    }

    /* draw outline */   
    vp_bgroup("frame");
    vp_color(framecol);
    vp_umove(min1, min2);

    if (cube) {
	if (flat) {
	    vp_udraw(max1,min2);
	    vp_udraw(max1,mid2);
	} else {
	    vp_udraw(mid1,min2);
	    vp_udraw(max1,max2-mid2);
	    vp_udraw(max1,max2);
	    vp_udraw(mid1,mid2);
	}
	vp_udraw(min1,mid2);

	vp_umove(min1,min2);
	if (flat) {
	    vp_udraw(min1,max2);
	    vp_udraw(mid1,max2);
	} else {
	    vp_udraw(min1,mid2);
	    vp_udraw(max1-mid1,max2);
	    vp_udraw(max1,max2);
	    vp_udraw(mid1,mid2);
	}
	vp_udraw(mid1,min2);
    } else { /* simple box */
	vp_udraw(min1, max2);
	vp_udraw(max1, max2);
	vp_udraw(max1, min2);
	vp_udraw(min1, min2);
    }

    vp_egroup();
}

void vp_framenum(float num)
/*< Outputting frame number >*/
{
    char string[80];
    float x, y;

    vp_bgroup("frame number");
    vp_umove(min1,min2);
    vp_where(&x,&y);

    sprintf (string, "%g", num);
    vp_tjust (TH_CENTER, TV_TOP);
    if ((NULL != title  &&  title->where == 'b') || 
	(NULL != label1 && label1->where == 'b'))
        vp_gtext (x, y-4.*labelsz, labelsz, 0., 0., labelsz, string);
    else
        vp_gtext (x+3.*labelsz, y-2.*labelsz, labelsz, 0., 0., labelsz, string);
    vp_egroup();
}

void vp_simplebarframe (void)
/*< Drawing simple frame for the bar label >*/
{
    float min, max;

    vp_bgroup("scalebar frame");
    vp_color(framecol);
    vp_fat (gridfat);

    if (vertbar) {
	min = orig2-0.5*inch2;
	max = orig2+0.5*inch2;

	vp_clip(barmin-0.5*inch1,min-0.5*inch2,
		barmax+0.5*inch1,max+0.5*inch2);

	vp_move(barmin,min);
	vp_draw(barmax,min);
	vp_draw(barmax,max);
	vp_draw(barmin,max);
	vp_draw(barmin,min);
    } else {
	min = orig1-0.5*inch1;
	max = orig1+0.5*inch1;

	vp_clip(min-0.5*inch1,barmin-0.5*inch2,
		max+0.5*inch1,barmax+0.5*inch2);

	vp_move(min,barmin);
	vp_draw(min,barmax);
	vp_draw(max,barmax);
	vp_draw(max,barmin);
	vp_draw(min,barmin);
    }

    vp_egroup();
}

void vp_frame(void)
/*< Drawing frame >*/
{
    int i;
    bool need;
    float num, xc, yc, vs, xp[4], yp[4], sx, sy;
    char string[32];

    vp_bgroup("plot frame");

    if (cube) {
	vp_bgroup("corner");
	if (flat) {
	    /* remove rectange */
	    vp_color(VP_BLACK);
	    xp[0] = mid1; yp[0] = mid2;
	    xp[1] = max1; yp[1] = mid2;
	    xp[2] = max1; yp[2] = max2;
	    xp[3] = mid1; yp[3] = max2;    
	    vp_ufill(xp,yp,4);
	} else {
	    /* remove two triangles */
	    vp_color(VP_BLACK);
	    xp[0] = min1-1;    yp[0] = mid2;
	    xp[1] = min1-1;    yp[1] = max2+1;
	    xp[2] = max1-mid1; yp[2] = max2+1;
	    vp_ufill(xp,yp,3);
	    xp[0] = mid1;      yp[0] = min2-1;
	    xp[1] = max1+1;    yp[1] = min2-1;
	    xp[2] = max1+1;    yp[2] = max2-mid2;
	    vp_ufill(xp,yp,3);
	}
	vp_egroup();
    }

    vp_simpleframe();
    
    if (NULL != label1) {
	vp_bgroup("axis 1");
	vp_fat (label1->fat);

	/* plot label */
	if (label1->where == 't') { 
	    vp_tjust (TH_CENTER, TV_BOTTOM);
	    vs = tickscale1*labelsz;
	} else {
	    vp_tjust (TH_CENTER, TV_TOP);
	    vs = -tickscale1*labelsz;
	}

	sy = label1->y;
	if (! axis1->parallel) 
	    sy += 2*(axis1->maxstrlen-1)*vs;

	vp_gtext(label1->x, sy, 
		 label1->xpath, label1->ypath, 
		 label1->xup, label1->yup, label1->text);

	/* plot tics */

	if (axis1->parallel) {
	    if (label1->where == 't') { 
		vp_tjust (TH_CENTER, TV_BASE);
	    } else {
		vp_tjust (TH_CENTER, TV_CAP);
	    }
	} else {
	    if (((label1->where != 'b') && labelrot) ||
		((label1->where == 'b') && !labelrot)) {
		vp_tjust (TH_RIGHT, TV_HALF);
	    } else {
		vp_tjust (TH_LEFT, TV_HALF);
	    }
	}

	for (i=0; i < axis1->ntic; i++) {
	    num = axis1->num0 + i*(axis1->dnum);
	    if (fabsf(axis1->dnum) > SF_EPS && 
		fabsf(num) < SF_EPS) num=0.;
	    
	    if (cube) {
		xc = (num-label1->min)*(mid1-min1)/(label1->max-label1->min);
		yc = min2;
	    } else {
		xc = num;
		yc = axis1->orig;
	    }

	    vp_umove (xc, yc);
	    vp_where (&xc, &yc);
	    vp_draw (xc, yc+vs);

	    snprintf (string,32,axis1->format, num);

	    if (axis1->parallel) {
		vp_gtext(xc, yc+ticksep*vs, 
			 labelsz, 0., 
			 0., labelsz, string);
	    } else {
		if (labelrot) {
		    vp_gtext(xc, yc+ticksep*vs, 
			     0., -labelsz, 
			     labelsz, 0., string);
		} else {
		    vp_gtext(xc, yc+ticksep*vs, 
			     0., labelsz, 
			     -labelsz, 0., string);
		}
	    }
	}
	vp_egroup();
    }    

    if (NULL != label2) {
	vp_bgroup("axis 2");
	vp_fat (label2->fat);

	/* plot label */
	if (((label2->where != 'l') && labelrot) ||
	    ((label2->where == 'l') && !labelrot)) {
	    vp_tjust (TH_CENTER, TV_BOTTOM);
	} else {
	    vp_tjust (TH_CENTER, TV_TOP);
	}

	vs = label2->where == 'l'? -tickscale2*labelsz: tickscale2*labelsz;

	sx = label2->x;
	if (! axis2->parallel) 
	    sx += 2*(axis2->maxstrlen-1)*vs;

	vp_gtext(sx, label2->y, 
		 label2->xpath, label2->ypath, 
		 label2->xup, label2->yup, label2->text);
	
        /* plot tics */

	if (axis2->parallel) {
	    if (((label2->where != 'l') && labelrot) ||
		((label2->where == 'l') && !labelrot)) {
		vp_tjust (TH_CENTER, TV_BASE);
	    } else {
		vp_tjust (TH_CENTER, TV_CAP);
	    }
	} else {
	    if (label2->where == 'l') {
		vp_tjust (TH_RIGHT, TV_HALF);
	    } else {
		vp_tjust (TH_LEFT, TV_HALF);
	    }
	}

	for (i=0; i < axis2->ntic; i++) {
	    num = axis2->num0 + i*(axis2->dnum);
	    if (fabsf(axis2->dnum) > SF_EPS && 
		fabsf(num) < SF_EPS) num=0.;

	    if (cube) {
		yc = mid2+(num-label2->max)*(mid2-min2)/
		    (label2->max-label2->min);
		xc = min1;
	    } else {
		yc = num;
		xc = axis2->orig;
	    }	    

	    /* draw the tick */
	    vp_umove (xc, yc);
	    vp_where (&xc, &yc);
	    vp_draw (xc+vs, yc);

	    /* tick label */
	    snprintf (string,32,axis2->format, num);

	    if (axis2->parallel) {
		if (labelrot) {
		    vp_gtext(xc+ticksep*vs, yc, 
			     0., -labelsz, 
			     labelsz, 0., string);
		} else {
		    vp_gtext(xc+ticksep*vs, yc, 
			     0., labelsz, 
			     -labelsz, 0., string);
		}
	    } else {
		vp_gtext(xc+ticksep*vs, yc, labelsz, 0., 0., labelsz, string);
	    }
	}
	vp_egroup();
    }

    if (NULL != label3) {
	vp_bgroup("axis 3");
	vp_fat (label3->fat);

	/* plot label */

	vp_tjust (TH_CENTER, TV_TOP);
	vs = -tickscale3*labelsz;

	sx = label3->x;
	sy = label3->y;
	if (! axis3->parallel) {
	    if (flat) {
		sy += 2*(axis3->maxstrlen-1)*vs;
	    } else {
		sy += 2*(axis3->maxstrlen-1)*vs*costh;	
		sx -= 2*(axis3->maxstrlen-1)*vs*sinth;
	    }
	}
	
	vp_gtext(sx, sy, 
		 label3->xpath, label3->ypath, 
		 label3->xup, label3->yup, label3->text);
	
        /* plot tics */
	if (axis3->parallel) {
	    vp_tjust (TH_CENTER, TV_CAP);
	} else {
	    if (labelrot) {
		vp_tjust (TH_LEFT, TV_HALF);
	    } else {
		vp_tjust (TH_RIGHT, TV_HALF);
	    }
	}

	for (i=0; i < axis3->ntic; i++) {
	    num = axis3->num0 + i*(axis3->dnum);
	    if (fabsf(axis3->dnum) > SF_EPS && 
		fabsf(num) < SF_EPS) num=0.;

	    if (flat) {
		xc = mid1+(num-label3->min)*(max1-mid1)/
		    (label3->max-label3->min);
		yc = min2;
	    
		vp_umove (xc, yc);
		vp_where (&xc, &yc);
		vp_draw (xc, yc+vs);

		snprintf (string,32,axis3->format, num);

		if (axis3->parallel) {
		    vp_gtext(xc, yc+ticksep*vs,
			     labelsz,0.,0.,labelsz,string);
		} else {
		    if (labelrot) {
			vp_gtext(xc, yc+ticksep*vs, 
				 0., -labelsz, 
				 labelsz, 0., string);
		    } else {
			vp_gtext(xc, yc+ticksep*vs, 
				 0., labelsz, 
				 -labelsz, 0., string);
		    }
		}

	    } else { /* not flat */

		xc = mid1+(num-label3->min)*(max1-mid1)/
		    (label3->max-label3->min);
		yc = min2+(num-label3->min)*(max2-mid2-min2)/
		    (label3->max-label3->min);
	    
		vp_umove (xc, yc);
		vp_where (&xc, &yc);
		vp_draw (xc-vs*sinth, yc+vs*costh);
		
		snprintf (string,32,axis3->format, num);

		if (axis3->parallel) {
		    vp_gtext(xc-ticksep*vs*sinth, 
			     yc+ticksep*vs*costh, 
			     labelsz*costh,labelsz*sinth,
			     -labelsz*sinth,labelsz*costh,string);
		} else {
		    if (labelrot) {
			vp_gtext(xc-ticksep*vs*sinth, 
				 yc+ticksep*vs*costh,
				 labelsz*sinth, -labelsz*costh, 
				 labelsz*costh,  labelsz*sinth, string);
		    } else {
			vp_gtext(xc-ticksep*vs*sinth, 
				 yc+ticksep*vs*costh,
				 -labelsz*sinth, labelsz*costh, 
				 -labelsz*costh, -labelsz*sinth, string);
		    }
		}
	    }
	}
	vp_egroup();
    }

    if (NULL != label4) {
	vp_bgroup("axis 4");
	vp_fat (label4->fat);

	/* plot label */
	if (((label4->where != 'l') && labelrot) ||
	    ((label4->where == 'l') && !labelrot)) {
	    vp_tjust (TH_CENTER, TV_BOTTOM);
	} else {
	    vp_tjust (TH_CENTER, TV_TOP);
	}

	vs = label4->where == 'l'? -tickscale4*labelsz: tickscale4*labelsz;

	sx = label4->x;
	if (! axis4->parallel) 
	    sx += 2*(axis4->maxstrlen-1)*vs;

	vp_gtext(sx, label4->y, 
		 label4->xpath, label4->ypath, 
		 label4->xup, label4->yup, label4->text);
	
        /* plot tics */

	if (axis4->parallel) {
	    if (((label4->where != 'l') && labelrot) ||
		((label4->where == 'l') && !labelrot)) {
		vp_tjust (TH_CENTER, TV_BASE);
	    } else {
		vp_tjust (TH_CENTER, TV_CAP);
	    }
	} else {
	    if (label4->where == 'l') {
		vp_tjust (TH_RIGHT, TV_HALF);
	    } else {
		vp_tjust (TH_LEFT, TV_HALF);
	    }
	}

	for (i=0; i < axis4->ntic; i++) {
	    num = axis4->num0 + i*(axis4->dnum);
	    if (fabsf(axis4->dnum) > SF_EPS && 
		fabsf(num) < SF_EPS) num=0.;

	    yc = mid2+(num-label4->min)*(max2-mid2)/
		(label4->max-label4->min);
	    xc = min1;	    

	    /* draw the tick */
	    vp_umove (xc, yc);
	    vp_where (&xc, &yc);
	    vp_draw (xc+vs, yc);

	    /* tick label */
	    snprintf (string,32,axis2->format, num);

	    if (axis4->parallel) {
		if (labelrot) {
		    vp_gtext(xc+ticksep*vs, yc, 
			     0., -labelsz, 
			     labelsz, 0., string);
		} else {
		    vp_gtext(xc+ticksep*vs, yc, 
			     0., labelsz, 
			     -labelsz, 0., string);
		}
	    } else {
		vp_gtext(xc+ticksep*vs, yc, labelsz, 0., 0., labelsz, string);
	    }
	}
	vp_egroup();
    }

    if (NULL != title) {
	vp_bgroup("title");
	vp_fat (title->fat);

	if ((title->where == 't') || 
	    ((title->where == 'r') && labelrot) ||
	    ((title->where == 'l') && !labelrot)) {
	    vp_tjust (TH_CENTER, TV_BOTTOM);
	} else {
	    vp_tjust (TH_CENTER, TV_TOP);
	} 

	vp_gtext(title->x, title->y, 
		 title->xpath, title->ypath, 
		 title->xup, title->yup, title->text);
	vp_egroup();
    }
    
    if (cube) {
	/* draw colored lines */
	vs = tickscale*labelsz;

	yc = mid2-(mid2-min2)*(frame1+0.5)*d1/(l2min-l2max);

	vp_bgroup("colored lines");
	vp_color(cubelinecol);
	vp_umove(min1,yc);
	
	if (flat) {
	    vp_udraw(max1,yc);
	} else {
	    vp_udraw(mid1,yc);
	    vp_udraw(max1,yc+max2-mid2);
	}

	if (!sf_getbool("framelabel",&need) && 
	    !sf_getbool("framelabel1",&need))
	    need = (bool) (NULL != label1);
	/* to put numbers at frame ends */

	if (need) {
	    vp_color(framelabelcol);
	    vp_where(&xc,&yc);
	    
	    num = l2max+(frame1+0.5)*d1;
	    if (fabsf(d1) > SF_EPS && 
		fabsf(num) < 0.001*fabsf(d1)) num=0.;
	    snprintf (string,32,"%1.5g",num);
	    
	    if (NULL != label2) {
		if (axis2->parallel) {
		    if (labelrot) {
			vp_tjust (TH_CENTER, TV_BASE);
			vp_gtext(xc+0.5*vs,yc,0.,-labelsz,labelsz,0.,string);
		    } else {
			vp_tjust (TH_CENTER, TV_CAP);
			vp_gtext(xc+0.5*vs,yc,0.,labelsz,-labelsz,0.,string);
		    }
		} else {
		    vp_tjust (TH_LEFT, TV_HALF);
		    vp_gtext(xc+0.5*vs, yc, labelsz, 0., 0., labelsz, string);
		}
	    }
	}

	xc = min1+(mid1-min1)*(frame2+0.5)*d2/(l1max-l1min);

	vp_color(cubelinecol);
	vp_umove(xc,min2);
	
	if (flat) {
	    vp_udraw(xc,max2);
	} else {
	    vp_udraw(xc,mid2);
	    vp_udraw(xc+max1-mid1,max2);
	}

	if (!sf_getbool("framelabel",&need) && 
	    !sf_getbool("framelabel2",&need))
	    need = (bool) (NULL != label2);
	/* to put numbers at frame ends */
	 
	if (need) {
	    vp_color(framelabelcol);
	    vp_where(&xc,&yc);
	    
	    num = l1min+(frame2+0.5)*d2;
	    if (fabsf(d2) > SF_EPS && 
		fabsf(num) < SF_EPS) num=0.;
	    snprintf (string,32,"%1.5g",num);
	    
	    vp_tjust (TH_CENTER, TV_BOTTOM);
	    vp_gtext(xc,yc+0.5*vs, labelsz, 0., 0., labelsz, string);
	}	    

	if (!sf_getbool("framelabel",&need) && 
	    !sf_getbool("framelabel3",&need))
	    need = (bool) (NULL != label3);
	/* to put numbers at frame ends */

	if (need) {
	    num = label3->min+(frame3+0.5)*d3;
	    if (fabsf(d3) > SF_EPS && 
		fabsf(num) < SF_EPS) num=0.;
	    snprintf (string,32,"%1.5g",num);
	}

	yc = mid2+(max2-mid2)*(frame3+0.5)*d3/(l3max-l3min);
	xc = mid1+(max1-mid1)*(frame3+0.5)*d3/(l3max-l3min);

	vp_color(cubelinecol);	
	if (flat) {
	    vp_umove(xc,min2);
	    vp_udraw(xc,mid2);
	} else {
	    vp_umove(xc,yc+min2-mid2);
	    vp_udraw(xc,yc); 
	}

	if (need && flat) {
	    vp_color(framelabelcol);
	    vp_where(&xc,&yc);

	    vp_tjust (TH_CENTER, TV_BOTTOM);
	    vp_gtext(xc,yc+0.5*vs, labelsz, 0., 0., labelsz, string);
	}

	yc = mid2+(max2-mid2)*(frame3+0.5)*d3/(l3max-l3min);
	xc = mid1+(max1-mid1)*(frame3+0.5)*d3/(l3max-l3min);

	vp_color(cubelinecol);
	if (flat) {
	    vp_umove(min1,yc);
	    vp_udraw(mid1,yc);
	} else {
	    vp_umove(xc,yc);
	    vp_udraw(xc+min1-mid1,yc);
	}

	if (need) {
	    vp_color(framelabelcol);
	    vp_where(&xc,&yc);

	    if (flat) {
		if (axis2->parallel) {
		    if (labelrot) {
			vp_tjust (TH_CENTER, TV_BASE);
			vp_gtext(xc+0.5*vs,yc,0.,-labelsz,labelsz,0.,string);
		    } else {
			vp_tjust (TH_CENTER, TV_CAP);
			vp_gtext(xc+0.5*vs,yc,0.,labelsz,-labelsz,0.,string);
		    }
		} else {
		    vp_tjust (TH_LEFT, TV_HALF);
		    vp_gtext(xc+0.5*vs, yc, labelsz, 0., 0., labelsz, string);
		}
	    } else {
		vp_gtext(xc-0.5*vs*sinth, yc+0.5*vs*costh, 
			 labelsz*costh,labelsz*sinth,
			 -labelsz*sinth,labelsz*costh,string);
	    }
	} /* label3 */

	vp_egroup();
    } /* if cube */
    
    vp_uclip (min1, min2, max1, max2);
    vp_egroup();
}

void vp_barframe(void)
/*< Drawing bar label frame >*/
{
    int i;
    float num, xc, yc, vs;
    char string[32];

    vp_simplebarframe();
    
    if (NULL != barlabel) {
	vp_fat (barlabel->fat);

	if (vertbar) {
	    if (((barlabel->where != 'l') && labelrot) ||
		((barlabel->where == 'l') && !labelrot)) { 
		vp_tjust (TH_CENTER, TV_BOTTOM);
	    } else {
		vp_tjust (TH_CENTER, TV_TOP);
	    }
	    vs = barlabel->where == 'l'? -tickscale1*labelsz: tickscale1*labelsz;
	    if (! baraxis->parallel) 
		barlabel->x += 2*(baraxis->maxstrlen-1)*vs;
	} else {
	    if (barlabel->where == 't') {
		vp_tjust (TH_CENTER, TV_BOTTOM);
		vs = tickscale2*labelsz;
	    } else {
		vp_tjust (TH_CENTER, TV_TOP);
		vs = -tickscale2*labelsz;
	    }
	    if (! baraxis->parallel) 
		barlabel->y += 2*(baraxis->maxstrlen-1)*vs;
	}

	vp_gtext(barlabel->x, barlabel->y, 
		 barlabel->xpath, barlabel->ypath, 
		 barlabel->xup, barlabel->yup, barlabel->text);

	/* plot tics */

	if (vertbar) {
	    if (baraxis->parallel) {
		if (((barlabel->where != 'l') && labelrot) ||
		    ((barlabel->where == 'l') && !labelrot)) {
		    vp_tjust (TH_CENTER, TV_BASE);
		} else {
		    vp_tjust (TH_CENTER, TV_CAP);
		}
	    } else {
		if (barlabel->where == 'l') {
		    vp_tjust (TH_RIGHT, TV_HALF);
		} else {
		    vp_tjust (TH_LEFT, TV_HALF);
		}
	    }
	} else {
	    if (baraxis->parallel) {
		if (barlabel->where == 't') {
		    vp_tjust (TH_CENTER, TV_BASE);
		} else {
		    vp_tjust (TH_CENTER, TV_CAP);
		}
	    } else {
		if (((barlabel->where != 'b') && labelrot) ||
		    ((barlabel->where == 'b') && !labelrot)) {
		    vp_tjust (TH_RIGHT, TV_HALF);
		} else {
		    vp_tjust (TH_LEFT, TV_HALF);
		}
	    }
	}

	for (i=0; i < baraxis->ntic; i++) {
	    num = baraxis->num0 + i*(baraxis->dnum);
	    if (fabsf(baraxis->dnum) > SF_EPS && 
		fabsf(num) < SF_EPS) num=0.;

	    snprintf (string,32,baraxis->format, num);	    
	    
	    if (vertbar) {
		yc = bar0 + i*dbar;
		xc = baraxis->orig;
		
		vp_move (xc, yc);
		vp_draw (xc+vs, yc);

		if (baraxis->parallel) {
		    if (labelrot) {
			vp_gtext(xc+ticksep*vs, yc, 
				 0., -labelsz, 
				 labelsz, 0., string);
		    } else {
			vp_gtext(xc+ticksep*vs, yc, 
				 0., labelsz, 
				 -labelsz, 0., string);
		    }
		} else {
		    vp_gtext(xc+ticksep*vs, yc, 
			     labelsz, 0., 
			     0., labelsz, string);
		}
	    } else { /* horizontal bar */
		xc = bar0 + i*dbar;
		yc = baraxis->orig;

		vp_move (xc, yc);
		vp_draw (xc, yc+vs);

		if (baraxis->parallel) {
		    vp_gtext(xc, yc+ticksep*vs, 
			     labelsz, 0., 
			     0., labelsz, string);
		} else {
		    if (labelrot) {
			vp_gtext(xc, yc+ticksep*vs, 
				 0., -labelsz, 
				 labelsz, 0., string);
		    } else {
			vp_gtext(xc, yc+ticksep*vs, 
				 0., labelsz, 
				 -labelsz, 0., string);
		    }
		}
	    }
	}
    }
}

void vp_barraster (int nbuf, unsigned char** buf /* buf[1][nbuf] */)
/*< Filling bar label with rasters >*/
{
    vp_barframe();
    if (vertbar) {
	vp_raster(buf, false, 256, nbuf, 1, 
		  barmin,orig2-0.5*inch2,barmax,orig2+0.5*inch2, 3);
    } else {
	vp_raster(buf, false, 256, nbuf, 1, 
		  orig1-0.5*inch1,barmin,orig1+0.5*inch1,barmax, 0);
    }
    vp_simplebarframe();
}

void vp_cuberaster(int n1, int n2, 
		   unsigned char** buf      /* buf[n2][n1] */ , 
		   int f1, int f2, int f3   /* frame numbers */) 
/*< Filling 3-D cube with rasters >*/
{
    vp_uraster (buf, false, 256, n1, n2, min1,min2,max1,max2, 3);
    vp_cubeframe(f1,f2,f3);
}

void vp_cubeframe(float f1, float f2, float f3   /* frame numbers */) 
/*< Drawing 3-D frame >*/
{
    frame1 = f1;
    frame2 = f2;
    frame3 = f3;
    vp_frame();    
}

void vp_barline (int nc     /* number of contours */, 
		 float *c   /* contours [nc] */, 
		 float cmin /* minimum value */, 
		 float cmax /* maximum value */)
/*< Drawing contour lines in a bar label >*/
{
    int ic;
    float level, c0, dc, ci;

    c0 = 0.5*(cmax+cmin);
    dc = cmax-cmin+SF_EPS;

    vp_barframe();
    for (ic = 0; ic < nc; ic++) {
	ci = (c[ic]-c0)/dc;
	if (fabsf(ci) > 0.5) continue;
	
	vp_plot_set (ic);
	if (vertbar) {
	    level = orig2 + inch2*ci;
	    vp_move(barmin,level);
	    vp_draw(barmax,level);
	} else {
	    level = orig1 + inch1*ci;
	    vp_move(level,barmin);
	    vp_draw(level,barmax);
	}
    }
    /*   vp_simplebarframe(); */
}

/* 	$Id$	 */
