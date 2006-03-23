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

#include <float.h>
#include <math.h>
#include <stdio.h>

#include <rsf.h>
/*^*/

#include "stdplot.h"
#include "axis.h"
#include "vplot.h"
#include "plot.h"

static float min1,min2, max1,max2, mid1,mid2, inch1,inch2, orig1,orig2, inch3;
static float labelsz, barlabelsz, barmin, barmax, bar0, dbar, sinth, costh;
static float d1, d2, d3, frame1;
static int framecol, frame2, frame3, gridcol, gridfat=1;
static int cubelinecol=VP_WHITE;
static bool labelrot, transp, wheretics, scalebar, vertbar, wherebartics;
static bool cube=false, flat;
static char blank[]=" ";
static const float aspect=0.8;

static struct Label {
    char where;
    int fat;
    float x, y, xpath, ypath, xup, yup, min, max;
    char* text;
}   /*@null@*/ *label1=NULL, 
    /*@null@*/ *label2=NULL, 
    /*@null@*/ *label3=NULL, 
    /*@null@*/ *title=NULL, 
    /*@null@*/ *barlabel=NULL;

static struct Axis {
    int ntic;
    float or, dnum, num0;
}   /*@null@*/ *axis1=NULL, 
    /*@null@*/ *axis2=NULL, 
    /*@null@*/ *axis3=NULL, 
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
    bool pad, set, xreverse, yreverse;
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
	    vertbar = (bartype[0] == 'v');
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
		yll -= 0.4*barwd;
		yur -= (0.04*screenht + barwd);
		barmax = 0.98*screenht;
		barmin = barmax-barwd;
	    } else {
		yur += 0.4*barwd;
		yll += (0.04*screenht + barwd);
		barmin = 0.12*screenht;
		barmax = barmin+barwd;
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

    vp_coordinates();
}

void vp_cubecoord (bool front /* front or side */,
		   float umin1, float umax1, float umin2, float umax2)
/*< setup the coordinate system for the cube front >*/
{
    float scale1, scale2;

    /* find scales and user origin */
    if (front) {
	scale1 = inch1 * (mid1 - min1) / (max1 - min1) / (umax1 - umin1);
	vp_orig (orig1-0.5*inch1,orig2-0.5*inch2);
    } else {
	scale1 = inch1 * (max1 - mid1) / (max1 - min1) / (umax1 - umin1);
	vp_orig (orig1+inch1*((mid1 - min1)/(max1 - min1)-0.5),orig2-0.5*inch2);
    }
    
    scale2 = inch2 * (mid2 - min2) / (max2 - min2) / (umax2 - umin2);

    /* setup the coordinate system */
    vp_scale (scale1, scale2);
    vp_uorig (umin1, umin2);
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
		       bool flat1                 /* flat flag */) 
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

    if (sf_getbool ("wantaxis", &want) && !want) {
	/* if draw axes */
	label1 = NULL;
	label2 = NULL;
	if (cube) label3 = NULL;
	return;
    }

    if (sf_getbool ("wantaxis1",&want) && !want) {
	/* if draw first axis */
	label1 = NULL;
    } else if (NULL == label1) { 
	label1 = (struct Label*) sf_alloc(1,sizeof(struct Label));
    }

    if (sf_getbool ("wantaxis2",&want) && !want) {
	/* if draw second axis */
	label2 = NULL;
	if (!cube && NULL == label1) return;
    } else if (NULL == label2) {
	label2 = (struct Label*) sf_alloc(1,sizeof(struct Label));
    }

    if (cube) {
	if (sf_getbool ("wantaxis3",&want) && !want) {
	    /* if draw third axis (in cube plots) */
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

    if (!sf_getfloat ("labelsz",&labelsz)) labelsz=8.;
    /* label size */
    if (cube) {
	labelsz *= 0.03; /* slightly smaller */
    } else {
	labelsz /= 33.;
    }
    vs = 2.5*labelsz;

    if (NULL != label1) {
	if (!cube && NULL != (where = sf_getstring("wherexlabel"))) {
	    /* where to put horizontal axis (top,bottom) */
	    label1->where = *where;
	} else {
	    label1->where = where1;
	}
	if (!sf_getint ("labelfat",&(label1->fat))) label1->fat=0;
	/* label fatness */
	if ((NULL == (labl=sf_getstring(transp? "label2":"label1"))) &&
	    (NULL == (labl=sf_histstring(in,
					 transp? "label2":"label1")))) {
	    label1->text = blank;
	} else if ((NULL == (unit=sf_getstring(transp? "unit2":"unit1"))) &&
		   (NULL == (unit=sf_histstring(in,
						transp? "unit2":"unit1")))) {
	    label1->text = labl;
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
	    if (!sf_histint(in,"n2",&n2)) n2=1;
	    if (!sf_histfloat(in,"o2",&o2)) o2=0.;
	    if (!sf_histfloat(in,"d2",&d2)) d2=1.;
	    label1->min = o2-0.5*d2;
	    label1->max = o2+(n2-0.5)*d2;
	} else {
	    label1->min = min1;
	    label1->max = max1;
	}
    }

    if (NULL != label3) {
	if (!sf_getint ("labelfat",&(label3->fat))) label3->fat=0;
	/* label fatness */
	if ((NULL == (labl=sf_getstring("label3"))) &&
	    (NULL == (labl=sf_histstring(in,"label3")))) {
	    label3->text = blank;
	} else if ((NULL == (unit=sf_getstring("unit3"))) &&
		   (NULL == (unit=sf_histstring(in,"unit3")))) {
	    label3->text = labl;
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

	n3 = sf_leftsize(in,2);
	if (!sf_histfloat(in,"o3",&o3)) o3=0.;
	if (!sf_histfloat(in,"d3",&d3)) d3=1.;
	label3->min = o3-0.5*d3;
	label3->max = o3+(n3-0.5)*d3;
    }

    if (NULL != label2) {
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
	} else if ((NULL == (unit=sf_getstring(transp? "unit1":"unit2"))) &&
		   (NULL == (unit=sf_histstring(in,
						transp? "unit1":"unit2")))) {
	    label2->text = labl;
	} else {
	    len = strlen(labl)+strlen(unit)+4;
	    label2->text = sf_charalloc(len);
	    snprintf(label2->text,len,"%s (%s)",labl,unit);
	    free(labl);
	    free(unit);
	}

	if (!sf_getbool ("labelrot",&labelrot)) labelrot = false;
	/* if rotate vertical label */
	if (labelrot) vs *= 1.7;

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
	    if (!sf_histint(in,"n1",&n1)) n1=1;
	    if (!sf_histfloat(in,"o1",&o1)) o1=0.;
	    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
	    label2->max = o1-0.5*d1;
	    label2->min = o1+(n1-0.5)*d1;
	} else {
	    label2->min = min2;
	    label2->max = max2;
	}
    } 
}

static void make_baraxis (float min, float max)
{
    char* where;

    if (barlabel == NULL) return;

    if (NULL == baraxis) 
	baraxis = (struct Axis*) sf_alloc(1,sizeof(struct Axis));

    if (!sf_getfloat ("axisorbar",&(baraxis->or))) {
	/* bar origin */
	if (vertbar) {
	    baraxis->or = (barlabel->where == 'l'? barmin: barmax);
	} else {
	    baraxis->or = (barlabel->where == 'b'? barmin: barmax);
	}
    }

    if (!sf_getint ("nbartic",&(baraxis->ntic))) baraxis->ntic = 1;
    /* bar ticmarks */
    if (!sf_getfloat ("dbarnum", &(baraxis->dnum)) ||
	!sf_getfloat ("obarnum", &(baraxis->num0))) 
	baraxis->ntic = vp_optimal_scale((vertbar? inch2: inch1)/
					 (aspect*labelsz), 
					 min, max, 
					 &(baraxis->num0), 
					 &(baraxis->dnum));
    bar0 = (baraxis->num0-0.5*(min+max))/(max-min+FLT_EPSILON);
    dbar = (baraxis->dnum)/(max-min+FLT_EPSILON);

    if (vertbar) {
	bar0 = orig2 + inch2*bar0;
	dbar *= inch2;
    } else {
	bar0 = orig1 + inch1*bar0;
	dbar *= inch1;
    }

    wherebartics = (NULL != (where = sf_getstring ("wherebartics"))) &&
	('a' == *where);
}	

static void make_axes (void)
{
    char* where;

    if (label1 != NULL) {
	if (NULL == axis1) 
	    axis1 = (struct Axis*) sf_alloc(1,sizeof(struct Axis));

	if (!sf_getfloat ("axisor1",&(axis1->or)))
	    axis1->or = (label1->where == 'b'? min2: max2);
	/* first axis origin */

	if (!sf_getint ("n1tic",&(axis1->ntic))) axis1->ntic = 1;
	/* first axis ticmarks */
	if (!sf_getfloat ("d1num", &(axis1->dnum)) ||
	    !sf_getfloat ("o1num", &(axis1->num0))) 
	    axis1->ntic = vp_optimal_scale(cube? 
					   inch1*(mid1-min1)/
					   ((max1-min1)*aspect*labelsz):
					   inch1/(aspect*labelsz), 
					   label1->min, label1->max, 
					   &(axis1->num0), 
					   &(axis1->dnum));
    }	
    
    if (label2 != NULL) {
	if (NULL == axis2)
	    axis2 = (struct Axis*) sf_alloc(1,sizeof(struct Axis));

	if (!sf_getfloat ("axisor2",&(axis2->or)))
	    axis2->or = (label2->where == 'l'? min1: max1);
	/* second axis origin */

	if (!sf_getint ("n2tic",&(axis2->ntic))) axis2->ntic = 1;
	/* second axis ticmarks */
	if (!sf_getfloat ("d2num", &(axis2->dnum)) ||
	    !sf_getfloat ("o2num", &(axis2->num0))) 
	    axis2->ntic = vp_optimal_scale(cube?
					   inch2*(mid2-min2)/
					   ((max2-min2)*aspect*labelsz):
					   inch2/(aspect*labelsz), 
					   label2->min, label2->max, 
					   &(axis2->num0), 
					   &(axis2->dnum));
    }

    if (label3 != NULL) {
	if (NULL == axis3)
	    axis3 = (struct Axis*) sf_alloc(1,sizeof(struct Axis));

	if (!sf_getfloat ("axisor3",&(axis3->or)))
	    axis3->or = min1;
	/* third axis origin */
	
	if (!sf_getint ("n3tic",&(axis3->ntic))) axis3->ntic = 1;
	/* third axis ticmarks */
	if (!sf_getfloat ("d3num", &(axis3->dnum)) ||
	    !sf_getfloat ("o3num", &(axis3->num0))) 
	    axis3->ntic = vp_optimal_scale(inch3/(aspect*labelsz), 
					   label3->min, label3->max, 
					   &(axis3->num0), 
					   &(axis3->dnum));
    }

    wheretics = (NULL != (where = sf_getstring ("wheretics"))) &&
	('a' == *where);
}


static void make_grid (bool grid)
{
    bool need;
    float num;

    if (axis1 != NULL) { 
	if (!sf_getbool("grid",&need) && !sf_getbool("grid1",&need))
	    need = transp? false: grid;
	/* to draw grid */

	if (need) {
	    grid1 = (struct Axis*) sf_alloc(1,sizeof(struct Axis));

	    grid1->num0 = axis1->num0;
	    grid1->or = axis1->or;

	    if (!sf_getfloat ("g1num",&(grid1->dnum))) {
		/* grid marks on first axis */
		grid1->dnum = axis1->dnum;
		grid1->ntic = axis1->ntic;
	    } else {
		grid1->ntic=0; 
		for (num=grid1->num0; num <= max1; num += grid1->dnum) {
		    grid1->ntic++;
		}
	    }
	}
    }

    if (axis2 != NULL) {
	if (!sf_getbool("grid",&need) && !sf_getbool("grid2",&need))
	    need = transp? grid: false;
	/* to draw grid */

	if (need) {
	    grid2 = (struct Axis*) sf_alloc(1,sizeof(struct Axis));
	    
	    grid2->num0 = axis2->num0;
	    grid2->or = axis2->or;
	    
	    if (!sf_getfloat ("g2num",&(grid2->dnum))) {
		/* grid marks on second axis */
		grid2->dnum = axis2->dnum;
		grid2->ntic = axis2->ntic;
	    } else {
		grid2->ntic=0; 
		for (num=grid2->num0; num <= max2; num += grid2->dnum) {
		    grid2->ntic++;
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
	grid3->or = axis3->or;
	grid3->dnum = axis3->dnum;
	grid3->ntic = axis3->ntic;
    }

    if (NULL != grid1 || NULL != grid2 || NULL != grid3) {
	if (!sf_getint("gridcol",&gridcol)) gridcol=grid? VP_RED: framecol;
	/* grid color */
	if (!sf_getint("gridfat",&gridfat)) gridfat=1;
	/* grid fatness */
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

    titlesz /= 33.;
    vs = titlesz * 0.6;

    if (title->where == 'l' || title->where == 'r') {	
	if (!sf_getbool ("labelrot",&labelrot)) labelrot = true;
	/* label rotation (for vertical labels */
	if (NULL != label2 && title->where == label2->where)
	    vs += 3.25*labelsz;

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
	    vs += 3.25*labelsz;
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
    
    if (!sf_getint ("barlabelfat",&(barlabel->fat))) barlabel->fat=0;
    /* bar label fatness */
    if (NULL == (labl=sf_getstring("barlabel")) && 
	NULL == (labl=sf_histstring(in,"label"))) { 
	/* bar label */
	barlabel->text = blank;
    } else if (NULL == (unit=sf_getstring("barunit")) && 
	       NULL == (unit=sf_histstring(in,"unit"))) {
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
	if (labelrot) vs *= 1.7;

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
    vp_set_dash (0);
}

void vp_simpleframe(void)
/*< Drawing simple frame >*/
{
    int i;
    float xc, yc, num;

    if (NULL != grid1 || NULL != grid2 || NULL != grid3) {
	vp_color (gridcol);
	vp_fat (gridfat);

	if (NULL != grid1) {
	    for (i=0; i < grid1->ntic; i++) {
		num = grid1->num0 + i*(grid1->dnum);
		if (fabsf(grid1->dnum) > FLT_EPSILON && 
		    fabsf(num) < FLT_EPSILON) num=0.;
		
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
		if (fabsf(grid2->dnum) > FLT_EPSILON && 
		    fabsf(num) < FLT_EPSILON) num=0.;

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
    }

    /* draw outline */   
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
}

void vp_framenum(float num)
/*< Outputting frame number >*/
{
    char string[80];
    float x, y;

    vp_umove(min1,min2);
    vp_where(&x,&y);

    sprintf (string, "%g", num);
    vp_tjust (TH_CENTER, TV_TOP);
    vp_gtext (x, y-4.*labelsz, labelsz, 0., 0., labelsz, string);
}

void vp_simplebarframe (void)
/*< Drawing simple frame for the bar label >*/
{
    float min, max;

/*    vp_clip(orig1+0.5*inch1,orig2-inch2,orig1+inch1,orig2+inch2); */

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
}

void vp_frame(void)
/*< Drawing frame >*/
{
    int i;
    float num, xc, yc, vs, xp[4], yp[4];
    char string[32];

    if (cube) {
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
	    xp[0] = min1;      yp[0] = mid2;
	    xp[1] = min1;      yp[1] = max2;
	    xp[2] = max1-mid1; yp[2] = max2;
	    vp_ufill(xp,yp,3);
	    xp[0] = mid1;      yp[0] = min2;
	    xp[1] = max1;      yp[1] = min2;
	    xp[2] = max1;      yp[2] = max2-mid2;
	    vp_ufill(xp,yp,3);
	}
    }

    vp_simpleframe();
    
    if (NULL != label1) {
	vp_fat (label1->fat);

	/* plot label */
	if (label1->where == 't') { 
	    vp_tjust (TH_CENTER, TV_BOTTOM);
	} else {
	    vp_tjust (TH_CENTER, TV_TOP);
	}
	vp_gtext(label1->x, label1->y, 
		 label1->xpath, label1->ypath, 
		 label1->xup, label1->yup, label1->text);

	/* plot tics */
	vs = label1->where == 't'? 0.5*labelsz: -0.5*labelsz;
	for (i=0; i < axis1->ntic; i++) {
	    num = axis1->num0 + i*(axis1->dnum);
	    if (fabsf(axis1->dnum) > FLT_EPSILON && 
		fabsf(num) < FLT_EPSILON) num=0.;
	    
	    if (cube) {
		xc = (num-label1->min)*(mid1-min1)/(label1->max-label1->min);
		yc = min2;
	    } else {
		xc = num;
		yc = axis1->or;
	    }

	    vp_umove (xc, yc);
	    vp_where (&xc, &yc);
	    vp_draw (xc, yc+vs);

	    snprintf (string,32,"%1.5g", num);
	    vp_gtext(xc, yc+1.5*vs, labelsz, 0., 0., labelsz, string);
	}
    }    

    if (NULL != label2) {
	vp_fat (label2->fat);

	/* plot label */
	if (label2->where == 'l') {
	    vp_tjust (TH_CENTER, TV_BOTTOM);
	} else {
	    vp_tjust (TH_CENTER, TV_TOP);
	}
	vp_gtext(label2->x, label2->y, 
		 label2->xpath, label2->ypath, 
		 label2->xup, label2->yup, label2->text);
	
        /* plot tics */
	vs = label2->where == 'l'? -0.5*labelsz: 0.5*labelsz;
	for (i=0; i < axis2->ntic; i++) {
	    num = axis2->num0 + i*(axis2->dnum);
	    if (fabsf(axis2->dnum) > FLT_EPSILON && 
		fabsf(num) < FLT_EPSILON) num=0.;

	    if (cube) {
		yc = mid2+(num-label2->max)*(mid2-min2)/
		    (label2->max-label2->min);
		xc = min1;
	    } else {
		yc = num;
		xc = axis2->or;
	    }	    

	    vp_umove (xc, yc);
	    vp_where (&xc, &yc);
	    vp_draw (xc+vs, yc);

	    snprintf (string,32,"%1.5g", num);

	    if (labelrot) {
		vp_gtext(xc+5.0*vs, yc, 0., -labelsz, labelsz, 0., string);
	    } else {
		vp_gtext(xc+1.5*vs, yc, 0., labelsz, -labelsz, 0., string);
	    }
	}
    }

    if (NULL != label3) {
	vp_fat (label3->fat);

	/* plot label */
	vp_tjust (TH_CENTER, TV_TOP);
	vp_gtext(label3->x, label3->y, 
		 label3->xpath, label3->ypath, 
		 label3->xup, label3->yup, label3->text);
	
        /* plot tics */
	vs = -0.5*labelsz;
	for (i=0; i < axis3->ntic; i++) {
	    num = axis3->num0 + i*(axis3->dnum);
	    if (fabsf(axis3->dnum) > FLT_EPSILON && 
		fabsf(num) < FLT_EPSILON) num=0.;

	    if (flat) {
		xc = mid1+(num-label3->min)*(max1-mid1)/
		    (label3->max-label3->min);
		yc = min2;
	    
		vp_umove (xc, yc);
		vp_where (&xc, &yc);
		vp_draw (xc, yc+vs);

		snprintf (string,32,"%1.5g", num);
		vp_gtext(xc, yc+1.5*vs,labelsz,0.,0.,labelsz,string);
	    } else {
		xc = mid1+(num-label3->min)*(max1-mid1)/
		    (label3->max-label3->min);
		yc = min2+(num-label3->min)*(max2-mid2-min2)/
		    (label3->max-label3->min);
	    
		vp_umove (xc, yc);
		vp_where (&xc, &yc);
		vp_draw (xc-vs*sinth, yc+vs*costh);
		
		snprintf (string,32,"%1.5g", num);
		vp_gtext(xc-1.5*vs*sinth, yc+1.5*vs*costh, 
			 labelsz*costh,labelsz*sinth,
			 -labelsz*sinth,labelsz*costh,string);
	    }
	}
    }
    
    if (NULL != title) {
	vp_fat (title->fat);
	if (title->where == 't' || 
	    (title->where == 'r' && labelrot) ||
	    (title->where == 'l' && !labelrot)) {
	    vp_tjust (TH_CENTER, TV_BOTTOM);
	} else {
	    vp_tjust (TH_CENTER, TV_TOP);
	}
	vp_gtext(title->x, title->y, 
		 title->xpath, title->ypath, 
		 title->xup, title->yup, title->text);
    }
    
    if (cube) {
	/* draw colored lines */
	vs = 0.5*labelsz;

	if (NULL != label2) {
	    yc = mid2-(mid2-min2)*(frame1+0.5)*d1/(label2->min-label2->max);

	    vp_color(cubelinecol);
	    vp_umove(min1,yc);
 
	    if (flat) {
		vp_udraw(max1,yc);
	    } else {
		vp_udraw(mid1,yc);
		vp_udraw(max1,yc+max2-mid2);
	    }

	    vp_color(VP_YELLOW);
	    vp_where(&xc,&yc);

	    num = label2->max+(frame1+0.5)*d1;
	    if (fabsf(d1) > FLT_EPSILON && 
		fabsf(num) < FLT_EPSILON) num=0.;
	    snprintf (string,32,"%1.5g",num);

	    vp_tjust (TH_CENTER, TV_TOP);
	    if (labelrot) {
		vp_gtext(xc+3.0*vs,yc,0.,-labelsz,labelsz,0.,string);
	    } else {
		vp_gtext(xc+0.5*vs, yc, 0., labelsz, -labelsz, 0., string);
	    }
	}

	if (NULL != label1) {
	    xc = min1+(mid1-min1)*(frame2+0.5)*d2/(label1->max-label1->min);

	    vp_color(cubelinecol);
	    vp_umove(xc,min2);
	    
	    if (flat) {
		vp_udraw(xc,max2);
	    } else {
		vp_udraw(xc,mid2);
		vp_udraw(xc+max1-mid1,max2);
	    }
	    
	    vp_color(VP_YELLOW);
	    vp_where(&xc,&yc);
	    
	    num = label1->min+(frame2+0.5)*d2;
	    if (fabsf(d2) > FLT_EPSILON && 
		fabsf(num) < FLT_EPSILON) num=0.;
	    snprintf (string,32,"%1.5g",num);
	    
	    vp_tjust (TH_CENTER, TV_BOTTOM);
	    vp_gtext(xc,yc+0.5*vs, labelsz, 0., 0., labelsz, string);
	}

	if (NULL != label3) {
	    yc = mid2+(max2-mid2)*(frame3+0.5)*d3/(label3->max-label3->min);
	    xc = mid1+(max1-mid1)*(frame3+0.5)*d3/(label3->max-label3->min);

	    vp_color(cubelinecol);

	    if (flat) {
		vp_umove(xc,min2);
		vp_udraw(xc,mid2);
	    } else {
		vp_umove(xc,yc+min2-mid2);
		vp_udraw(xc,yc); 
	    }

	    if (flat) {
		vp_umove(mid1,yc);
		vp_udraw(min1,yc);
	    } else {
		vp_umove(xc,yc);
		vp_udraw(xc+min1-mid1,yc);
	    }

	    vp_color(VP_YELLOW);
	    vp_where(&xc,&yc);

	    num = label3->min+(frame3+0.5)*d3;
	    if (fabsf(d3) > FLT_EPSILON && 
		fabsf(num) < FLT_EPSILON) num=0.;
	    snprintf (string,32,"%1.5g",num);

	    vp_tjust (TH_CENTER, TV_BOTTOM);

	    if (flat) {
		if (labelrot) {
		    vp_gtext(xc-4.0*vs, yc, 0.,-labelsz, labelsz,0.,string);
		} else {
		    vp_gtext(xc-0.5*vs, yc, 0., labelsz,-labelsz,0.,string);
		}
	    } else {
		vp_gtext(xc-0.5*vs*sinth, yc+0.5*vs*costh, 
			 labelsz*costh,labelsz*sinth,
			 -labelsz*sinth,labelsz*costh,string);
	    }
	} /* label3 */
    } /* if cube */

    vp_uclip (min1, min2, max1, max2);
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

	if (barlabel->where == 't' || 
	    barlabel->where == 'l') { 
	    vp_tjust (TH_CENTER, TV_BOTTOM);
	} else {
	    vp_tjust (TH_CENTER, TV_TOP);
	}

	vp_gtext(barlabel->x, barlabel->y, 
		 barlabel->xpath, barlabel->ypath, 
		 barlabel->xup, barlabel->yup, barlabel->text);

	/* plot tics */
	if (vertbar) {
	    vs = barlabel->where == 'l'? -0.5*labelsz: 0.5*labelsz;
	} else {
	    vs = barlabel->where == 't'? 0.5*labelsz: -0.5*labelsz;
	}

	for (i=0; i < baraxis->ntic; i++) {
	    num = baraxis->num0 + i*(baraxis->dnum);
	    if (fabsf(baraxis->dnum) > FLT_EPSILON && 
		fabsf(num) < FLT_EPSILON) num=0.;

	    snprintf (string,32,"%1.5g", num);	    
	    
	    if (vertbar) {
		yc = bar0 + i*dbar;
		xc = baraxis->or;
		
		vp_move (xc, yc);
		vp_draw (xc+vs, yc);

		if (labelrot) {
		    vp_gtext(xc+5.0*vs, yc, 0., 
			     -labelsz, labelsz, 0., string);
		} else {
		    vp_gtext(xc+1.5*vs, yc, 0.,
			     labelsz, -labelsz, 0., string);
		}
	    } else {
		xc = bar0 + i*dbar;
		yc = baraxis->or;

		vp_move (xc, yc);
		vp_draw (xc, yc+vs);
		
		vp_gtext(xc, yc+1.5*vs, labelsz, 0., 0., labelsz, string);
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
    cubelinecol = VP_YELLOW;
    vp_cubeframe(f1,f2,f3);
}

void vp_cubeframe(float f1, int f2, int f3   /* frame numbers */) 
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
    dc = cmax-cmin+FLT_EPSILON;

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
