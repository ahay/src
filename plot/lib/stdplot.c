#include <float.h>
#include <math.h>
#include <stdio.h>

#include <rsf.h>

#include "stdplot.h"
#include "axis.h"
#include "vplot.h"

static float min1, min2, max1, max2, inch1, inch2, orig1, orig2;
static float labelsz, barlabelsz, barmin, barmax;
static int framecol;
static bool labelrot, transp, wheretics, scalebar, vertbar;
static char blank[]=" ";

static struct Label {
    char where;
    int fat;
    float x, y, xpath, ypath, xup, yup;
    char* text;
} *label1=NULL, *label2=NULL, *title=NULL, *barlabel=NULL;

static struct Axis {
    int ntic;
    float or, dnum, num0;
} *axis1=NULL, *axis2=NULL;

static void make_title (sf_file in, char wheret);
static void make_barlabel (void);
static void make_labels (sf_file in, char where1, char where2);
static void make_axes (void);
static void swap(float*a,float*b);

static void swap(float*a,float*b)
{
    float f;

    f = *a; *a = *b; *b=f;
}

void vp_stdplot_init (float umin1, float umax1, float umin2, float umax2,
		      bool transp1, bool xreverse1, bool yreverse1, bool pad1)
{
    bool pad, set, xreverse, yreverse;
    float mid, off, scale1, scale2, uorig1, uorig2, crowd, barwd;
    float xll, xur, yll, yur, screenratio, screenht, screenwd, marg;
    char* bartype;

    transp = transp1;
    if (!sf_getbool ("xreverse",&xreverse)) xreverse = xreverse1;
    if (!sf_getbool ("yreverse",&yreverse)) yreverse = yreverse1;
    if (!sf_getbool ("pad",&pad)) pad = pad1;
    if (!sf_getbool ("wantscalebar",&scalebar)) scalebar = false;

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
    if (!sf_getfloat ("min2",&min2)) min2=umin2;
    if (!sf_getfloat ("max1",&max1)) max1=umax1;
    if (!sf_getfloat ("max2",&max2)) max2=umax2;

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
    if (!sf_getfloat ("screenht",&screenht))
	screenht = VP_STANDARD_HEIGHT;
    if (!sf_getfloat ("screenwd",&screenwd))
	screenwd = screenht / screenratio;

    /* get inches */
    if (!sf_getfloat ("crowd",&crowd)) crowd  = 0.75;
    
    if (!sf_getfloat ("xinch",&inch1)) {
	if (!sf_getfloat ("crowd1",&inch1)) inch1 = crowd;
	inch1 *= screenwd;
    }
    if (!sf_getfloat ("yinch",&inch2)) {
	if (!sf_getfloat ("crowd2",&inch2)) inch2 = crowd;
	inch2 *= screenht;
    }

    /* get corners */
    set = sf_getfloat ("xll",&xll);
    if (!sf_getfloat ("xur",&xur)) {
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
    if (!sf_getfloat ("yur",&yur)) {
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
	vertbar = (NULL == (bartype = sf_getstring("bartype"))) ||
	    (bartype[0] == 'v');
	if (!sf_getfloat("barwidth",&barwd)) barwd = 0.36;
	if (vertbar) {
	    xur -= (0.07*screenwd + barwd);
	    xll -= 0.5*barwd;
	    barmax = 0.88*screenwd;
	    barmin = barmax-barwd;
	} else {
	    yur += 0.4*barwd;
	    yll += (0.04*screenht + barwd);
	    barmin = 0.12*screenht;
	    barmax = barmin+barwd;
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
    
    /* find scales and user origin */
    scale1 = inch1 / (max1 - min1);
    scale2 = inch2 / (max2 - min2);
    uorig2 = (min2 + max2)*0.5;
    uorig1 = (min1 + max1)*0.5;
    
    /* setup the coordinate system */
    vp_scale (scale1, scale2);
    vp_orig (orig1, orig2);
    vp_uorig (uorig1, uorig2);

    if (!sf_getint ("axiscol",&framecol)) framecol=7;
}

static void make_labels (sf_file in, char where1, char where2)
{
    float vs, xc, yc;
    bool want;
    char* where;
    struct Label *label;

    if (sf_getbool ("wantaxis", &want) && !want) {
	label1 = NULL;
	label2 = NULL;
	return;
    }

    if (sf_getbool ("wantaxis1",&want) && !want) {
	label1 = NULL;
    } else if (NULL == label1) { 
	label1 = (struct Label*) sf_alloc(1,sizeof(struct Label));
    }

    if (sf_getbool ("wantaxis2",&want) && !want) {
	label2 = NULL;
	if (NULL == label1) return;
    } else if (NULL == label2) {
	label2 = (struct Label*) sf_alloc(2,sizeof(struct Label));
    }

    if (transp) {
	label = label1;
	label1 = label2;
	label2 = label;
    }

    if (!sf_getfloat ("labelsz",&labelsz)) labelsz=8.;
    labelsz /= 33.;
    vs = 2.5*labelsz;

    if (NULL != label1) {
	if (NULL != (where = sf_getstring("wherexlabel"))) {
	    label1->where = *where;
	} else {
	    label1->where = where1;
	}
	if (!sf_getint ("labelfat",&(label1->fat))) label1->fat=0;
	if ((NULL == (label1->text=sf_getstring(transp? "label2":"label1"))) &&
	    (NULL == (label1->text=sf_histstring(in,
						 transp? "label2":"label1"))))
	    label1->text = blank;
	
	label1->xpath = labelsz;
	label1->xup = 0.;
	label1->ypath = 0.;
	label1->yup = labelsz;

	xc = 0.5*(max1 + min1);
	yc = (label1->where == 't') ? max2: min2;
	vp_umove (xc, yc);
	vp_where (&xc, &yc);

	label1->x = xc;
	label1->y = (label1->where == 't') ? yc+vs: yc-vs;
    }

    if (NULL != label2) {
	if (NULL != (where = sf_getstring("whereylabel"))) {
	    label2->where = *where;
	} else {
	    label2->where = where2;
	}
	if (!sf_getint ("labelfat",&(label2->fat))) label2->fat=0;
	if ((NULL == (label2->text=sf_getstring(transp? "label1":"label2"))) &&
	    (NULL == (label2->text=sf_histstring(in,
						 transp? "label1":"label2"))))
	    label2->text = blank;

	if (!sf_getbool ("labelrot",&labelrot)) labelrot = true;
	if (labelrot) vs *= 1.7;

	label2->ypath = labelrot? -labelsz: labelsz;
	label2->yup = 0.;
	label2->xpath = 0.;
	label2->xup = labelrot? labelsz: -labelsz;

	xc  = (label2->where == 'l')? min1: max1;
	yc = 0.5*(min2 + max2);

	vp_umove (xc, yc);
	vp_where (&xc, &yc);

	label2->y = yc;	
	label2->x = (label2->where == 'l')? xc-vs: xc+vs;
    }
}

static void make_axes (void)
{
    const float aspect=0.8;
    char* where;

    if (label1 != NULL) {
	if (NULL == axis1) 
	    axis1 = (struct Axis*) sf_alloc(1,sizeof(struct Axis));

	if (!sf_getfloat ("axisor1",&(axis1->or)))
	    axis1->or = (label1->where == 'b'? min2: max2);

	if (!sf_getint ("n1tic",&(axis1->ntic))) axis1->ntic = 1;
	if (!sf_getfloat ("d1num", &(axis1->dnum)) ||
	    !sf_getfloat ("o1num", &(axis1->num0))) 
	    axis1->ntic = vp_optimal_scale(inch1/(aspect*labelsz), 
					   min1, max1, 
					   &(axis1->num0), 
					   &(axis1->dnum));
    }	
    
    if (label2 != NULL) {
	if (NULL == axis2)
	    axis2 = (struct Axis*) sf_alloc(1,sizeof(struct Axis));

	if (!sf_getfloat ("axisor2",&(axis2->or)))
	    axis2->or = (label2->where == 'l'? min1: max1);
	
	if (!sf_getint ("n2tic",&(axis2->ntic))) axis2->ntic = 1;
	if (!sf_getfloat ("d2num", &(axis2->dnum)) ||
	    !sf_getfloat ("o2num", &(axis2->num0))) 
	    axis2->ntic = vp_optimal_scale(inch2/(aspect*labelsz), 
					   min2, max2, 
					   &(axis2->num0), 
					   &(axis2->dnum));
    }

    wheretics = (NULL != (where = sf_getstring ("wheretics"))) &&
	('a' == *where);
}

static void make_title (sf_file in, char wheret)
{
    bool want;
    char* where;
    float titlesz, vs, xc, yc;
    
    if (sf_getbool ("wanttitle",&want) && !want) {
	title = NULL;
	return;
    }

    if (NULL == title)
	title = (struct Label*) sf_alloc(1,sizeof(struct Label));

    if (NULL != (where = sf_getstring("wheretitle"))) {
	title->where = *where;
    } else {
	title->where = wheret;
    }

    if (!sf_getint ("titlefat",&(title->fat))) title->fat=0;

    if (!sf_getfloat ("titlesz",&titlesz)) titlesz=10.;

    if (NULL == (title->text = sf_getstring("title")) &&
	NULL == (title->text = sf_histstring(in,"title")) &&
	NULL == (title->text = sf_histstring(in,"in"))) 
	title->text = blank;

    titlesz /= 33.;
    vs = titlesz * 0.6;

    if (title->where == 'l' || title->where == 'r') {	
	if (!sf_getbool ("labelrot",&labelrot)) labelrot = true;
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

static void make_barlabel (void)
{
    float vs, xc, yc;
    bool want;
    char* where;

    if (sf_getbool ("wantbaraxis", &want) && !want) {
	barlabel = NULL;
	return;
    }

    barlabel = (struct Label*) sf_alloc(1,sizeof(struct Label));

    if (!sf_getfloat ("barlabelsz",&barlabelsz)) {
	barlabelsz=labelsz;
    } else {
	barlabelsz /= 33.;
    }
    vs = 2.5*barlabelsz;
    if (vertbar && labelrot) vs *= 1.7;

    if (!sf_getint ("barlabelfat",&(barlabel->fat))) barlabel->fat=0;
    if (NULL == (barlabel->text=sf_getstring("barlabel"))) 
	barlabel->text = blank;

    if (NULL != (where = sf_getstring("wherebarlabel"))) {
	barlabel->where = *where;
    } else {
	barlabel->where = vertbar? 'r':'b';
    }

    if (vertbar) {
	barlabel->ypath = labelrot? -labelsz: labelsz;
	barlabel->yup = 0.;
	barlabel->xpath = 0.;
	barlabel->xup = labelrot? labelsz: -labelsz;

	xc  = (barlabel->where == 'l')? min1: max1;
	yc = 0.5*(min2 + max2);

	vp_umove (xc, yc);
	vp_where (&xc, &yc);

	barlabel->y = yc;	
	barlabel->x = (barlabel->where == 'l')? xc-vs: xc+vs;
    } else {
	barlabel->xpath = labelsz;
	barlabel->xup = 0.;
	barlabel->ypath = 0.;
	barlabel->yup = labelsz;

	xc = 0.5*(max1 + min1);
	yc = (barlabel->where == 't') ? max2: min2;
	vp_umove (xc, yc);
	vp_where (&xc, &yc);

	barlabel->x = xc;
	barlabel->y = (barlabel->where == 't')? yc+vs: yc-vs;
    }
}

void vp_frame_init (sf_file in, const char* where)
{
    make_labels(in,where[0],where[1]);
    if (scalebar) make_barlabel();
    make_axes();
    make_title(in,where[2]);
}

void vp_simpleframe(void)
{
    vp_color(framecol);
    vp_umove(min1, min2);
    vp_udraw(min1, max2);
    vp_udraw(max1, max2);
    vp_udraw(max1, min2);
    vp_udraw(min1, min2);
}

void vp_barframe (void)
{
    float min, max;

/*    vp_clip(orig1+0.5*inch1,orig2-inch2,orig1+inch1,orig2+inch2); */

    vp_color(framecol);
    if (vertbar) {
	min = orig2-0.5*inch2;
	max = orig2+0.5*inch2;

	vp_clip(barmin-0.1*inch1,min-0.1*inch2,
		barmax+0.1*inch1,max+0.1*inch2);

	vp_move(barmin,min);
	vp_draw(barmax,min);
	vp_draw(barmax,max);
	vp_draw(barmin,max);
	vp_draw(barmin,min);
    } else {
	min = orig1-0.5*inch1;
	max = orig1+0.5*inch1;

	vp_clip(min-0.1*inch1,barmin-0.1*inch2,
		max+0.1*inch1,barmax+0.1*inch2);

	vp_move(min,barmin);
	vp_draw(min,barmax);
	vp_draw(max,barmax);
	vp_draw(max,barmin);
	vp_draw(min,barmin);
    }
}

void vp_frame(void)
{
    int i;
    float num, xc, yc, vs;
    char string[32];

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
	    
	    xc = num;
	    yc = axis1->or;

	    vp_umove (xc, yc);
	    vp_where (&xc, &yc);
	    vp_draw (xc, yc+vs);

	    sprintf (string, "%1.5g", num);
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

	    yc = num;
	    xc = axis2->or;
	    
	    vp_umove (xc, yc);
	    vp_where (&xc, &yc);
	    vp_draw (xc+vs, yc);

	    sprintf (string, "%1.5g", num);

	    if (labelrot) {
		vp_gtext(xc+5.0*vs, yc, 0., -labelsz, labelsz, 0., string);
	    } else {
		vp_gtext(xc+1.5*vs, yc, 0., labelsz, -labelsz, 0., string);
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
    
    vp_uclip (min1, min2, max1, max2);
}
