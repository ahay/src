#include <rsf.h>

#include "stdplot.h"
#include "vplot.h"

static float min1, min2, max1, max2, labelsz;
static bool labelrot;
static int framecol;
static char blank[]=" ";

static struct Label {
    char where;
    int fat;
    float x, y, xpath, ypath, xup, yup;
    char* text;
} *label1=NULL, *label2=NULL, *title=NULL;

static void make_title (sf_file in);
	
void vp_stdplot_init (float umin1, float umax1, float umin2, float umax2)
{
    float scale1, scale2, orig1, orig2, uorig1, uorig2, crowd;
    float xll, xur, yll, yur, screenratio, screenht, screenwd;
    float inch1, inch2, marg;
    bool set;

    vp_style (STANDARD);

    /* get max and min */
    if (!sf_getfloat ("min1",&min1)) min1=umin1;
    if (!sf_getfloat ("min2",&min2)) min2=umin2;
    if (!sf_getfloat ("max1",&max1)) max1=umax1;
    if (!sf_getfloat ("max2",&max2)) max2=umax2;

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

static void make_title (sf_file in)
{
    bool want;
    char* where;
    float titlesz, vs, xc, yc;
    
    if (sf_getbool ("wanttitle",&want) && !want) {
	title = NULL;
	return;
    }

    title = (struct Label*) sf_alloc(1,sizeof(*title));

    if (NULL != (where = sf_getstring("wheretitle"))) {
	title->where = *where;
    } else {
	title->where = 'b';
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

void vp_frame_init (sf_file in, int col)
{
    framecol = col;
    make_title(in);
}

void vp_frame(void)
{
    vp_color(framecol);
    vp_umove(min1, min2);
    vp_udraw(min1, max2);
    vp_udraw(max1, max2);
    vp_udraw(max1, min2);
    vp_udraw(min1, min2);
   
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
