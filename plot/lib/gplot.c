#include <math.h>
#include <string.h>

#include <rsf.h>

#include "gplot.h"
#include "axis.h"
#include "vplot.h"

static int axisfat, axiscol, labelfat, *fat, *col, titlefat;
static float labelsz, min1, max1, min2, max2, inch1, inch2, titlesz; 
static float backcol[3],fillcol[3], xur, xll, yur, yll, *dashtype;
static float screenratio, screenht, screenwd;
static char *title, blank[]=" ";
static bool transp, yreverse, xreverse, where1, where2, labelrot, wheretics;
static bool fmin1, fmax1, fmin2, fmax2, wheretitle, verttitle, wanttitle;

struct Axis {
    int ntic;
    float or, dnum, num0, tic0, dtic;
    char *label;
    bool want;
} axis1, axis2;

static void getscl(struct Axis *axis);
static void labelaxis(struct Axis *axis);
static void tjust(bool where);
static void makelabel (float num, float num2, struct Axis *axis);
static void maketic (float num, struct Axis *axis);

/* initialize the length of axis and the number of tics */
void vp_axis_init (const sf_file in)
{
    float ftmp, num;
    bool want;
    char *stmp, *label, *unit;
    size_t len;
    struct Axis axis;

    if (!sf_getint ("axisfat",&axisfat)) axisfat=0;
    if (!sf_getint ("axiscol",&axiscol)) axiscol=7;
    if (!sf_getint ("labelfat",&labelfat)) labelfat=0;
    if (!sf_getfloat ("labelsz",&labelsz)) labelsz=8.;
    
    if ((NULL == (label=sf_getstring("label1"))) &&
	(NULL == (label=sf_histstring(in,"label1")))) {  
	axis1.label = blank;
    } else if ((NULL == (unit=sf_getstring("unit1"))) &&
	       (NULL == (unit=sf_histstring(in,"unit1")))) {
	axis1.label = label;
    } else {
	len = strlen(label)+strlen(unit)+4;
	axis1.label = sf_charalloc(len);
	snprintf(axis1.label,len,"%s (%s)",label,unit);
	free(label);
	free(unit);
    }

    if ((NULL == (label=sf_getstring("label2"))) &&
	(NULL == (label=sf_histstring(in,"label2")))) {
	axis2.label = blank;
    } else if ((NULL == (unit=sf_getstring("unit2"))) &&
	       (NULL == (unit=sf_histstring(in,"unit2")))) {
	axis2.label = label;           
    } else {
	len = strlen(label)+strlen(unit)+4;
	axis2.label = sf_charalloc(len);
	snprintf(axis2.label,len,"%s (%s)",label,unit);
	free(label);
	free(unit);
    }

    where1 = (NULL != (stmp = sf_getstring ("wherexlabel")) &&
	      'b' == *stmp);
    where2 = (NULL != (stmp = sf_getstring ("whereylabel")) &&
	      'r' == *stmp);
    
    /* checking to see if wantaxis is fetched */
    if (!sf_getbool ("wantaxis", &want)) {
	/* setting the default to put the axes on the plot */
        want = true;
	axis1.want = true;
	axis2.want = true;
    } else if (!want) {
	axis1.want = false;
	axis2.want = false;
    } else {
	if (!sf_getbool ("wantaxis1",&(axis1.want)))
	    axis1.want = true;
	if (!sf_getbool ("wantaxis2",&(axis2.want)))
	    axis2.want = true;
    }

    /* frame or axis */
    wheretics = NULL != (stmp = sf_getstring ("wheretics")) &&
	'a' == *stmp;

    if (!sf_getfloat ("axisor1",&(axis1.or)))
	axis1.or = where1? min2: max2;

    if (!sf_getfloat ("axisor2",&(axis2.or)))
	axis2.or = where2? min1: max1;

    if (!sf_getint ("n1tic",&(axis1.ntic))) axis1.ntic = 1;
    if (!sf_getint ("n2tic",&(axis2.ntic))) axis2.ntic = 1;

    if (!sf_getfloat ("d1num", &(axis1.dnum))) getscl (&axis1);
    if (!sf_getfloat ("d2num", &(axis2.dnum))) getscl (&axis2);

    if (0. == axis1.dnum) sf_error("%s: zero d1num",__FILE__);
    if (0. == axis2.dnum) sf_error("%s: zero d2num",__FILE__);

    if (!sf_getfloat("o1num",&(axis1.num0))) {
	ftmp = (min1 < max1)? min1: max1;
	for (num = floorf(ftmp / axis1.dnum) * axis1.dnum - axis1.dnum; 
	     num < ftmp; num += axis1.dnum) ;
	axis1.num0 = num;
    }

    if (!sf_getfloat("o2num",&(axis2.num0))) {
	ftmp = (min2 < max2)? min2: max2;
	for (num = floorf(ftmp / axis2.dnum) * axis2.dnum - axis2.dnum; 
	     num < ftmp; num += axis2.dnum) ;
	axis2.num0 = num;
    }

    axis1.dtic = axis1.dnum / axis1.ntic ;
    ftmp = (min1 < max1)? min1: max1;
    for (num = axis1.num0 - axis1.ntic * axis1.dtic;
	 num < ftmp; num += axis1.dtic);
    axis1.tic0 = num;

    axis2.dtic = axis2.dnum / axis2.ntic ;
    ftmp = (min2 < max2)? min2: max2;
    for (num = axis2.num0 - axis2.ntic * axis2.dtic;
	 num < ftmp; num += axis2.dtic);
    axis2.tic0 = num;

    if (transp) { /* swap */
	axis = axis1;
	axis1 = axis2;
	axis2 = axis;
    }
    
    if(yreverse) axis1.or = (min2+max2)-axis1.or;    
    if(xreverse) axis2.or = (min1+max1)-axis2.or;
}

void vp_color_init (void)
{
    if (!sf_getfloats ("backcol",backcol,3))
	backcol[0] = backcol[1] = backcol[2] = 0.;
    if (!sf_getfloats ("fillcol",fillcol,3))
	fillcol[0] = fillcol[1] = fillcol[2] = 0.;
    vp_coltab (0, backcol[0], backcol[1], backcol[2]);
}

void vp_coord_init (bool transp1, bool yreverse1)
{
    float crowd, marg;
    bool set;

    if (!sf_getfloat ("screenratio",&screenratio))
	screenratio = VP_SCREEN_RATIO;
    if (!sf_getfloat ("screenht",&screenht))
	screenht = VP_STANDARD_HEIGHT;
    if (!sf_getfloat ("screenwd",&screenwd))
	screenwd = screenht / screenratio;

    if (!sf_getfloat ("crowd",&crowd))   crowd  = 0.75;

    if (!sf_getfloat ("xinch",&inch1)) {
	if (!sf_getfloat ("crowd1",&inch1)) inch1 = crowd;
	inch1 *= screenwd;
    }
    if (!sf_getfloat ("yinch",&inch2)) {
	if (!sf_getfloat ("crowd2",&inch2)) inch2 = crowd;
	inch2 *= screenht;
    }

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
	    marg =  screenht - inch1;
	    yll = marg * 1./2.;
	    yur = screenht - marg * 1./2.;
	}
    } else if (!set) {
	yll = yur - inch2;
    }

    if (!sf_getbool ("transp",&transp)) transp = transp1;
    if (!sf_getbool ("xreverse",&xreverse)) xreverse = false;
    if (!sf_getbool ("yreverse",&yreverse)) yreverse = yreverse1;
    if (!sf_getbool ("labelrot",&labelrot)) labelrot = false;
}

/*
 *	0 continuous   DEFAULT
 *	1 fine dash
 *	2 fine dot
 *	3 dash
 *	4 large dash
 *	5 dot dash
 *	6 large dash small dash
 *	7 double dot
 *	8 double dash
 *	9 loose dash  The part after the decimal point determines 
 *                     the pattern repetition interval
 */
void vp_dash_fig (int type)
{
    float size;
    float dash[2], gap[2];
    
    switch (type) {
	case 1:
	    size = .3;
	    dash[0] = dash[1] = 1.;
	    gap[0] = gap[1] = 2.;
	    break;
	case 2:
	    size = .2;
	    dash[0] = dash[1] = 1.;
	    gap[0] = gap[1] = 6.;
	    break;
	case 3:
	    size = .4;
	    dash[0] = dash[1] = 4.;
	    gap[0] = gap[1] = 2.;
	    break;
	case 4:
	    size = .6;
	    dash[0] = dash[1] = 3.;
	    gap[0] = gap[1] = 2.;
	    break;
	case 5:
	    size = .5;
	    dash[0] = .3;
	    dash[1] = 3.;
	    gap[0] = gap[1] = 1.;
	    break;
	case 6:
	    size = .6;
	    dash[0] = 4.;
	    dash[1] = 2.;
	    gap[0] = gap[1] = 1.;
	    break;
	case 7:
	    size = .4;
	    dash[0] = dash[1] = 1.;
	    gap[0] = 2.;
	    gap[1] = 4.;
	    break;
	case 8:
	    size = .8;
	    dash[0] = dash[1] = 5.;
	    gap[0] = 2.;
	    gap[1] = 4.;
	    break;
	case 9:
	    size = .6;
	    dash[0] = dash[1] = 1.;
	    gap[0] = gap[1] = 1.;
	    break;
	case 0:
	    dash[0] = 0.;
	    gap[0] = 0.;
	    dash[1] = 0.;
	    gap[1] = 0.;
	    break;
	default:
	    dash[0] = 0.;
	    gap[0] = 0.;
	    dash[1] = 0.;
	    gap[1] = 0.;
	    break;
    }
    if (dash[0] + dash[1] + gap[0] + gap[1] != 0.) {
	/*If not default case then find the decimal part of dash->dashtype*/
	type = fmodf(type,1.);
	if (type == 0.) type = .4;
	size *= type / (.4 * (dash[0] + dash[1] + gap[0] + gap[1]));
	dash[0] *= size;
	dash[1] *= size;
	gap[0] *= size;
	gap[1] *= size;
    }
    vp_setdash (dash, gap, 2);
}

void vp_fillin (void)
{
    float xp[4], yp[4];
    vp_coltab (8, fillcol[0], fillcol[1], fillcol[2]);
    vp_color (8);
    xp[0] = min1;
    xp[1] = max1;
    xp[2] = max1;
    xp[3] = min1;
    yp[0] = min2;
    yp[1] = min2;
    yp[2] = max2;
    yp[3] = max2;
    vp_uarea (xp, yp, 4, 0, 1, 1);
} 

void vp_framenum (int i3, float d3, float o3, 
		  float xmin, float ymin, float labelsz)
{    
    float ch, xc, yc, f;
    char  string[80];

    vp_fat (0);
    f = i3 * d3 + o3;
    sprintf (string, (fmodf(f,1.) == 0.)? "%.0f": "%.3f", f);

    vp_umove (xmin, ymin);
    vp_where (&xc, &yc);
    ch = labelsz / 33.;
    tjust (true);
    vp_gtext (xc, yc - 4 * ch, ch, 0., 0., ch, string);
}

static void getscl(struct Axis *axis)
{
    float min, max, inch;

    if (axis == &axis2)  {
	min = min1;
	max = max1;
	inch = inch1;     
    } else {
	min = min2;
	max = max2;
	inch = inch2;
    }

    if (inch <= 0) sf_error("%s: negative inch %g",__FILE__,inch);

    /*******************************************************
    axis->dnum = vp_opttic(min,max,inch,axis->num0,labelsz);
    ********************************************************/
}

void vp_invmassage (float* min, float* max, float mid, float dev)
{
    float tmin, tmax;
    
    tmin = *min;
    tmax = *max;
    if (tmin == tmax) {
	*min = 0.;
	*max = 2*(tmax * dev + mid);
    }
}

static void labelaxis(struct Axis *axis)
{
    float vs, ch, xc, yc, pos1, pos2;
    float x, y, xup, yup, xpath, ypath;

    ch = labelsz / 33.;
    vs = 0.75 * ch;

    if (axis == &axis2) {
	pos1 = 0.5*(max1 + min1);
	pos2 = where1? min2: max2;

	x = xc;
	xpath = ch ;
	xup = 0.;
	ypath = 0.;
	yup = ch;

	vp_umove (pos1, pos2);
	vp_where (&xc, &yc);
	y = where1?
	    yc + ch + 2*vs:
	    yc - ch - 2*vs;
    } else {
	pos1 = where2? min1: max1;
	pos2 = 0.5*(min2 + max2);	
	
	y = yc;
	ypath = labelrot? ch: -ch;
	yup = 0.;
	xpath = 0.;
	xup = labelrot? -ch: ch;

	vp_umove (pos1, pos2);
	vp_where (&xc, &yc);

	if (labelrot) {
	    x = where2? 
		xc - ch - 2*vs: 
		xc + ch + 2*vs;
	} else {
	    x += where2?
		- (ch+vs):
		+ (ch+vs);
	}
    }

    vp_fat (labelfat);

    tjust((axis == &axis2)? where1: where2);
    vp_gtext (x, y, xpath, ypath, xup, yup, axis->label);
}

static void labeltic (struct Axis *axis)
{
    float ftmp, num, num0, dnum;
    float lmin1, lmax1, lmin2, lmax2, nmax;

    vp_fat (labelfat);

    if (axis == &axis2) {
	lmin1 = min2;
	lmax1 = max2;
	lmin2 = min1;
	lmax2 = max1;
    } else {
	lmax1 = max1;
	lmin1 = min1;
	lmax2 = max2;
	lmin2 = min2;
    }
    
    if (lmax1 < lmin1) {
	ftmp = lmin1;
	lmin1 = lmax1;
	lmax1 = ftmp;
    }
    
    num0 = axis->num0;    
    dnum = axis->dnum;
    
    if (num0 >= lmax1) num0 = lmin1 + (lmax1 - lmin1) / 2.;
    
    if ((xreverse && axis == &axis1) || 
	(yreverse && axis == &axis2)) {
        nmax = num0;
        num0 = lmin1 + lmax1 - num0;
	for (num = num0; num >= lmin1; num -= dnum) {
	    makelabel (num, nmax, axis); 
            nmax += dnum ;  
	}
    } else {
	for (num = num0; num <= lmax1; num += dnum) {
	    makelabel (num,  num, axis);
	}
    }
}

static void makelabel (float num, float num2, struct Axis *axis)
{
    float xc, yc, lmin1, lmin2, lmax1, lmax2, point1, point2;
    float x, y, xup, yup, xpath, ypath, ftmp, num1, ch, vs;
    const float eps=.0000001;
    char string[80];

    ch = labelsz / 33.;
    vs = 0.75 * ch;
    
    if (axis == &axis2) {
	lmin1 = min2;
	lmax1 = max2;
	lmin2 = min1;
	lmax2 = max1;
    } else {
	lmax1 = max1;
	lmin1 = min1;
	lmax2 = max2;
	lmin2 = min2;
    }


    if (lmax1 < lmin1 || 
	(xreverse  && axis == &axis1) || 
	(yreverse  && axis == &axis2)) {
	ftmp = lmin1;
	lmin1 = lmax1;
	lmax1 = ftmp;
    }

    if (fabsf (num) < (lmax1 - lmin1) / 10000) num = 0.0;

    if (!wheretics) {
	if (axis == &axis2) {
	    point1 = num;
	    point2 = where1? min2: max2;
	} else {
	    point1 = where2? max2: min2;
	    point2 = num;
	}
    } else {
	if (axis == &axis2) {
	    point1 = num;
	    point2 = axis->or;
	} else {
	    point1 = axis->or;
	    point2 = num;
	}
    }

    vp_umove (point1, point2);
    vp_where (&xc, &yc);

    num1 = (fabsf(num2) < eps)? 0.0: num2;
    sprintf (string, "%1.5g", num1);

    if (axis == &axis1) {
	y = yc;
	xpath = 0.;
	ypath =  labelrot * ch;
	xup =  -labelrot * ch;
	yup = 0;
	if (!where2) {
	    x = xc - vs;
	    if (!labelrot) x -= (ch + vs);
	} else {
	    x = xc + vs;
	    if (!labelrot) x += (ch + vs);
        }
    } else {
	x = xc;
	yup = ch;
	xpath = ch;
	ypath = 0.;
	xup = 0.;
	y = where1? yc - vs: yc + vs;
    }

    tjust ((axis == &axis2)? where1: where2);
    vp_gtext (x, y, xpath, ypath, xup, yup, string);
}

void vp_massage (float *min, float *max, float *mid, float *dev)
{
    float lmin, lmax, lmid, ldev;

    lmin = *min;
    lmax = *max;
    lmid = (lmin + lmax)/2.0;
    ldev=fabsf(lmin-lmid);

    if (0. != ldev) {
	*min = (lmin - lmid) / ldev;
	*max = (lmax - lmid) / ldev;
    } else {
	*min = lmin - lmid;
	*max = lmax - lmid;
    }
    
    *mid = lmid;
    *dev = ldev;
}

void vp_minmax (float emin1, float emin2, float emax1, float emax2)
{
    bool btmp;
    float ftmp;

    fmin1 = (!sf_getfloat ("min1",&min1)); if(fmin1) min1=emin1;
    fmin2 = (!sf_getfloat ("min2",&min2)); if(fmin2) min2=emin2;
    fmax1 = (!sf_getfloat ("max1",&max1)); if(fmax1) max1=emax1;
    fmax2 = (!sf_getfloat ("max2",&max2)); if(fmax2) max2=emax2;

    if (transp) {
	ftmp=min1; min1=min2; min2=ftmp;
	btmp=fmin1; fmin1=fmin2; fmin2=btmp;

	ftmp=max1; max1=max2; max2=ftmp;
	btmp=fmax1; fmax1=fmax2; fmax2=btmp;
    }
}

void vp_pad_init(bool pad, bool npad)
{
    float mid1, mid2, dev1, dev2;

    if (pad) {
        vp_massage(&min1, &max1, &mid1, &dev1);
        if (!fmin1 || npad) min1 *= 1.04;
	if (!fmax1 || npad) max1 *= 1.04;
        vp_invmassage(&min1, &max1, mid1, dev1);

	vp_massage(&min2, &max2, &mid2, &dev2);
        if (!fmin2 || npad) min2 *= 1.04;
	if (!fmax2 || npad) max2 *= 1.04;
        vp_invmassage(&min2, &max2, mid2, dev2);
    }
    
    if (min1 == 0 && max1 == 0) {
	min1 = -1.;
	max1 = 1.;
    }
    
    if (min2 == 0 && max2 == 0) {
	min2 = -1.;
	max2 = 1.;
    }

    if (min1 == max1) {
        max1 *= 1.04;
        min1 *= 0.96;
    }

    if (min2 == max2) {
        max2 *= 1.04;
        min2 *= 0.96;
    }
}

static void plotaxis(struct Axis *axis)
{
    vp_fat(axisfat);
    if (axis == &axis1) {
	vp_umove(min1, axis->or);
	vp_udraw(max1, axis->or);
    } else {
	vp_umove(axis->or, min2);
	vp_udraw(axis->or, max2);
    }
} 

void vp_plotframe(int col)
{
    vp_color(col);
    vp_umove(min1, min2);
    vp_udraw(min1, max2);
    vp_udraw(max1, max2);
    vp_udraw(max1, min2);
    vp_udraw(min1, min2);
}

void vp_plot_init(int n2)
{
    int  i, len;

    dashtype = sf_floatalloc(n2);
    if (!sf_getfloats ("dash",dashtype,n2)) {
	/*  line dash type	
	    0 continuos (default)
	    1 fine dash
	    2 fine dot
	    3 dash
	    4 large dash
	    5 dot dash
	    6 large dash small dash
	    7 double dot
	    8 double dash
	    9 loose dash  The part after the decimal point determines the pattern repetition interval */	   
	for (i = 0; i < n2; i++) 
	    dashtype[i] = 0.;
    }

    fat = sf_intalloc(n2);
    if (!sf_getints("plotfat",fat,n2)) {
	/* line fatness */
	for (i = 0; i < n2; i++)
	    fat[i] = 0;
    }
    
    col = sf_intalloc(n2);
    if (!sf_getints("plotcol",col,n2)) {
	/* line color 
	   7 white
	   6 yellow (default)
	   5 cyan
	   4 green
	   3 magenta
	   2 red
	   1 blue
	   0 black */
	for (i = 0; i < n2; i++)
	    col[i] = 6 - (i % 6);
    }
}

void vp_plot_param(void)
{
    if (fillcol[0] != backcol[0] || 
	fillcol[1] != backcol[1] || 
	fillcol[2] != backcol[2]) vp_fillin();
    vp_coltab(0, backcol[0], backcol[1], backcol[2]);
}

static void plottic (struct Axis * axis)
{
    float ftmp, lmin1, lmin2, lmax1, lmax2, num, num0, dnum;

    num0 = axis->num0;
    dnum = axis->dnum;
    vp_fat (axisfat);

    if (axis == &axis2) {
	lmin1 = min2;
	lmax1 = max2;
	lmin2 = min1;
	lmax2 = max1;
    } else {
	lmin1 = min1;
	lmax1 = max1;
	lmin2 = min2;
	lmax2 = max2;
    }

    if (lmax1 < lmin1){
	ftmp = lmin1;
	lmin1 = lmax1;
	lmax1 = ftmp;
    } 
    
    if (num0 >= lmax1) num0 = lmin1 + (lmax1-lmin1)/2.; 

    if ((xreverse && axis == &axis1) || 
	(yreverse && axis == &axis2)) {
	ftmp = lmin1;
	lmin1 = lmax1;
	lmax1 = ftmp;
        num0 = min1 + max1 - num0;
	dnum = - dnum;
    }

    for (num = num0; num >= max1; num += dnum) {
	maketic (num, axis);
    }
}

static void maketic (float num, struct Axis *axis)
{
    float xc, yc, pos1, pos2, ch;
    float lmin1, lmax1, lmin2, lmax2, tic, num0, dtic, dnum, vs;

    ch = labelsz / 33.;
    vs = 0.5 * ch;

    if (axis == &axis2) {
	lmin1 = min2;
	lmin2 = min1;
	lmax1 = max2;
	lmax2 = max1;
    } else {
	lmin1 = min1;
	lmin2 = min2;
	lmax1 = max1;
	lmax2 = max2;
    }

    if (fabsf (num) < (max1 - min1)/10000.) num = 0.;

    if (axis == &axis2) {
	pos1 = num;
	if (wheretics) {
	    pos2 = axis->or;
	} else {
	    pos2 = where2? lmin2: lmax2;
	}
    } else {
	if (wheretics) {
	    pos1 = axis->or;
	} else {
	    pos1 = where1? lmin2: lmax2;
	}
	pos2 = num;
    }
    
    vp_umove (pos1, pos2);
    vp_where (&xc, &yc);
    
    if (axis == &axis2) {
	pos1 = xc;
	pos2 = where2? yc-vs: yc+vs;
    } else {
	pos1 = where1? xc-vs: xc+vs;
	pos2 = yc;
    }
    
    vp_draw (pos1, pos2);

    num0 = axis->num0;
    dnum = axis->dnum;
    dtic = axis->dtic;

    if (dtic != dnum) {
	for (tic = num0; tic <= lmax1; tic += dtic) {
	    if (axis == &axis2) {
		pos1 = tic;
		pos2 = where2? lmin2: lmax2;
	    } else {
		pos1 = where1? lmin2: lmax2;
		pos2 = tic;
	    }

	    vp_umove (pos1, pos2);
	    vp_where (&xc, &yc);

	    if (axis == &axis2) {
		pos1 = xc;
		if (axis->or == lmin2) {
		    pos2 = where2? yc-0.5*vs: yc+0.5*vs;
		} else {
		    pos2 = where2? yc-0.25*vs: yc+0.25*vs;
		}
	    } else {
		if (axis->or == lmin2) {
		    pos1 = where1? xc-0.5*vs: xc+0.5*vs;
		} else {
		    pos1 = where1? xc-0.25*vs: xc+0.25*vs;
		}
		pos2 = yc;
	    }

	    vp_draw (pos1, pos2);
	}
    }
}

void vp_plottitle(void)
{
    float labelvs, labelch, ch, vs, xc, yc, pos1, pos2;
    float x, y, xup, yup, xpath, ypath;

    ch = titlesz / 33.;
    vs = ch * 0.6;

    vp_fat (titlefat);
    vp_color(axiscol);
    
    if (( verttitle && wheretitle != where2) ||
	(!verttitle && wheretitle != where1)) {
	labelch = labelsz / 33.;
	labelvs = 0.75 * labelch;
    } else {
	labelch = 0.;
	labelvs = 0.; 
    }
     
    if (!verttitle) {
	pos1 = 0.5*(max1+min1);
	pos2 = wheretitle? max2: min2;

	ypath = 0.;
	yup = ch;
	xpath = ch;
	xup = 0.;

	vp_umove (pos1, pos2);
	vp_where (&xc, &yc);
	x = xc;
	y = wheretitle? 
	    yc + vs + labelch + 3*labelvs:
	    yc - vs - labelch - 3*labelvs; 
    } else {
	pos2 = 0.5*(max2+min2);
	pos1 = wheretitle? max1: min1;

	ypath = labelrot? 0.: -ch;
	yup = 0.;
	xpath = 0.;
	xup = labelrot? 0.: ch;

	vp_umove (pos1, pos2);
	vp_where (&xc, &yc);
	y = yc;
	if (labelrot) {
	    x = wheretitle?
		xc - vs - labelch - 3*labelvs:
		xc + vs + labelch + 3*labelvs;
	} else {
	    x += wheretitle?
		-(ch + vs + labelvs):
		+(ch + vs + labelvs);
	}
    }	   
 
    tjust(wheretitle);
    vp_gtext(x, y, xpath, ypath, xup, yup ,title);
}

void vp_stdplot (int fastplt, bool wantframe, bool wantframenum,
		 int n3, float o3, float d3)
{
    if (fastplt < 20) {
	vp_clip (-VP_MAX, -VP_MAX, VP_MAX, VP_MAX);
	vp_color (axiscol);
	vp_fat (axisfat);

	if (wantframe && fastplt < 16) vp_plotframe (axiscol);

	if (fastplt < 14) {
	    if (axis1.want) {
		if (!wantframe || axis1.or != min1) plotaxis (&axis1);
		if (fastplt < 12) {
		    if (axis1.dnum != 0.) {
			if (axis1.ntic != 0) plottic (&axis1);
			labeltic (&axis1);
		    }
		    labelaxis (&axis1);
		}
	    }
	    if (axis2.want) {
		if (!wantframe || axis2.or != min2) plotaxis (&axis2);
		if (fastplt < 12) {
		    if (axis2.dnum != 0.) {
			if (axis2.ntic != 0) plottic (&axis2);
			labeltic (&axis2);
		    }
		    labelaxis (&axis2);
		}
	    }
	}

	if (fastplt < 3 && wanttitle) vp_plottitle ();

	if (fastplt < 13 && wantframenum && n3 > 1) {
	    vp_framenum (n3, d3, o3, min1, min2, labelsz);
	}
    }
}

void vp_title_init(sf_file file)
{
    char *where;

    if (NULL != (where = sf_getstring("wheretitle"))) {
	verttitle = (*where == 'l' || *where == 'r');
	wheretitle = (*where == 'l' || *where == 't');
    } else {
	verttitle = false;
	wheretitle = true;
    }
    if (!sf_getbool ("wanttitle",&wanttitle)) wanttitle=true;
    if (!sf_getint ("titlefat",&titlefat)) titlefat=0;
    if (!sf_getfloat ("titlesz",&titlesz)) titlesz=10.;
    if (NULL == (title = sf_getstring("title")) &&
	NULL == (title = sf_histstring(file,"title")) &&
	NULL == (title = sf_histstring(file,"in"))) 
	title = blank;
}

static void tjust(bool where)
{
    if (where) {
	vp_tjust (TH_CENTER, TV_TOP);
    } else {
	vp_tjust (TH_CENTER, TV_BOTTOM);
    }

/*	vp_tjust (TH_SYMBOL, TV_SYMBOL); */
}

void vp_vplot_init (void)
{
    float scale1, scale2, orig1, orig2, uorig1, uorig2;

    vp_style (STANDARD);
    if (min2 == max2) sf_error ("min2 is equal to max2");
    if (min1 == max1) sf_error ("min1 is equal to max1");

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

    scale1 = inch1 / (max1 - min1);
    scale2 = inch2 / (max2 - min2);
    uorig2 = (min2 + max2)*0.5;
    uorig1 = (min1 + max1)*0.5;
    
    vp_scale (scale1, scale2);
    vp_orig (orig1, orig2);
    vp_uorig (uorig1, uorig2);

    vp_uclip (min1, min2, max1, max2);
}

void vp_rotate (int n, float* x, float* y)
{
    int i;

    if (yreverse) {
	for (i=0; i < n; i++) {
	    y[i] = (min2 + max2) - y[i];
	}
    } 

    if (xreverse) {
	for (i=0; i < n; i++) {
	    x[i] = (min1 + max1) - x[i];
	}
    } 
}
