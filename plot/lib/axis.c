#include <math.h>
#include <stdio.h>
#include <string.h>
#include <float.h>

#include "axis.h"
#include "vplot.h"

static float power (float f, int ie);
static float nice_number (float d);
static int optimal_scale(int n, float min, float max, 
			 float *onum, float *dnum);

/* axis goes from (x1,y1) to (x2,y2) in vplot coordinates
 * its scale ranges from num1 to num2, ltic is the size of tick marks 
 * if dnum=0., an optimal linear scale is estimated */
void vp_simple_axis (float x1, float y1, 
		     float x2, float y2, 
		     float num1, float num2,
		     float onum, float dnum, 
		     float ltic, char* label, float size)
{
    int nopt, i;
    float xpath, ypath, xup, yup, dist, costh, sinth, dtic, otic;
    float dxtic, dytic, loc, num, xpos, ypos;
    char string[32];
    const float pad = 0.15; /* between tic and number, number and label */
    const float aspect = 0.8;

    vp_color (7); /* white */

    /* draw axis */
    vp_move(x1,y1);
    vp_draw(x2,y2);

    /* compute sines and cosines of axis angle */
    dist = hypotf(x2-x1,y2-y1);
    costh = (x2-x1)/dist;
    sinth = (y2-y1)/dist;

    /* x and y sizes of tic marks */
    dxtic = ltic * sinth;
    dytic = -ltic * costh;

    if (x1 <= x2) {
        vp_tjust(TH_CENTER,TV_TOP);
        xpath = size * costh;
        ypath = size * sinth;
        xup = -size * sinth;
        yup = size * costh;
    } else {	
        vp_tjust(TH_CENTER,TV_BOTTOM);
        xpath = -size * costh;
        ypath = -size * sinth;
        xup = size * sinth;
        yup = -size * costh;
    } 

    if (dnum == 0.) { /* figure out onum and dnum */
	nopt = vp_optimal_scale(dist/(aspect*size), num1, num2, &onum, &dnum);
    } else { /* use client's scale */
	nopt = 1+ (int) (floor) (num2-onum)/dnum;
	if (onum+(nopt-1)*dnum > num2) nopt--;
    }

    /* figure out the tic mark spacing */
    otic = (onum-num1)*dist/(num2-num1);
    dtic = dnum*dist/(num2-num1); 

    /* move to each tic mark location, draw the tic and the number */
    for (i=0; i < nopt; i++) {
	num = onum + i*dnum;
        sprintf(string,"%1.5g",num);
	loc = otic + i*dtic;
	xpos = x1 + loc * costh;
	ypos = y1 + loc * sinth;
        vp_move(xpos,ypos);
	vp_draw(xpos+dxtic,ypos+dytic);
	vp_gtext(xpos+dxtic+pad*sinth,
		 ypos+dytic-pad*costh,xpath,ypath,xup,yup,string);
    }
    
    /* now the axis label */
    xpos = x1 + dist/2. * costh + dxtic + size * sinth + 2. * pad * sinth;
    ypos = y1 + dist/2. * sinth + dytic - size * costh - 2. * pad * costh;
    vp_gtext(xpos,ypos,xpath,ypath,xup,yup,label);
}

int vp_optimal_scale(float chars, float min, float max, 
		     float *onum, float *dnum)
{
    int i, ntics, nopt;
    float num;
    char string[32];

    for (ntics = chars; ntics >= 1; ntics--) {
	nopt = optimal_scale(ntics, min, max, onum, dnum);
	
	for (i=0; i < nopt; i++) {
	    num = *onum + i*(*dnum);
	    if (fabsf(*dnum) > FLT_EPSILON && 
		fabsf(num) < FLT_EPSILON) num=0.;

	    sprintf(string,"%1.5g",num);
	    if ((strlen(string)+2)*nopt > chars) break;
	}
	    
	if (i == nopt) break;
	if (ntics > nopt) ntics = nopt;
    } 
    
    return nopt;
}

/* Algorithm for internal labeling from Tom Steppe, "Composing
 * well-tempered linear scales", Computer Language, September 1989,
 * 49-65. */
static int optimal_scale(int n, float min, float max, float *onum, float *dnum)
{
    int lo, hi;
    float d, nice;

    d = fabsf(max-min)/n;
    
    nice = nice_number(d);

    if (min <= max) {
	lo = (int) ceilf (min/nice);
	if ((lo-1)*nice+FLT_EPSILON >= min) lo--;
	
	hi = (int) floorf (max/nice);
	if ((hi+1)*nice <= max+FLT_EPSILON) hi++;

	*onum = lo * nice;
	*dnum = nice;
    } else {
	lo = (int) ceilf (max/nice);
	if ((lo-1)*nice+FLT_EPSILON >= max) lo--;
	
	hi = (int) floorf (min/nice);
	if ((hi+1)*nice <= min+FLT_EPSILON) hi++;

	*onum = hi * nice;
	*dnum = -nice;
    }
    
    return (hi-lo+1);
}

/* smallest nice number >= d */
static float nice_number (float d)
{
    float p, nice;
    int i, ie;
    const float set[] = {1.0, 2.0, 5.0};

    ie = (int) floor (log10f(d));
    nice = p = power (10.,ie);

    for (i=1; nice < d; i++) {
	if (i >= 3) {
	    i=0;
	    p *= 10.;
	}
	
	nice = set[i]*p;
    }

    return nice;
}

static float power (float f, int ie)
{
    float p;

    if (ie < 0) {
	f = 1./f;
	ie = -ie;
    }

    p = (ie & 1)? f: 1.; 
    
    for (ie >>= 1; ie > 0; ie >>=1) {
	f *= f;
	if (ie & 1) p *= f;
    }

    return p;
}

    
