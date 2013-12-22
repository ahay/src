/* Axis operations. */
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

#include "axis.h"
#include "vplot.h"

static float power (float f, int ie);
static float nice_number (float d);
static int optimal_scale(int n, float min, float max, 
			 /*@out@*/ float *onum, 
			 /*@out@*/ float *dnum);

void vp_simple_axis (float x1, float y1,     /* starting coordinates */
		     float x2, float y2,     /* ending coordinates */
		     float num1, float num2, /* number range */
		     float onum, float dnum, /* sampling */
		     float ltic              /* tick size */, 
		     char* label             /* axis label */, 
		     float size              /* symbol size */)
/*< Draw an axis.
...
Axis goes from (x1,y1) to (x2,y2) in vplot coordinates
its scale ranges from num1 to num2, ltic is the size of tick marks 
if dnum=0., an optimal linear scale is estimated
>*/
{
    int nopt, i, maxstrlen;
    float xpath, ypath, xup, yup, dist, costh, sinth, dtic, otic;
    float dxtic, dytic, loc, num, xpos, ypos;
    char string[32];
    const float pad = 0.15; /* between tic and number, number and label */
    const float aspect = 0.8;

    vp_bgroup("simple axis");
    vp_color (VP_WHITE); /* white */

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
	nopt = dist/(aspect*size);
    } else { /* use client's scale */
	nopt = 1+ (int) (floor) (num2-onum)/dnum;
	if (onum+(nopt-1)*dnum > num2) nopt--;
    }

    nopt = vp_optimal_scale(nopt,(bool) (dnum==0.),true, "%1.5g",
			    num1, num2, &onum, &dnum, &maxstrlen);

    /* figure out the tic mark spacing */
    otic = (onum-num1)*dist/(num2-num1);
    dtic = dnum*dist/(num2-num1); 

    vp_bgroup("tic marks");
    /* move to each tic mark location, draw the tic and the number */
    for (i=0; i < nopt; i++) {
	num = onum + i*dnum;
	if (fabsf(dnum) > SF_EPS && 
	    fabsf(num)  < SF_EPS) num=0.;

        snprintf(string,32,"%1.5g",num);
	loc = otic + i*dtic;
	xpos = x1 + loc * costh;
	ypos = y1 + loc * sinth;
        vp_move(xpos,ypos);
	vp_draw(xpos+dxtic,ypos+dytic);
	vp_gtext(xpos+dxtic+pad*sinth,
		 ypos+dytic-pad*costh,xpath,ypath,xup,yup,string);
    }
    vp_egroup();
    
    vp_bgroup("axis label");
    xpos = x1 + dist/2. * costh + dxtic + size * sinth + 2. * pad * sinth;
    ypos = y1 + dist/2. * sinth + dytic - size * costh - 2. * pad * costh;
    vp_gtext(xpos,ypos,xpath,ypath,xup,yup,label);
    vp_egroup();

    vp_egroup(); /* simple axis */
}

int vp_optimal_scale(int chars                /* characters */, 		
		     bool modify              /* modify the scale */,
		     bool parallel            /* parallel to the axis */,
		     const char* format       /* string format */,
		     float min, float max     /* scale range */, 
		     /*@out@*/ float *onum    /* output origin */,
		     /*@out@*/ float *dnum    /* output scaling */,
		     /*@out@*/ int *maxstrlen /* longest tick label */)
/*< Find an optimal scale. Returns the number of tics >*/
{
    int i, ntics, nopt, len;
    float num;
    char string[1024];

    nopt = 0;
    *maxstrlen = 1;

    for (ntics = chars; ntics >= 1; ntics--) {
	nopt = modify? optimal_scale(ntics, min, max, onum, dnum):chars;
	/* nopt - optimal number of tickmarks */

	*maxstrlen = 1;
	for (i=0; i < nopt; i++) {
	    num = *onum + i*(*dnum);
            /* number on a tickmark */

	    if (fabsf(*dnum) > SF_EPS && 
		fabsf(num) < SF_EPS) num=0.;
	    /* set it to zero if it is close to zero */

	    snprintf(string,1024,format,num);
	    /* put it in a string */
	    
	    len = strlen(string);
	    if (len > *maxstrlen) *maxstrlen = len;

	    /* see if we can fit it in */
	    if (modify) {
		if (parallel) {
		    if ((len+2)*nopt > chars) break;
		} else {
		    if (2.5*nopt > chars) break;
		}	    
	    }
	}
	    
	if (i == nopt) break;
	/* everything fits in, we are done */
	
	if (ntics > nopt) ntics = nopt;
    } 
    
    return nopt;
}

static int optimal_scale(int n, float min, float max, 
			 /*@out@*/ float *onum, /*@out@*/ float *dnum)
/* Algorithm for internal labeling from Tom Steppe, "Composing
 * well-tempered linear scales", Computer Language, September 1989,
 * 49-65. */
{
    int lo, hi;
    float d, nice;

    d = fabsf(max-min)/n;
    
    nice = nice_number(d);

    if (min <= max) {
	lo = (int) ceilf (min/nice);
	if ((lo-1)*nice+SF_EPS >= min) lo--;
	
	hi = (int) floorf (max/nice);
	if ((hi+1)*nice <= max+SF_EPS) hi++;

	*onum = lo * nice;
	*dnum = nice;
    } else {
	lo = (int) ceilf (max/nice);
	if ((lo-1)*nice+SF_EPS >= max) lo--;
	
	hi = (int) floorf (min/nice);
	if ((hi+1)*nice <= min+SF_EPS) hi++;

	*onum = hi * nice;
	*dnum = -nice;
    }
    
    return (hi-lo+1);
}

static float nice_number (float d)
/* smallest nice number >= d */
{
    float p, nice;
    int i, ie;
    const float set[] = {1.0, 2.0, 5.0};

    if (d == 0.0f) d=SF_EPS;

    ie = (int) floorf (log10f(d));
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
/* power function */
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

    
