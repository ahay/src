#include <math.h>
#include <stdio.h>

#include "axis.h"
#include "vplot.h"

float vp_opttic(float min, float max, float inch, float num0, float size) 
{
    const float eps=.000001, bignum=10000.;
    float num, dnum, dmax, div, b;
    int length, counter, len, i, div2;
    char string[100];
    const float sqr[3] = {1.414214, 3.162278, 7.071068};
    const float vint[4] = {1., 2., 5., 10.};
    
    if (fabsf(max - min) < eps) return 0.5;
    
    div = fabsf(max-min)/inch;

    if (div < eps) return (max > min)? 0.5: -0.5;
    
    div2 = log10f (div);
    if (div < 1.) div2--;
    b = div / powf (10.,(float) div2);
    
    for (i = 0; i < 3; i++) {
	if (sqr[i] > b) break;
    }
    
    for (dnum = vint[i] * powf (10.,div2); dnum < dmax; dnum *= 2.) {
	length = 0;
	counter = (int) (fabsf(max - min)/dnum);
	for (i= 0; i < counter; i++) {
	    num = num0 + i * dnum;
	    if (fabsf(num) < fabsf(max - min)/bignum) num = 0.;
	    sprintf (string, "%1.5g", num);
	    len = strlen (string);
	    if (len > length) length = len;
	}
	dmax = (length + 1.5) * fabsf(max - min)*size/(33.* inch);
    }

    return dnum;
}

void vp_simpleaxis (float x1, float y1, 
		    float x2, float y2, 
		    float num1, float num2,
		    float dnum, float ltic, char* label, float size)
{
    float ch, xpath, ypath, xup, yup, dist, costh, sinth, dtic, xpos, ypos;
    float dxtic, dytic, loc, num, pad, dnumtic;
    char string[10];

    /* pad is the space skipped between tic and number, and number and label */
    pad = 0.15;
    ch = size / 33.;

    vp_color (7);
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
        xpath = ch * costh;
        ypath = ch * sinth;
        xup = -ch * sinth;
        yup = ch * costh;
    } else {	
        vp_tjust(TH_CENTER,TV_BOTTOM);
        xpath = -ch * costh;
        ypath = -ch * sinth;
        xup = ch * sinth;
        yup = -ch * costh;
    } 

    dnumtic = dnum;
    if (0. == dnumtic) dnumtic = vp_opttic(num1,num2,dist,num1,size);
    if (num1 > num2) dnumtic = -dnumtic;

    /* figure out the tic mark spacing */
    dtic = dnumtic/(num2-num1) * dist; 
    if (dtic < 0.) dtic *= -1.;

    /* move to each tic mark location, draw the tic and the number */
    for (loc=0., num=num1; loc < dist; loc+=dtic, num+=dnumtic) {
        sprintf(string,"%1.5g",num);
	xpos = x1 + loc * costh;
	ypos = y1 + loc * sinth;
        vp_move(xpos,ypos);
	vp_draw(xpos+dxtic,ypos+dytic);
	vp_gtext(xpos+dxtic+pad*sinth,
		 ypos+dytic-pad*costh,xpath,ypath,xup,yup,string);
    }

    /* now the axis label */
    xpos = x1 + loc/2. * costh + dxtic + ch * sinth + 2. * pad * sinth;
    ypos = y1 + loc/2. * sinth + dytic - ch * costh - 2. * pad * costh;
    vp_gtext(xpos,ypos,xpath,ypath,xup*1.5,yup*1.5,label);
}
