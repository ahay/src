/* */
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

#include <rsf.h>
#include <rsfplot.h>

#define NARRAY 41

int main (void)
{
    float xval[NARRAY], zval[NARRAY];
    int plotfat=1, ith, iflag, emph;
    bool wantframe=false,option=true,animate=false;
    float xll=1.705,yll=1.37,xur=11.945,yur=8.87;
    float h=1.0,y=0.0,time=1.6,v=2.0,theta=30., xstart, zstart;
    float a, b, xs, xg, surface,slop,xmin,zmin,xmax,zmax;
    float xscale, zscale, dth, thtmp, xsp;
    float tanth, sinth, costh, d, xr, zr, x,z, x1,z1, xend,zend;

    vp_init();

    /*				set parameters */
    a = v * time /2.;
    b = sqrtf(a*a-h*h);
    xs = y - h;
    xg = y + h;
    surface = 0.;

    slop = 0.05*a;
    xmin = -2.;
    zmin = -1.;
    xmax =  2.;
    zmax = 2.;

    /* set vplot scales */
    vp_orig(xll,yur);
    vp_uorig(xmin,zmin);
    xscale = (xur-xll)/(xmax-xmin);
    zscale = -(yur-yll)/(zmax-zmin);
    vp_scale(xscale,zscale);

/* 					draw and annotate surface */
    vp_fat(2*plotfat);

    vp_color(7);   
    vp_umove(xmin,surface);
    vp_udraw(xmax,surface);

    vp_color(5);
    vp_fat(plotfat);
    vp_tjust(4,6);
    vp_utext(xs-0.5*slop,surface-slop,8,0,"S");
    vp_utext(xg,surface-slop,8,0,"G");
    vp_utext(y,surface-slop,8,0,"M");

    dth = 180./(NARRAY-1);

    vp_color(4);
    /* LOOP on theta (ith) to generate ellipse */
    for (ith = 0; ith < NARRAY; ith++) {
	thtmp = 90. - ith * dth;
        thtmp = thtmp*SF_PI/180.;
        tanth = tanf(thtmp);
        sinth = sinf(thtmp);
	costh = cosf(thtmp);

	d = hypotf(a*sinth,b*costh);

	xr  = y - d*sinth - h*h*costh*costh*sinth/d;
	zr  = surface + d*costh - h*h*costh*sinth*sinth/d;

	xval[ith] = xr;
	zval[ith] = zr;
    }				
	
/*					Draw ellipse */
    vp_color(5);
    vp_upline(xval,zval,NARRAY);

    if (option) {
/*					Loop on theta to draw SRG raypaths */

	iflag = 0;
	emph  = 3;

	for (ith = 0; ith < NARRAY; ith +=5) {
	    iflag++;
	    thtmp = 90. - ith * dth;
	    thtmp = thtmp*SF_PI/180.;
	    tanth = tanf(thtmp);
	    sinth = sinf(thtmp);
	    costh = cosf(thtmp);

	    d = hypotf(a*sinth,b*costh);

	    xr  = y - d*sinth - h*h*costh*costh*sinth/d;
	    zr  = surface + d*costh - h*h*costh*sinth*sinth/d;

/*					Shot to reflection point */
	    vp_color(4);

	    x = xs;
	    z = surface;
	    x1 = xr;
	    z1 = zr;
	    if(iflag == emph) vp_color(6);
	    if(ith != 1 && ith != NARRAY) vp_uarrow(x,z,x1,z1,0.02*zmax);
	    if(iflag == emph) vp_color(4);
/*					Annotate beta */
	    xend = (xs + xr)/2.;
	    zend = (surface + zr)/2.;

	    if(iflag == emph) {
		vp_color(6);
		vp_utext(xend-slop,  zend+0.8*slop,8,0,"\\F9 b");
		vp_color(4);
	    }

/*					reflection point to geophone */
	    x = xr;
	    z = zr;
	    vp_umove(x,z);
	    x1 = xg;
	    z1 = surface;
	    if(iflag == emph) vp_color(6);

	    if(ith != 1 && ith != NARRAY) vp_uarrow(x,z,x1,z1,0.02*zmax);

	    if(iflag == emph) vp_color(6);
/*					Annotate alpha */
	    xend = (xg + xr)/2.;
	    zend = (surface + zr)/2.;

	    if(iflag == emph) {
		vp_color(6);
		vp_utext(xend+slop,  zend+0.5*slop,8,0,"\\F9 a");
		vp_color(4);
	    }
	    if(animate) {
		vp_erase();	      /* clear screen, then redraw basic stuff */
/* 					redraw and annotate surface */
		vp_fat(2*plotfat);
		vp_color(7);
		vp_umove(xmin,surface);		
		vp_udraw(xmax,surface);

		vp_color(5);
		vp_fat(plotfat);
		vp_tjust(4,6);
		vp_utext(xs-0.5*slop,surface-slop,8,0,"S");
		vp_utext(xg,surface-slop,8,0,"G");
		vp_utext(y,surface-slop,8,0,"M");

/*					redraw ellipse */
		vp_color(5);
		vp_upline(xval,zval,NARRAY);
	    }					   
	}				/* end loop on theta */

    } else {
/*					Begin work for drawing dipping layer */
        theta *= SF_PI/180.;
        tanth = tanf(theta);
        sinth = sinf(theta);
	costh = cosf(theta);

	d = hypotf(a*sinth,b*costh);

	xr  = y - d*sinth - h*h*costh*costh*sinth/d;
	zr  = surface + d*costh - h*h*costh*sinth*sinth/d;

	xsp = y - h   - 2*(d-h*sinth)*sinth;
	xstart = xsp;
	zstart = surface + d/costh + tanth * (xstart-y);
	xend = 0.5*(xg + y);
	zend = surface + d/costh + tanth * (xend-y);


/* 					Draw dipping reflector */
	vp_color(2);
	vp_umove(xstart,zstart);
	vp_udraw(xend,zend);
	vp_dash(0.,0.,0.,0.);
/*					indicate angle theta */
	vp_color(7);
	vp_dash(.1,.1,.1,.1);
	vp_umove(xend,zend);
	vp_udraw(xend-0.5*h,zend);
	vp_dash(.0,.0,.0,.0);

/* 					finite-offset raypaths */
	vp_color(4);

/*					Shot to reflection point */
 	vp_uarrow(xs,surface,xr,zr,0.02*zmax);
	
/*					reflection point to geophone */
	vp_uarrow(xr,zr,xg,surface,0.02*zmax);
	
/* 					text */


        vp_color(5);
        vp_fat(plotfat);
        vp_tjust(4,6);
        vp_utext(xr,zr+2.0*slop,8,0,"R");
 	vp_utext(xend-5*slop,  zend-1.*slop,8,0,"\\F9 q");
    }

    if(wantframe) {
	vp_dash(0.,0.,0.,0.);
	vp_color(7);
	vp_fat(1);
	vp_umove(xmin,zmin);
	vp_udraw(xmin,zmax);
	vp_udraw(xmax,zmax);
	vp_udraw(xmax,zmin);
	vp_udraw(xmin,zmin);
    }

    return 0;
}

/* 	$Id$	 */
