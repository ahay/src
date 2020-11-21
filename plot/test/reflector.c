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

#include <stdlib.h>
#include <math.h>

#include <rsf.h>
#include <rsfplot.h>

#define NANGLES 11
#define NARCS2 6

int main(void)
{
    bool wantframe=false;
    float xll=1.705,yll=1.37,xur=11.945,yur=8.87, plotfat=1.;
    float xstart,zstart,xend,zend,xscale,zscale, y=0.5;
    float xmin,zmin,xmax,zmax,slop,surface=0.;
    float d=1.4,xs,zs,xr,zr,xc;
    float degrad,theta=30.,tanth,sinth,costh;

#ifdef REFLKINE
    float x,z, zc=0.;
#endif

#ifdef TUCHEL
    float xp,zp;
#endif

#ifdef REFLEXPT
    int narcs=3, iarc,iang;
    float drad,arcmin,phi,dth,radius,arcln,zr2,xr2,d2,zs2,xs2,xsp,zsp;
    float xval[NARCS2][NANGLES], zval[NARCS2][NANGLES];
#endif

#if defined(TUCHEL) || defined(REFLEXPT)
    float xs1, zs1, d1, xr1, zr1;
#endif

    /*		 			Origin is at midpoint */
    degrad = acosf(-1.)/180.;
    theta *= degrad;
    tanth = tanf(theta);
    sinth = sinf(theta);
    costh = cosf(theta);

#ifdef REFLEXPT
    arcln=45.;
    dth = arcln*degrad/(NANGLES-1);
    arcmin = -arcln / 2. * degrad;
    drad = d / (narcs + 1);
#endif

    xs = y;					/* Shotpoint */
    zs = surface;
    xr = xs - d*sinth; 			        /* Reflection point */
    zr = zs + d*costh; 
    xc = xr;					/* Migrated location */

#if defined(TUCHEL) || defined(REFLEXPT)
    xs1= y + d/2.;				/* Second shotpoint */
    zs1= surface;
    d1 = d + (xs1-xs)*sinth;			/* Second Reflection point */
    xr1= xs1 - d1*sinth; 			
    zr1= zs1 + d1*costh; 
#endif

#ifdef REFLEXPT
    xsp = y  - 2*d*sinth;			/* Image shotpoint */
    zsp = surface + 2*d*costh;

    xs2= y - d/2.;				/* Third shotpoint */
    zs2= surface;

    d2 = d + (xs2-xs)*sinth;			/* Third Reflection point */

    xr2= xs2 - d2*sinth; 			
    zr2= zs2 + d2*costh; 
#endif

#ifdef TUCHEL
    /* For Tuchel's Law */
    xp = xs + (xs1-xs)*costh*costh;
    zp = zs + (xs1-xs)*costh*sinth;
#endif

    slop = 0.05*d;
    xmin = -2.;
    zmin = 0.;
    xmax =  2.;
    zmax = 3.;
/*
    xcenter=0.5*(xmin+xmax);
    zcenter=0.5*(zmin+zmax);
*/

/*						dipping-reflector parameters */
    xstart = xc - d;
    if (xstart < xmin) xstart = xmin;
    xend = xc + d;
    if (xend > xmax) xend = xmax;
    zstart = surface + d/costh + tanth * (xstart-y);
    zend = surface + d/costh + tanth * (xend-y);

/*						Pre-compute arcs of circles */

/*						Downgoing wavefronts */
#ifdef REFLEXPT
    radius = 0.;
    for (iarc = 0; iarc < narcs; iarc++) {
	radius += drad;
	for (iang = 0; iang < NANGLES; iang++) {
	    phi = arcmin + iang * dth;
	    xval[iarc][iang] = xs - radius * sin(theta+phi);
	    zval[iarc][iang] = zs + radius * cos(theta+phi);
	}
    }
/*						Upcoming wavefronts */
    radius = d + drad/2.;
    for (iarc = narcs; iarc < NARCS2; iarc++) {
	radius += drad;
	for (iang = 0; iang < NANGLES; iang++) {
	    phi = arcmin + iang * dth;
	    xval[iarc][iang] = xsp + radius * sin(theta+phi);
	    zval[iarc][iang] = zsp - radius * cos(theta+phi);
	}
    }
#endif

    vp_init();

/* 						set vplot scales */
    vp_orig(xll,yur);
    vp_uorig(xmin,zmin);
    xscale = (xur-xll)/(xmax-xmin);
    zscale = -(yur-yll)/(zmax-zmin);
    vp_scale(xscale,zscale);
/*						Text centering */
    vp_tjust(4,1);

/* 					draw surface */
    vp_fat(2*plotfat);

    vp_color(7);
    vp_umove(xmin,surface);
    vp_udraw(xmax,surface);

/*					draw dipping reflector */
    vp_color(2);   
    vp_umove(xstart,zstart);
    vp_udraw(xend,zend);
    vp_dash(0.,0.,0.,0.);

/*					indicate angle theta */
    vp_color(5);
    vp_utext(xend-5*slop,zend-0.5*slop,8,0,"\\F9 q");

    vp_color(7);
    vp_dash(.1,.1,.1,.1);
    vp_umove(xend,zend);
    vp_udraw(xend-0.5*d,zend);
    vp_dash(.0,.0,.0,.0);

/* 					main zero-offset raypath */
    vp_color(4);
/*					Shot to reflection point */
    vp_uarrow(xs,surface,xr,zr,0.02*zmax);
	
/*					reflection point to geophone */
    vp_uarrow(xr,zr,xs,surface,0.02*zmax);
	
#ifdef REFLEXPT
    vp_color(5);
    vp_fat(plotfat);
    vp_utext(xs, zs -slop,8,0,"S\\s60 \\_2");
    vp_utext(xs1,zs1-slop,8,0,"S\\s60 \\_3");
    vp_utext(xs2,zs2-slop,8,0,"S\\s60 \\_1");
    vp_utext(xr,zr+3.0*slop,8,0,"R\\s60 \\_2");
    vp_utext(xr1,zr1+3.0*slop,8,0,"R\\s60 \\_3");
    vp_utext(xr2,zr2+3.0*slop,8,0,"R\\s60 \\_1");
/*					  second shotpoint */
    vp_color(4);
    vp_uarrow(xs1,zs1,xr1,zr1,0.02*zmax);
    vp_uarrow(xr1,zr1,xs1,zs1,0.02*zmax);

/*					  third shotpoint */
    vp_color(4);
    vp_uarrow(xs2,zs2,xr2,zr2,0.02*zmax);
    vp_uarrow(xr2,zr2,xs2,zs2,0.02*zmax);

/*					  Draw Wavefronts */
    vp_color(6);
    for (iarc=0; iarc < NARCS2; iarc++) {
	vp_upline(xval[iarc],zval[iarc],NANGLES);
    }

#endif

#ifdef REFLKINE				      
/* reflector kinematics */
    vp_color(5);
    vp_fat(plotfat);
    vp_utext(xs,surface-slop,8,0,"S");
    vp_utext(xc,surface-slop,8,0,"C");
    vp_utext(xr,zr+3.0*slop,8,0,"R");
    x=0.5*(xr+y);
    z=0.5*(zr+surface);
    vp_utext(x-slop,z,8,0,"d");
    vp_dash(.1,.1,.1,.1);	/* dashed vertical line */
    vp_umove(xc,zc);
    vp_udraw(xr,zr);
    x=xc; 
    z=0.5*(zc+zr);
    vp_utext(x-slop,z,8,0,"z");
#endif

#ifdef TUCHEL					
/* Tuchel's law */
    vp_color(5);
    vp_fat(plotfat);
    vp_utext(xs, zs -slop,8,0,"S");
    vp_utext(xs1,zs1-slop,8,0,"S\\F15 \\v131 \\-");
    vp_utext(xr,zr+3.0*slop,8,0,"R");
    vp_utext(xr1,zr1+3.0*slop,8,0,"R\\F15 \\v131 \\-");
    vp_utext(xp+1.0*slop,zp+1.0*slop,8,0,"P");
    x=0.5*(xr+y);
    z=0.5*(zr+surface);
    vp_utext(x-slop,z,8,0,"d");
/*					  second shotpoint */
    vp_color(4);
    vp_uarrow(xs1,zs1,xr1,zr1,0.02*zmax);
    vp_uarrow(xr1,zr1,xs1,zs1,0.02*zmax);

    vp_dash(.1,.1,.1,.1);	/* dashed wavefront line */
    vp_color(6);
    vp_umove(xs,zs);
    vp_udraw(xp,zp);
#endif

/* draw a frame */
    if (wantframe) {
	vp_dash(0.,0.,0.,0.);
	vp_color(7);
	vp_fat(1);
	vp_umove(xmin,zmin);
	vp_udraw(xmin,zmax);
	vp_udraw(xmax,zmax);
	vp_udraw(xmax,zmin);
	vp_udraw(xmin,zmin);
    }

    exit(0);
}

/* 	$Id$	 */
