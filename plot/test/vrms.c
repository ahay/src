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
#include <stdlib.h>

#include <rsfplot.h>

int main(void)
{
    bool wantframe=false;
    float xll=1.705,yll=1.37,xur=11.945,yur=8.87,theta3=20.;
    float sn[3],cs[3],vv[3]={1.5,2.0,2.5},dz[3]={1.0,2.0,1.5};
    float zbottom, h, t, x, z, x1,z1, xtxt, ztxt;
    float xmin, zmin, xmax, zmax, xscale, zscale;
    int i, plotcol=6,plotfat=1;

    vp_init();

/* precompute quantities related to vplot scales */
    zbottom = 0.;
    for (i = 0; i < 3; i++) {
	zbottom += dz[i];
    }
    theta3 *= acosf(-1.)/180.;
    sn[2]=sinf(theta3);

    h=0.; 
    t=0.;
    for (i=2; i >= 0; i--) {
	if (i != 2)  sn[i] = sn[i+1]*vv[i]/vv[i+1]; /* snell's law */
	cs[i] = sqrtf(1. - sn[i]*sn[i]);
	h += dz[i]*sn[i]/cs[i];
	t += 2*dz[i]/(cs[i]*vv[i]);
    }

    xmin = -h - 0.1*zbottom; if (xmin > -2.) xmin=-2.;
    zmin = -0.1 * zbottom;
    xmax = -xmin;
    zmax = 1.1*zbottom; if (zmax < 3.) zmax=3.;

/* set vplot scales */
    vp_orig(xll,yur);
    vp_uorig(xmin,zmin);
    xscale = (xur-xll)/(xmax-xmin);
    zscale = -(yur-yll)/(zmax-zmin);
    vp_scale(xscale,zscale);

/* draw flat reflectors */
    vp_color(plotcol);
    vp_fat(plotfat);

    z=0.;
    for (i=0; i < 3; i++) {
	z += dz[i];
	vp_umove(xmin,z);
	vp_udraw(xmax,z);
    }

/* draw vertical zero-offset ray */
    vp_color(7);
    vp_umove(0.,0.);
    vp_udraw(0.,zbottom);

/* draw ray paths from bottom reflector to geophone */
    vp_color(2);

    x=0; 
    z=zbottom;
    for (i=2; i >= 0; i--) {	
	if (i != 0) {
	    vp_umove(x,z);
	    x += dz[i]*sn[i]/cs[i];
	    z -= dz[i];
	    vp_udraw(x,z);
	} else  {
	    x1 = x; 
	    z1 = z; 
	    x += dz[i]*sn[i]/cs[i];
	    z -= dz[i]; 
	    vp_uarrow(x1,z1,x,z,0.02*zbottom);
	}  
    }

/* draw ray paths from source to bottom reflector */
    vp_color(2);
    x=0; 
    z=zbottom;
    for (i=2; i >= 0; i--) {
	if (i != 2) {
	    vp_umove(x,z);
	    x -= dz[i]*sn[i]/cs[i];
	    z -= dz[i];
	    vp_udraw(x,z);
	} else  {
	    x1 = x; 
	    z1 = z; 
	    x -= dz[i]*sn[i]/cs[i];
	    z -= dz[i]; 
	    vp_uarrow(x,z,x1,z1,0.02*zbottom);
	}  
    }

/* draw surface */
    vp_color(7);
    vp_fat(2*plotfat);
    vp_umove(xmin,0.);
    vp_udraw(xmax,0.);

/* draw text */
    vp_color(5);
    vp_fat(plotfat);
    vp_tjust(2,3);

    x = 0.;
    z = zbottom * 1.05;
    vp_utext(x,z,8,0,"R");

    x = -h;
    z = 0.- zbottom *.05;
    vp_utext(x,z,8,0,"S");

    x = h;
    z = 0. - zbottom *.05;
    vp_utext(x,z,8,0,"G");

    x = -h;
    z = 0. + dz[0]/2.;
    vp_utext(x,z,8,0,"v\\s60 \\_1");

    z += dz[0]/2. + dz[1]/2.;
    vp_utext(x,z,8,0,"v\\s60 \\_2");

    z += dz[1]/2. + dz[2]/2.;
    vp_utext(x,z,8,0,"v\\s60 \\_3");

    x = 0.;
    z = 0. - zbottom *.05;
    vp_utext(x,z,8,0,"2h");
    x = 0. + 0.05*zbottom;
    vp_umove(x,z);
    x1 = h - 0.05 * zbottom;
    vp_uarrow(x,z,x1,z,-0.02*zbottom);
    x = 0. - 0.05*zbottom;
    vp_umove(x,z);
    x1 = -h + 0.05 * zbottom;
    vp_uarrow(x,z,x1,z,-0.02*zbottom);
	
/*						Highlight the trigonometry */
    vp_color(7);
    vp_fat(plotfat);
/*						bottom of triangle */
    z=zbottom-dz[2];
    x=dz[2]*sn[2]/cs[2];
    vp_umove(x,z);
    x1 = x + dz[1]*sn[1]/cs[1];
    z1 = z;
    vp_udraw(x1,z1);
    xtxt = (x+x1)/2.;
    ztxt = z + 0.05 * zbottom;
    vp_utext(xtxt,ztxt,8,0,"\\F9 D\\F3 x\\_\\s60 i");

/*						vertical side of triangle */
    x = x1;
    z = z1 - dz[1];
    vp_umove(x1,z1);
    vp_udraw(x,z);
    xtxt = x + 0.05 * zbottom;
    ztxt = (z+z1)/2.;
    vp_utext(xtxt,ztxt,8,0,"\\F9 D\\F3 z\\_\\s60 i");
    xtxt = x - sn[1]/cs[1]*dz[1]/8.;
    ztxt = z + dz[1]/4.;
    vp_utext(xtxt,ztxt,8,0,"\\F9 q\\F3 \\s60 \\_i");

/*						vertical side in tau space */
    xtxt = 0. - 0.05*zbottom;
    ztxt = z + dz[1]/2.;
    vp_utext(xtxt,ztxt,8,90,
	     "\\F9 Dt\\F3 \\s60 \\_i\\^\\s100 v\\s60 \\_i\\^\\s100 /2");

/*						hypotenuse of triangle */
    z1=zbottom-dz[2];
    x1=dz[2]*sn[2]/cs[2];
    vp_umove(x,z);
    vp_udraw(x1,z1);
    xtxt = (x+x1)/2. - 0.05 *zbottom;
    ztxt = (z+z1)/2.;
    vp_utext(xtxt,ztxt,8,60,
	     "\\F9 D\\F3 t\\s60 \\_i\\^\\s100 v\\s60 \\_i\\^\\s100 /2");

/* draw a frame */

    if(wantframe) {
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

/* 	$Id: vrms.c 7107 2011-04-10 02:04:14Z ivlad $	 */
