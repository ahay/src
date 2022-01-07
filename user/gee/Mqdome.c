/* 3-D synthetic image from Jon Claerbout. */
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

#include "random.h"

int main (int argc, char* argv[])
{
    bool impedance;
    float ranget,o1,d1,o2,d2,o3,d3;
    float gaussvel, thro, rand, t1, t2;
    float frac1, frac2, f1, zrad, xnew;
    float thetac, dr, x, y, z, rmax;
    float xp, yp, rcap, rad2, dz;
    float r, t, frac, zcap, ynew, interval;
    int large, n1, n2, n3, endtaper;
    int i1,  i2, i3, n1a, n1b, it1, it2, i,slicei;
    float *imp, *refl, *refl1, *refl2, ***earth;
    sf_triangle tr;
    sf_file mod, trace;

    sf_init(argc, argv);
    mod = sf_output("out");
    sf_setformat(mod,"native_float");

    if (!sf_getfloat("ranget",&ranget)) ranget=3.;
    if (!sf_getint("n1",&n1)) n1=400;
    if (!sf_getint("n2",&n2)) n2=100;
    if (!sf_getint("n3",&n3)) n3=50;
    sf_putint(mod,"n1",n1);
    sf_putint(mod,"n2",n2);  
    sf_putint(mod,"n3",n3);

    if (!sf_getint("large",&large)) large = (int) n1*ranget;

    if (!sf_getfloat("d1",&d1)) d1=0.004; 
    if (!sf_getfloat("d2",&d2)) d2=0.01;  
    if (!sf_getfloat("d3",&d3)) d3=0.02;  
    
    if (!sf_getfloat("o1",&o1)) o1=0.;    
    if (!sf_getfloat("o2",&o2)) o2=-3.*d2;
    if (!sf_getfloat("o3",&o3)) o3=-3.*d3;

	sf_putfloat(mod,"d1",d1);
	sf_putfloat(mod,"d2",d2);
	sf_putfloat(mod,"d3",d3);

	sf_putfloat(mod,"o1",o1);
	sf_putfloat(mod,"o2",o2);
	sf_putfloat(mod,"o3",o3);
		
    if (!sf_getfloat("gaussvel",&gaussvel)) gaussvel=2.5;
    if (!sf_getfloat("throw",&thro)) thro=0.01;
    if (!sf_getint("endtaper",&endtaper)) endtaper=20;
    if (!sf_getint("slicei",&slicei)) slicei=40;

    if (!sf_getbool("impedance",&impedance)) impedance=false;

    sf_putstring(mod,"label1","Time (s)");
    sf_putstring(mod,"label2","West-East (km)");
    sf_putstring(mod,"label3","South-North (km)");

    if (NULL != sf_getstring("trace")) {
	/* file to optionally output the master trace */
	trace = sf_output("trace");
	sf_putint(trace,"n1",large);
	sf_putfloat(trace,"d1",d1);
	sf_putfloat(trace,"o1",o1);
	sf_setformat(trace,"native_float");
    } else {
	trace = NULL;
    }

    imp = sf_floatalloc (large);
    earth = sf_floatalloc3 (n1,n2,n3);
    refl = sf_floatalloc (n1);

    random_init (19921995);

    imp[0] = 0.;
    for (i1=1; i1 < large; i1++) {
	rand = random0 ();
	if (rand > .3) { 
	    imp[i1] = imp[i1-1];
	} else { 
	    rand = random0 ();
	    imp[i1] = .1 * rand;
	}
    }

    for (i1=1; i1 < large; i1++) {
	imp[i1] += imp[i1-1];
    }

    /* optionally output the master trace */
    if (NULL != trace) sf_floatwrite(imp,large,trace);

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		earth[i3][i2][i1] = 1.0;
	    }
	}
    }

    thetac = 45 * SF_PI/180.;
    dr = fabsf(d2) + fabsf(d1);
    x = o2+d2*(n2-1);
    y = o3+d3*(n3-1);
    rmax = hypotf(x,y) + dr;
    xp = 0.75 * x;
    yp = 0.75 * y;
    rcap = 0.3 * rmax;
    /* rad1 = (xp*xp + yp*yp) - rcap; */
    rad2 = (xp*xp + yp*yp);
    n1a = .15 * n1;
    n1b = .8  * n1;
    dz = 1./(n1b-n1a);

    for (i3=0; i3 < n3; i3++) {
	y = o3+d3*i3;
	for (i2=0; i2 < n2; i2++) {
	    x = o2+d2*i2;
	    r = hypotf(x,y)  + dr/2.;

	    /* shallow beds */      
	    for (i1=0; i1 < n1a-1; i1++) {
		earth[i3][i2][i1] = 1. + (imp[i1+1] - imp[i1]);
	    }
	    
	    /* deep beds */
	    for (i1=n1b; i1 < n1; i1++) {
		t = i1+1. +.3 * (i2+1.) +.2 * (i3+1.);
		it1  = t;
		it2  = t+1;
		frac = t - it1;
		interval = 
		    ((imp[it2-1] * (1.-frac) + imp[it2] * frac) - 
		     (imp[it1-1] * (1.-frac) + imp[it1] * frac));
		earth[i3][i2][i1] = 1.0   + interval;
	    }

	    for (i1=n1a-1; i1 < n1b; i1++) {
		z = ( i1+1. - n1a) / (n1b - n1a); /* 0 to 1 */
		zcap = z * gaussvel /4.0;
		z = dz + (1-2*dz) * z;   /* dz to 1-dz */

		t1 = (large/gaussvel) * (r/rmax) / sqrtf( -logf( z-dz/2));
		t2 = (large/gaussvel) * (r/rmax) / sqrtf( -logf( z+dz/2));

		t1 = (large/gaussvel) * (r/rmax) * sqrtf( -logf( z+dz/2));
		t2 = (large/gaussvel) * (r/rmax) * sqrtf( -logf( z-dz/2));

		if( t1 > t2) sf_error("misconception");

		xnew = cosf(thetac) * (x-xp) + sinf(thetac) * (y-yp);
		xnew = xnew * 0.75;
		ynew = -sinf(thetac) * (x-xp) + cosf(thetac) * (y-yp);
		if ( rcap*rcap > ( xnew*xnew + ynew*ynew + zcap*zcap )  ) {
		    if (r < rad2) 
			zrad = (rad2 - r)/rcap; 
		    else     
			zrad = 0.0;  /* no faulting after center */ 
		    t1  +=  thro * n1 * zrad; 
		    t2  +=  thro * n1 * zrad; 
		}
		it1  = t1;
		it2  = t2;
		if( it1 >=1 && it1+1 < large) {
		    frac1 = t1 - it1;
		    if( it2 >=1 && it2+1 < large) {
			frac2 = t2 - it2;
			interval = 
			    ((imp[it2-1] * (1.-frac2) + imp[it2] * frac2) - 
			     (imp[it1-1] * (1.-frac1) + imp[it1] * frac1));
			earth[i3][i2][i1] =  1.0   + interval / (t2 - t1);
		    }
		}
	    }
	}
    }

    /*   Form reflectivity from impedance. */
    if (!impedance) {
	refl[0] = 0.;
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		refl2 = earth[i3][i2]+1;
		refl1 = earth[i3][i2];
		for (i1=0; i1 < n1-1; i1++) {
		    refl[i1+1] = 
			(refl1[i1] - refl2[i1]) / 
			(refl1[i1] + refl2[i1]);
		}
		for (i1=0; i1 < n1; i1++) {
		    earth[i3][i2][i1] = refl[i1];
		}
	    }
	}
    }

    /* temporal smoothing */
    for (i = 2; i <= 4; i++) {
	tr = sf_triangle_init (i, n1,false);
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		sf_smooth2(tr,0,1,impedance,earth[i3][i2]);
	    }
	}
	sf_triangle_close (tr);
    }

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    for (i1=n1; i1 >= 1; i1--) {
		earth[i3][i2][i1] -= earth[i3][i2][i1-1];
	    }
	    earth[i3][i2][0] = 0.;
	}
    }
    
    f1 = (float) endtaper;
    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < endtaper; i1++) {
		earth[i3][i2][i1] *= (i1+1.)/f1;
		earth[i3][i2][n1-1-i1] *= (i1+1.)/f1;
	    }
	    sf_floatwrite(earth[i3][i2],n1,mod);
	}
    }

    exit (0);
}

/* 	$Id$	 */
