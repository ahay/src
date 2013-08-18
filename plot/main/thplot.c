/* Hidden-line surface plot.

Takes: > plot.vpl */
/*
  Copyright (C) 1987 The Board of Trustees of Stanford University
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
/* Modified from the original version by Shuki Ronen. */

#include <math.h>
#include <float.h>

#include <rsf.h>
#include <rsfplot.h>

int main(int argc, char* argv[])
{
    int nx, ny, n3, titlsz, titlefat, plotfat, plotcolup, plotcoldn;
    int i, i3, i2, i1, axissz, axisfat, gainstep;
    float *max, *min, ***fff, *ff, alpha, alpha2, zc, sz, tt, ee, gg, dx, dy;
    float pclip, clip, pbias, cosa, sina, xc, scalex, xs, ratio;
    float ax[4], ay[4], dclip, xlen, zlen, zmax, zmin, y, fm, scalez, ys;
    float dz, ox, oy, x2, x1, f2, f1, old, z, x, f, r, y2;
    const float eps=1.e-5;
    size_t len;
    char *label, *unit, *labels[3], *title, key[8];
    bool uflag, dflag, norm, axis, axis1, axis2, axis3, wanttitle;
    sf_file in=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    vp_init();

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&ny)) ny=1;
    n3 = sf_leftsize(in,2);

    max = sf_floatalloc (nx+1);
    min = sf_floatalloc (nx+1);
    fff = sf_floatalloc3 (nx,ny,n3);
    ff = fff[0][0];

    sf_floatread(ff,nx * ny * n3,in);

    if (!sf_getbool("uflag",&uflag)) uflag=true;
    /* if y, plot upper side of the surface */
    if (!sf_getbool("dflag",&dflag)) dflag=true;
    /* if y, plot down side of the surface */

    if (!sf_getfloat("alpha",&alpha)) alpha=45.;
    /* apparent angle in degrees, |alpha| < 89 */
    alpha2 = alpha * SF_PI / 180.;
    cosa = cosf (alpha2);
    sina = sinf (alpha2);

    if (!sf_getint("titlsz",&titlsz)) titlsz=9;
    /* title size */
    if (!sf_getint("axissz",&axissz)) axissz=6;
    /* axes size */
    if (!sf_getint("plotfat",&plotfat)) plotfat=0;
    /* line fatness */
    if (!sf_getint("titlefat",&titlefat)) titlefat=2;
    /* title fatness */
    if (!sf_getint("axisfat",&axisfat)) axisfat=2;
    /* axes fatness */
    if (!sf_getint("plotcolup",&plotcolup)) plotcolup=VP_YELLOW;
    /* color of the upper side */
    if (!sf_getint("plotcoldn",&plotcoldn)) plotcoldn=VP_RED;
    /* color of the lower side */

    for (i=0; i < 3; i++) {
	snprintf(key,8,"label%d",i+1);
	if (NULL == (label = sf_getstring(key))) 
	    label = sf_histstring(in,key);
        /*< label# label on #-th axis >*/ 
	snprintf(key,7,"unit%d",i+1);
	if (NULL == (unit = sf_getstring(key))) 
	    unit = sf_histstring(in,key);
	/*< unit# unit on #-th axis >*/
	if (NULL != label && NULL != unit) {
	    len = strlen(label)+strlen(unit)+4;
	    labels[i] = sf_charalloc(len);
	    snprintf(labels[i],len,"%s (%s)",label,unit);
	    free(label);
	    free(unit);
	} else {
	    labels[i] = label;
	}
    }

    if (!sf_getbool("wanttitle",&wanttitle)) wanttitle=true;

    if (wanttitle && 
	NULL == (title = sf_getstring("title")) && 
	NULL == (title = sf_histstring(in,"title"))) {
	title = sf_histstring(in,"in");
    } else {
	title = NULL;
    }

    if (sf_getfloat("tpow",&tt) && 0. != tt) {
	/*< tpow=0 time power gain >*/
	for (i2=0; i2 < ny*n3; i2++) {
	    for (i1=0; i1 < nx; i1++) {
		fff[0][i2][i1] *= powf((float) i1,tt);
	    }
	}
    }

    if (sf_getfloat("epow",&ee) && 0. != ee) {
	/*< epow=0 exponential gain >*/
	for (i2=0; i2 < ny*n3; i2++) {
	    for (i1=0; i1 < nx; i1++) {
		fff[0][i2][i1] *= expf(i1*ee);
	    }
	}
    }

    if (sf_getfloat("gpow",&gg) && 1. != gg) {
	/*< gpow=1 power gain >*/
	for (i1=0; i1 < nx*ny*n3; i1++) {
	    ff[i1] = powf(ff[i1],gg);
	}
    }

    if (!sf_getbool("axis",&axis)) axis=true;
    if (!sf_getbool("axis1",&axis1)) axis1=true;
    if (!sf_getbool("axis2",&axis2)) axis2=true;
    if (!sf_getbool("axis3",&axis3)) axis3=true;
    /* plot axis */

    if (!sf_getfloat ("clip",&clip)) clip = 0.;
    /* data clip */

    if (0. == clip) {
	if (!sf_getfloat ("pclip",&pclip)) pclip = 100.;
	/* data clip percentile */
	if (pclip <=0. || pclip > 100.)
	    sf_error("pclip=%g should be > 0 and <= 100",pclip);

	if (!sf_getint("gainstep",&gainstep)) gainstep=0.5+nx/256.;
	/* subsampling for gpow and clip estimation */
	if (gainstep <= 0) gainstep=1;

	if (!sf_getfloat("bias",&pbias)) pbias=0.;
	/* subtract bias from data */

	vp_gainpar (NULL,fff[0],nx,ny,gainstep,
		    pclip,100.,&clip,&gg,false,&pbias,n3,-2,0);
    }

    if (!sf_getfloat ("dclip",&dclip)) dclip=1.;
    /* change the clip: clip *= dclip */
    clip *= dclip;

    if (0.==clip) clip=FLT_EPSILON;

    if (!sf_getbool("norm",&norm)) norm=true;
    /* normalize by the clip */

    for (i1=0; i1 < nx*ny*n3; i1++) {
	ff[i1] /= clip;
	if (norm) {
	    if (ff[i1] > 1.) {
		ff[i1] = 1.; 
	    } else if (ff[i1] < -1.) {
		ff[i1] = -1.;
	    }
	} 
    }

    zlen=VP_STANDARD_HEIGHT;
    xlen=VP_STANDARD_HEIGHT/VP_SCREEN_RATIO;

    if (!sf_getfloat("xc",&xc)) xc=1.5;
    if (!sf_getfloat("zc",&zc)) zc=3;
    /* lower left corner of the plot */

    if (!sf_getfloat("ratio",&ratio)) ratio=5.;
    /* plot adjustment */

    dy = SF_MIN((xlen-2*xc)/((nx+ny)*cosa),
		(zlen-2*zc)/(ny*sina));
    dx = 1./ny;
    scalex = dy*cosa/dx;

    zmax = -SF_HUGE;
    zmin = SF_HUGE;
    for (i3 = 0; i3 < n3; i3++) {
	for (i2 = 0; i2 < ny; i2++) {
	    for (i1 = 0; i1 < nx; i1++) {
		y = i2 * dy;
		fm = fff[i3][i2][i1] + y * sina;
		if (fm < zmin) {
		    zmin = fm;
		} 
		if (fm > zmax) {
		    zmax = fm;
		}
	    }
	}
    }
    sf_getfloat("zmax",&zmax);
    sf_getfloat("zmin",&zmin);

    if (zmax == zmin) zmax = zmin+FLT_EPSILON;

    if (!sf_getfloat("sz",&sz)) sz=6.;
    /* vertical scale */

    scalez = sz / (zmax - zmin);
    dz = dy * sina / scalez;
    vp_scale (scalex, scalez);

    vp_style (VP_STANDARD);
    vp_fat (plotfat);

    for (i3 = 0; i3 < n3; i3++) {
	if (i3 > 0) vp_erase();

	/* initial a  */
	for (i1 = 0; i1 <= nx; i1++) {
	    max[i1] = -SF_HUGE;
	    min[i1] = SF_HUGE;
	}

	/* make a pseudo 3-D plot */
	for (i2 = 0; i2 < ny; i2++) {
	    y = i2 * dy;

	    oy = zc + y * sina;
	    ox = xc + y * cosa;
	    vp_orig (ox, oy);

	    if (i2 == 0) {
		vp_umove (0., 0.);
		vp_where (ax, ay);
		vp_umove (nx*dx, 0.);
		vp_where (ax + 2, ay + 2);
		vp_umove (0., 1.);
		vp_where (ax + 3, ay + 3);
	    } else if (i2 == ny - 1) {
		vp_umove (0., 0.);
		vp_where (ax + 1, ay + 1);
	    }

	    if (dflag) {
		vp_color(plotcoldn);

		for (i1 = 0; i1 < nx; i1++) {
		    min[i1] = min[i1+1] - dz;
		}

		f2 = fff[i3][i2][0];
		x2 = 0.;
		vp_umove (f2, x2);
		
		old = f2;
		for (i1 = 1; i1 < nx; i1++) {
		    x1 = x2;
		    f1 = f2;
		    x2 = dx * i1;
		    f2 = fff[i3][i2][i1];

		    if (f1 < min[i1-1] + eps && 
			f2 < min[i1] + eps) {
			/* see both */
			old = min[i1];
			min[i1] = f2;
			vp_udraw (x2, f2);
		    } else if (f1 > min[i1-1] && f2 > min[i1]) {
			/* see none */
			vp_umove (x2, f2);
		    } else if (f1 < old && f2 > min[i1]) {
			/* see left only */
			z = (f1 - old)/
			    (f1 - old + min[i1] - f2);
			x = (1.- z) * x1 + z * x2;
			f = (1.- z) * f1 + z * f2;
			vp_udraw (x, f);
			vp_umove (x2, f2);
		    } else if (f1 > min[i1-1] && f2 < min[i1]) {
			/* see right only */
			z = (f1 - min[i1-1]) / 
			    (f1 - min[i1-1] + min[i1] - f2);
			x = (1.- z) * x1 + z * x2;
			f = (1.- z) * f1 + z * f2;
			old = min[i1];
			min[i1] = f2;
			vp_umove (x, f);
			vp_udraw (x2, f2);
		    } else {
			/* see none */
			vp_umove (x2, f2);
		    }
		}
	    }

	    if (uflag) {
		vp_color(plotcolup);
		
		/* update the max */
		for (i1 = 0; i1 < nx; i1++) {
		    max[i1] = max[i1+1] - dz;
		}

		f2 = fff[i3][i2][0];
		x2 = 0.; 
		vp_umove (x2,f2);
			
		old = f2;
		for (i1 = 1; i1 < nx; i1++) {
		    f1 = f2;
		    x1 = x2;
		    x2 = dx * i1;
		    f2 = fff[i3][i2][i1];

		    if (f1 > max[i1-1] - eps && 
			f2 > max[i1] - eps) {
			/* see both */
			old = max[i1];
			max[i1] = f2;
			vp_udraw (x2, f2);
		    } else if (f1 < max[i1-1] && f2 < max[i1]) { 
			/* see none */
			vp_umove (x2, f2);
		    } else if (f1 > old && f2 < max[i1]) {
			/* see left */
			z = (f1 - old)/
			    (f1 - old + max[i1] - f2);
			x = (1.- z) * x1 + z * x2;
			f = (1.- z) * f1 + z * f2;
			vp_udraw (x, f);
			vp_umove (x2, f2);
		    } else if (f1 < max[i1-1] && f2 > max[i1]) {
                        /* see right */
			z = (f1-max[i1-1])/
			    (f1-max[i1-1] + max[i1] - f2);
			x = (1.- z) * x1 + z * x2;
			f = (1.- z) * f1 + z * f2;
			old = max[i1];
			max[i1] = f2;
			vp_umove (x, f);
			vp_udraw (x2, f2);
		    } else {
			/* see none */
			vp_umove (x2, f2);
		    }
		}
	    } /* uflag */

	} /* i2 */


	if (axis) {
	    vp_fat (axisfat);
	    vp_color (VP_WHITE);
	    r = 0.05 * fabsf (ay[1] - ay[0]);

	    if (axis1) {
		vp_arrow(ax[0],ay[0],
			 ax[0]+1.2*(ax[2]-ax[0]),
			 ay[0]+1.2*(ay[2]-ay[0]),r);
		if (NULL != labels[0]) 
		    vp_text (ax[2] - 0.5, ay[2] - 0.5, axissz, 0, labels[0]);
	    }

	    if (axis2) {
		x2 = ax[0]+1.2*(ax[1]-ax[0]);
		y2 = ay[0]+1.2*(ay[1]-ay[0]);
		xs = ax[0] - .15 * (ax[1]-ax[0]);
		ys = ay[0] - .15 * (ay[1]-ay[0]);
		vp_arrow( xs, ys, x2, y2, r);
		if (NULL != labels[1]) 
		    vp_text (ax[1] - 0.5*cosa, ay[1] + 0.5*sina, axissz, alpha, labels[1]);
	    }

	    if (axis3) {
		vp_arrow (ax[0], ay[0],
			  ax[0]+1.2*(ax[3]-ax[0]),
			  ay[0]+1.2*(ay[3]-ay[0]), r);
		if (NULL != labels[2]) 
		    vp_text (ax[0], ay[0] - 0.1, axissz, 0, labels[2]);
	    }
	}

	if (NULL != title) {
	    vp_fat (titlefat);
	    vp_color (VP_WHITE);

	    vp_text (0.4, VP_STANDARD_HEIGHT-0.8, titlsz, 0, title);
	}
    } /* i3 */


    exit (0);
}
