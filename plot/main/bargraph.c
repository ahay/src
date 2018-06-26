/* Bar plot.

Takes: > plot.vpl

Run "sfdoc stdplot" for more parameters.
*/
/*
  Copyright (C) 2016 University of Texas at Austin
  
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


#include <string.h>
#include <float.h>
#include <math.h>

#include <rsf.h>
#include <rsfplot.h>

static int n;
static float *t, pclip;
static void getminmax(const float* f, float* min, float* max);

int main(int argc, char* argv[])
{
    bool transp, framenum;
    int n1, n2, n3, i1, i2, i3;
    float min1, max1, min2, max2, o3, d3, o1, d1, xi, yi, yi1, tt, wd, dx=0.0f;
    float **x, **y, **tmp, xp[4], yp[4];    
    float ***data=NULL;
    sf_datatype type;
    sf_file in; 

    sf_init(argc,argv);
    in = sf_input("in");
    vp_init();

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    n = n1*n2;
    n3 = sf_leftsize(in,2);
    if (n3 > 1) {
	if (!sf_histfloat(in,"o3",&o3)) o3=0.;
	if (!sf_histfloat(in,"d3",&d3)) d3=1.;
    }

    x = sf_floatalloc2(n1,n2);
    y = sf_floatalloc2(n1,n2);
    t = sf_floatalloc(n);

    if (!sf_getbool("wantframenum",&framenum)) framenum = (bool) (n3 > 1);
    /* if y, display third axis position in the corner */
    
    if (!sf_getfloat("pclip",&pclip)) pclip=100.; /* clip percentile */

    if (!sf_getfloat("width",&wd)) wd=0.8; /* bar width */
    wd *= 0.5;

    type = sf_gettype(in);
    switch (type) {
	case SF_FLOAT:
	    if (!sf_histfloat(in,"o1",&o1)) o1=0.;
	    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
	    
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    x[i2][i1] = o1 + i1*d1;
		}
	    }
	    dx = wd*d1;
	    break;
	case SF_COMPLEX:
	    data = sf_floatalloc3(2,n1,n2);
	    break;
	default:
	    sf_error("Wrong data type (need float or complex)");
    }

    vp_plot_init(n2);

    if (!sf_getbool ("transp",&transp)) transp=false;
    /* if y, transpose the axes */
 
    for (i3 = 0; i3 < n3; i3++) {
	if (SF_COMPLEX == type) {
	    sf_floatread(data[0][0],2*n,in);
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    x[i2][i1] = data[i2][i1][0];
		    y[i2][i1] = data[i2][i1][1];
		}
	    }
	    getminmax(x[0],&min1,&max1);
	    dx = wd*(max1-min1)/(n1-1);
	} else {
	    sf_floatread(y[0],n,in);
	    min1=o1;
	    max1=o1+(n1-1)*d1;
	}
	getminmax(y[0],&min2,&max2);
	
	vp_stdplot_init (min1, max1, min2, max2,
			 transp,false,false,true);
	vp_frame_init(in,"blt",false);

	if (transp) {
	    tmp=x; x=y; y=tmp;
	    tt=max1; max1=max2; max2=tt;
	    tt=min1; min1=min2; min2=tt;
	}

	if (i3 > 0) vp_erase();

	if (framenum) vp_framenum(o3+i3*d3);
	vp_frame();

	for (i2=0; i2 < n2; i2++) {
	    vp_plot_set (i2);

	    for (i1=0; i1 < n1; i1++) {
		xi = x[i2][i1];
		yi = y[i2][i1];
		yi1 = y[i2+1][i1];

		if (isfinite(xi) && 
		    isfinite(yi)) {

		    xp[0]=xi-dx;  yp[0]=0;
		    xp[1]=xi-dx;  yp[1]=yi;
		    xp[2]=xi+dx;  yp[2]=yi;
		    xp[3]=xi+dx;  yp[3]=0;

		    vp_ufill(xp,yp,4);
		} 
	    }
	}
	
	if (transp) {
	    tmp=x; x=y; y=tmp;
	    tt=max1; max1=max2; max2=tt;
	    tt=min1; min1=min2; min2=tt;
	}
    } 
   

    exit(0);
}

static void getminmax(const float* f, float* min, float* max)
{
    int i, m, nc;
    float fmin, fmax, fi, fbig, fsml, fdif;

    m=0;
    fmin=+FLT_MAX;
    fmax=-FLT_MAX;
    for (i=0; i < n; i++) {
	fi = f[i];
	if (isfinite(fi)) {
	    t[m] = fi;
	    m++;
	    if (fmin > fi) fmin=fi;
	    if (fmax < fi) fmax=fi;
	}
    }

    nc = (1.-0.01*pclip)*m;
    if (nc < 0) nc=0;
    if (nc >= m) nc=m-1;
    
    if (nc > 0 && nc < m-1) {
	fsml = sf_quantile(nc,m,t);
	fbig = sf_quantile(m-nc-1, m,t);
	fdif = fbig - fsml;

	*min = SF_MAX(fmin-fdif/8,fsml-fdif);
	*max = SF_MIN(fmax+fdif/8,fbig+fdif);
    } else {
	*min = fmin;
	*max = fmax;
    }
}


/* 	$Id$	 */
