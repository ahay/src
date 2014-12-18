/* Plot data with wiggly traces. 

Takes: > plot.vpl

Run "sfdoc stdplot" for more parameters.

June 2013 program of the month:
http://www.ahay.org/rsflog/index.php?/archives/340-Program-of-the-month-sfwiggle.html
*/
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

#include <rsf.h>
#include <rsfplot.h>

static void checkit(float *x, float *y);
static bool transp;

int main(int argc, char* argv[])
{
    int n1, n2, n3, i3, i2, i1, i, fatp, xmask, ymask;
    float o1, o2, o3, d1, d2, d3, **xp, pclip, zplot, zdata, x1, x2, y1, y2;
    float gpow, scale, x0, y0, zero, xmax, xmin, bias=0.0;
    float **pdata=NULL, *q=NULL, *x=NULL, *y=NULL, *px=NULL, *py=NULL, *tmp;
    bool poly, polyneg, seemean, verb, xreverse, yreverse;
    sf_file in, xpos;

    sf_init(argc,argv);
    in = sf_input("in");

    vp_init();

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    n3 = sf_leftsize(in,2);

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_histfloat(in,"o1",&o1)) o1=0.;

    if (NULL != sf_getstring("xpos")) {
	/* optional header file with trace positions */

	xpos = sf_input("xpos");	
	if (SF_FLOAT != sf_gettype(xpos)) sf_error("Need float xpos");

	xp = sf_floatalloc2(n2,n3);

	sf_floatread(xp[0],n2*n3,xpos);
	sf_fileclose(xpos);	
    } else {
	xpos = NULL;

	if (!sf_histfloat(in,"d2",&d2)) d2=1.;
	if (!sf_histfloat(in,"o2",&o2)) o2=0.;

	xp = sf_floatalloc2(n2,1);

	for (i2=0; i2 < n2; i2++) {
	    xp[0][i2] = o2+i2*d2;
	}
    }

    if (!sf_getfloat("xmax",&xmax)) {
	/* maximum trace position (if using xpos) */
	if (NULL == xpos) {
	    xmax=xp[0][n2-1];
	} else {
	    xmax=-FLT_MAX;
	    for (i3=0; i3 < n3; i3++) {
		for (i2=0; i2 < n2; i2++) {
		    if (xp[i3][i2] > xmax) xmax=xp[i3][i2];
		}
	    }
	}
    }

    if (!sf_getfloat("xmin",&xmin)) {
	/* minimum trace position (if using xpos) */
	if (NULL == xpos) {
	    xmin=xp[0][0];
	} else {
	    xmin=FLT_MAX;
	    for (i3=0; i3 < n3; i3++) {
		for (i2=0; i2 < n2; i2++) {
		    if (xp[i3][i2] < xmin) xmin=xp[i3][i2];
		}
	    }
	}
    }

    if (!sf_histfloat(in,"d3",&d3)) d3=1.;
    if (!sf_histfloat(in,"o3",&o3)) o3=0.;

    if (!sf_getbool("poly",&poly)) poly=false;
    /* if draw polygons */

    if (!sf_getbool("polyneg",&polyneg)) polyneg=false;
    /* if polygons for negative values */

    if (poly) {
	if (!sf_getint("fatp",&fatp)) fatp = 1;
	/* polygon border fatness */
	if (!sf_getint("xmask",&xmask)) xmask = 1;
	/* polygon filling */
	if (!sf_getint("ymask",&ymask)) ymask = 1;
	/* polygon filling */
    }

    if (!sf_getfloat("pclip",&pclip)) pclip = 98.;
    /* clip percentile */

    if (!sf_getfloat("zplot",&zplot)) zplot = 0.75;
    /* vertical separation */
    d2 = (n2>1)? (xmax-xmin)/(n2-1): 1.0;
    zplot *= d2;

    if (!sf_getfloat("clip",&zdata)) zdata = 0.;
    /* data clip (estimated from pclip by default */

    x1 = o1;
    x2 = o1+(n1-1)*d1;

    y1 = xmin-zplot;
    y2 = xmax+zplot;

    if (!sf_getbool("seemean",&seemean)) seemean = false;
    /* if y, plot mean lines of traces */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */
    if (!sf_getbool ("transp",&transp)) transp=false;
    /* if y, transpose the axes */
    if (!sf_getbool("yreverse",&yreverse)) yreverse=false;
    /* if y, reverse the vertical axis */
    if (!sf_getbool("xreverse",&xreverse)) xreverse=false;
    /* if y, reverse the horizontal axis */

    pdata = sf_floatalloc2 (n1,n2);
    x = sf_floatalloc (n1);
    y = sf_floatalloc (n1);

    if (poly) {
	px = sf_floatalloc (n1+2);
	py = sf_floatalloc (n1+2);
    }

    vp_stdplot_init (x1, x2, y1, y2, transp, false, yreverse, false);
    vp_frame_init(in,"blt",true);
    vp_plot_init(1);

    for (i3=0; i3 < n3; i3++) {
	if (zdata <= 0.) {
	    vp_gainpar (in,pdata,n1,n2,1,
			pclip,pclip,&zdata,&gpow,false,&bias,
			n3,0,n3);
	    if (verb) sf_warning("clip=%g",zdata);
	} else {	    
	  sf_floatread(pdata[0],n1*n2,in);
	}

	if (zdata > 0.) {
	  scale = zplot/zdata;
	} else {
	  scale = 0.;
	  for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
	      scale += fabsf(pdata[i2][i1]);
	    }
	  }
	}

	if (i3 > 0) vp_erase();

	vp_frame();
	vp_plot_set(0);

	/* draw traces */
	for (i2=0; i2 < n2; i2++) {
	    q = pdata[i2];
	    zero = (NULL==xpos)? xp[0][i2]: xp[i3][i2];

	    for (i1=0; i1 < n1; i1++) {
		x[i1] = o1 + i1*d1;
		y[i1] = zero + q[i1] * scale;
	    }

	    /* test? */

	    if (transp) {tmp=x; x=y; y=tmp;}

	    /* yreverse and xreverse */

	    vp_bgroup("wiggles");
	    vp_umove(x[0],y[0]);	    
	    for (i1=1; i1 < n1; i1++) {
		vp_udraw(x[i1],y[i1]);
	    }
	    vp_egroup();

	    if (seemean) {
		vp_bgroup("mean");

		x0 = x1;
		y0 = zero;
		checkit(&x0,&y0);

		vp_umove(x0,y0);

		x0 = x2;
		y0 = zero;
		checkit(&x0,&y0);

		vp_udraw(x0,y0);
		vp_egroup();
	    }

	    if (poly) {
		vp_bgroup("polygons");
		i = 0;
		if (( polyneg && q[0] < 0.) || 
		    (!polyneg && q[0] > 0.)) { /* start polygon */
		    x0 = x1;
		    y0 = zero;
		    checkit(&x0,&y0);

		    px[i] = x0;
		    py[i] = y0;
		    i++;
		    px[i] = x[0];
		    py[i] = y[0];
		    i++;
		}

		for (i1 = 0; i1 < n1; i1++) {
		    if (0 == i) {
			if (( polyneg && q[i1] < 0.) ||
			    (!polyneg && q[i1] > 0.)) { /* start polygon */
			    x0 = o1+(i1 - q[i1]/(q[i1] - q[i1-1]))*d1;
			    y0 = zero;
			    checkit(&x0,&y0);

			    px[i] = x0;
			    py[i] = y0;
			    i++;
			    px[i] = x[i1];
			    py[i] = y[i1];
			    i++;
			}
		    } else {
			if (( polyneg && q[i1] < 0.) ||
			    (!polyneg && q[i1] > 0.)) { /* continue polygon */
			    px[i] = x[i1];
			    py[i] = y[i1];
			    i++;
			} else { /* end polygon */
			    x0 = o1+(i1 - q[i1]/(q[i1] - q[i1-1]))*d1;
			    y0 = zero;
			    checkit(&x0,&y0);

			    px[i] = x0;
			    py[i] = y0;
			    i++;

			    vp_uarea (px, py, i, fatp, xmask, ymask);
			    i = 0;
			}
		    }
		} /* i1 */

		if (0 != i) { /* end polygon */
		    x0 = x2;
		    y0 = zero;
		    checkit(&x0,&y0);

		    px[i] = x0;
		    py[i] = y0;
		    i++;
		    vp_uarea (px, py, i, fatp, xmask, ymask);
		}
		vp_egroup();
	    } /* if poly */
	} /* i2 */

	vp_simpleframe();

    } /* i3 */

    exit(0);
}

static void checkit(float *x, float *y)
{
    float t;

    if (transp) { t=*x; *x=*y; *y=t; }
}
