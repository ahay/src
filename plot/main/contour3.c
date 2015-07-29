/* Generate 3-D contour plot.

Takes: > plot.vpl
*/
/*
  Copyright (C) 2008 University of Texas at Austin

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

int main(int argc, char* argv[])
{
    bool flat, scalebar, nomin=true, nomax=true, barreverse, hasc, hasdc, hasc0, yreverse;
    char *cfilename;
    int n1,n2,n3, frame1,frame2,frame3, i1,i2,i3, iframe, nc, nc0, ic, maxstr;
    int n1pix, n2pix, m1pix, m2pix, n1front,n2front, movie, nframe=1, dframe; 
    float point1, point2, barmin, barmax, o1,o2,o3, d1,d2,d3, *c, dc, c0, zi;
    float min1,min2,min3,max1,max2,max3, zmin, zmax;
    float **front, **top, **side;
    vp_contour cfront, ctop, cside;
    sf_file in=NULL, cfile=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    vp_init();

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    n3 = sf_leftsize(in,2);

    if (!sf_histfloat(in,"o1",&o1)) o1=0.;
    if (!sf_histfloat(in,"o2",&o2)) o2=0.;  
    if (!sf_histfloat(in,"o3",&o3)) o3=0.;
 
    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;
    if (!sf_histfloat(in,"d3",&d3)) d3=1.;

    if (!sf_getfloat("min1",&min1)) min1=o1;
    if (!sf_getfloat("min2",&min2)) min2=o2;
    if (!sf_getfloat("min3",&min3)) min3=o3;

    if (!sf_getfloat("max1",&max1)) max1=o1+(n1-1)*d1;
    if (!sf_getfloat("max2",&max2)) max2=o2+(n2-1)*d2;
    if (!sf_getfloat("max3",&max3)) max3=o3+(n3-1)*d3;
    /* data window to plot */

    if (!sf_getbool ("yreverse",&yreverse)) yreverse=true;
    /* if y, reverse the first axis */
    if (yreverse) {
	zi = min1;
	min1 = max1;
	max1 = zi;
    }

    hasc = (bool) (NULL != (cfilename = sf_getstring("cfile")));
    /* contours in a file */

    if (hasc) {
	cfile = sf_input(cfilename);
	nc = sf_filesize(cfile);
	nc0 = nc;
    } else {
	if (!sf_getint("nc",&nc0)) nc0=50;
	/* number of contours */
	nc=nc0;
    }

    c = sf_floatalloc(nc);

    if (hasc) {
	sf_floatread(c,nc,cfile);
	sf_fileclose(cfile);

	hasdc = false;
	hasc0 = false;
    } else {
	hasc = sf_getfloats("c",c,nc);

	hasdc = sf_getfloat("dc",&dc);
	/* contour increment */
	hasc0 = sf_getfloat("c0",&c0);
	/* first contour */
    }

    if (!sf_getfloat("point1",&point1)) point1=0.5;
    /* fraction of the vertical axis for front face */
    if (!sf_getfloat("point2",&point2)) point2=0.5;
    /* fraction of the horizontal axis for front face */

    if (!sf_getint("frame1",&frame1)) frame1=0;
    if (!sf_getint("frame2",&frame2)) frame2=n2-1;
    if (!sf_getint("frame3",&frame3)) frame3=0;
    /* frame numbers for cube faces */

    /* sanity check */
    if (frame1 < 0) frame1 = 0;
    if (frame2 < 0) frame2 = 0;
    if (frame3 < 0) frame3 = 0;
    if (frame1 >= n1) frame1 = n1-1;
    if (frame2 >= n2) frame2 = n2-1;
    if (frame3 >= n3) frame3 = n3-1;

    if (!sf_getint("movie",&movie)) movie=0;
    /* 0: no movie, 1: movie over axis 1, 2: axis 2, 3: axis 3 */

    if (!sf_getint("dframe",&dframe)) dframe=1;
    /* frame increment in a movie */

    switch (movie) {
	case 0:
	    nframe = 1;
	    break;
	case 1:
	    nframe = (n1-frame1)/dframe;
	    break;
	case 2:
	    nframe = (n2-frame2)/dframe;
	    break;
	case 3:
	    nframe = (n3-frame3)/dframe;
	    break;
	default:
	    sf_error("movie=%d is outside [0,3]",movie);
	    break;
    }

    if (sf_getint("nframe",&iframe) && iframe < nframe) nframe=iframe;
    /* number of frames in a movie */

    if (nframe < 1) nframe=1;

    if (!sf_getint("n1pix",&n1pix)) n1pix = n1/point1+n3/(1.-point1);
    /* number of vertical pixels */
    if (!sf_getint("n2pix",&n2pix)) n2pix = n2/point2+n3/(1.-point2);
    /* number of horizontal pixels */

    m1pix = VP_STANDARD_HEIGHT*72;
    m2pix = VP_STANDARD_HEIGHT*72/VP_SCREEN_RATIO;

    if (n1pix < m1pix) n1pix=m1pix;
    if (n2pix < m2pix) n2pix=m2pix;

    n1front = n1pix*point1;
    n2front = n2pix*point2;

    if (n1front <= 0) n1front=1;
    if (n2front <= 0) n2front=1;
    if (n1front >= n1pix) n1front=n1pix-1;
    if (n2front >= n2pix) n2front=n2pix-1;

    front = sf_floatalloc2(n1,n2);
    top = sf_floatalloc2(n2,n3);
    side = sf_floatalloc2(n1,n3);

    if (!sf_getbool("flat",&flat)) flat=true;
    /* if n, display perspective view */

    if (!sf_getbool ("wantscalebar",&scalebar) && 
	!sf_getbool ("scalebar",&scalebar)) scalebar = false;
    /* if y, draw scalebar */
    if (scalebar) {
	nomin = (bool) (!sf_getfloat("minval",&barmin));
	/* minimum value for scalebar (default is the data minimum) */
	nomax = (bool) (!sf_getfloat("maxval",&barmax));
	/* maximum value for scalebar (default is the data maximum) */

	if (!sf_getbool("barreverse",&barreverse)) barreverse=false;
	/* if y, go from small to large on the bar scale */
    }

    sf_unpipe(in,n1*n2*n3*sizeof(float));

    if (!hasc) {
	if (!hasdc || !hasc0) {
	    zmin=SF_HUGE;
	    zmax=-SF_HUGE;

	    for (i3=0; i3 < n3; i3++) {
		sf_floatread(front[0],n1*n2,in);

		for (i2=0; i2 < n2; i2++) {
		    for (i1=0; i1 < n1; i1++) {
			zi= front[i2][i1];
			if      (zi < zmin) zmin=zi;
			else if (zi > zmax) zmax=zi;
		    }
		}
	    }
	    
	    if (hasdc) {
		for (c0 = floorf(zmin/dc) * dc - dc; c0 < zmin; c0 += dc) ;
	    } else if (hasc0) {		
		nc = vp_optimal_scale(nc0, true, false, "%g",
				      zmin-c0, zmax-c0,   &zi, &dc, &maxstr);
	    } else {
		nc = vp_optimal_scale(nc0, true, false, "%g",
				      zmin,    zmax,      &c0, &dc, &maxstr);
	    }
	}
	for (ic=0; ic < nc; ic++) {
	    c[ic] = c0 + dc*ic;
	}
    }

    vp_cubeplot_init (n1pix, n2pix, n1front, n2front, flat, false); 
    vp_frame_init (in,"blt",false);
    if (scalebar && !nomin && !nomax) vp_barframe_init (in,barmin,barmax);

    cfront = vp_contour_init(true,
			     n1,o1,d1,0.,
			     n2,o2,d2,0.);
    cside = vp_contour_init(true,
			    n1,o1,d1,0.,
			    n3,o3,d3,0.);
    ctop = vp_contour_init(false,
			   n2,o2,d2,0.,
			   n3,o3,d3,0.);

    vp_plot_init(nc);

    for (iframe=0; iframe < nframe; iframe++) {
	if (iframe > 0) vp_erase (); 
	
	if (0 == iframe || 3 == movie) { 
	    sf_seek(in,(off_t) frame3*n1*n2*sizeof(float),SEEK_SET);
	    sf_floatread(front[0],n1*n2,in);
	}

	vp_cubecoord(3,min2,max2,min1,max1);

	for (ic=0; ic < nc; ic++) {
	    vp_plot_set (ic);
	    vp_contour_draw(cfront,false,front,c[ic]);
	}

	if (0 == iframe || 2 == movie) {
	    for (i3=0; i3 < n3; i3++) {
		sf_seek(in,(off_t) (i3*n1*n2+frame2*n1)*sizeof(float),SEEK_SET);
		sf_floatread(side[i3],n1,in);
	    }
	}

	vp_cubecoord(2,min3,max3,min1,max1);

	for (ic=0; ic < nc; ic++) {
	    vp_plot_set (ic);
	    vp_contour_draw(cside,false,side,c[ic]);
	}

	if (0 == iframe || 1 == movie) {
	    for (i3=0; i3 < n3; i3++) {
		for (i2=0; i2 < n2; i2++) {
		    sf_seek(in,(off_t) (i3*n1*n2+i2*n1+frame1)*sizeof(float),
			    SEEK_SET);
		    sf_floatread(&top[i3][i2],1,in);
		}
	    }
	}
	
	vp_cubecoord(1,min2,max2,min3,max3);
	
	for (ic=0; ic < nc; ic++) {
	    vp_plot_set (ic);
	    vp_contour_draw(ctop,false,top,c[ic]);
	}
	
	vp_plot_unset();
	vp_coordinates();
	vp_cubeframe(frame1,frame2,frame3);
	
	if (scalebar) {
	    if (barreverse) {
		vp_barframe_init (in,barmax,barmin);
	    } else {
		vp_barframe_init (in,barmin,barmax);
	    }
	}

	switch (movie) {
	    case 1:
		frame1 += dframe;
		break;
	    case 2:
		frame2 += dframe;
		break;
	    case 3:
		frame3 += dframe;
		break;
	    default:
		break;
	}

	vp_purge(); 
    } /* frame loop */


    exit(0);
}
