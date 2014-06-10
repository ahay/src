/* Generate 3-D cube plot for surfaces.

Takes: > plot.vpl
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

int main(int argc, char* argv[])
{
    int n1,n2,n3, frame2,frame3, i1,i2,i3, iframe, np=3, orient;
    int n1pix,n2pix, m1pix,m2pix, n1front,n2front, movie, nframe=1; 
    float point1, point2, **front, **side, **top, **topt, *x, *y, o1, d1, o2, d2;    
    float min, max, f, frame1, dframe, oo1, dd1;
    bool nomin, nomax, yreverse;
    char *label1, *label2, *unit1, *unit2;
    off_t esize;
    bool flat;
    vp_contour cnt;
    sf_file in=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    vp_init();

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    n3 = sf_leftsize(in,2);

    esize = sf_esize(in);
    sf_unpipe(in,n1*n2*n3*esize);

    if (!sf_getint("orient",&orient)) orient=1;
    /* function orientation */

    if (!sf_getbool("yreverse",&yreverse)) yreverse=false;

    front = sf_floatalloc2(n1,1);
    side = sf_floatalloc2(1,n2);
    top = sf_floatalloc2(n1,n2);
    topt = (2==orient)? sf_floatalloc2(n2,n1): NULL;

    nomin = (bool) !sf_getfloat("min",&min);
    /* minimum function value */
    nomax = (bool) !sf_getfloat("max",&max);
    /* maximum function value */

    if (nomin) min = +FLT_MAX;
    if (nomax) max = -FLT_MAX;
    if (nomin || nomax) {
	for (i3=0; i3 < n3; i3++) {
	    sf_floatread(top[0],n1*n2,in);
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    f = top[i2][i1];
		    if (nomax && f > max) max=f;
		    if (nomin && f < min) min=f;
		}
	    }
	}
    }


    if (min == 0. && max == 0.) {
	min = -1.;
	max = 1.;
    } else if (min == max) {
        max *= 1.04;
        min *= 0.96;
    }

    if (!sf_histfloat(in,"o1",&o1)) o1=0.;
    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_histfloat(in,"o2",&o2)) o2=0.;
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;

    switch(orient) {
	case 1:
	default:    
	    dd1 = (min-max)/np;
	    oo1 = max+0.5*dd1;

	    sf_putint(in,"n1",np);
	    sf_putfloat(in,"o1",oo1);
	    sf_putfloat(in,"d1",dd1);

	    sf_putint(in,"n2",n1);
	    sf_putfloat(in,"o2",o1);
	    sf_putfloat(in,"d2",d1);

	    sf_putint(in,"n3",n2);
	    sf_putfloat(in,"o3",o2);
	    sf_putfloat(in,"d3",d2);

	    break;
	case 2:
	    dd1 = (max-min)/np;
	    oo1 = min+0.5*dd1;

	    sf_putint(in,"n2",np);
	    sf_putfloat(in,"o2",oo1);
	    sf_putfloat(in,"d2",dd1);

	    sf_putint(in,"n3",n2);
	    sf_putfloat(in,"o3",o2);
	    sf_putfloat(in,"d3",d2);

	    break;
	case 3:
	    dd1 = (max-min)/np;
	    oo1 = min+0.5*dd1;

	    sf_putint(in,"n1",n2);
	    sf_putfloat(in,"o1",o2);
	    sf_putfloat(in,"d1",d2);

	    sf_putint(in,"n2",n1);
	    sf_putfloat(in,"o2",o1);
	    sf_putfloat(in,"d2",d1);

	    sf_putint(in,"n3",np);
	    sf_putfloat(in,"o3",oo1);
	    sf_putfloat(in,"d3",dd1);

	    break;
    }

    label1 = sf_histstring(in,"label1");
    label2 = sf_histstring(in,"label2");

    if (NULL != label1) {
	sf_putstring(in,"label1","");
	sf_putstring(in,"label2",label1);
	free(label1);
    }
    if (NULL != label2) {
	sf_putstring(in,"label3",label2);
	free(label2);
    }

    unit1 = sf_histstring(in,"unit1");
    unit2 = sf_histstring(in,"unit2");

    if (NULL != unit1) {
	sf_putstring(in,"unit1","");
	sf_putstring(in,"unit2",unit1);
	free(unit1);
    }
    if (NULL != unit2) {
	sf_putstring(in,"unit3",unit2);
	free(unit2);
    }

    if (!sf_getfloat("point1",&point1)) point1=0.5;
    /* fraction of the vertical axis for front face */
    if (!sf_getfloat("point2",&point2)) point2=0.5;
    /* fraction of the horizontal axis for front face */

    if (!sf_getfloat("frame1",&frame1)) frame1=0.5*(min+max);
    if (!sf_getint("frame2",&frame2)) frame2=n1-1;
    if (!sf_getint("frame3",&frame3)) frame3=0;
    /* frame numbers for cube faces */

    /* sanity check */
    if (min < max) {
	if (frame1 < min) frame1 = min;
	if (frame1 > max) frame1 = max;
    } else {
	if (frame1 > min) frame1 = min;
	if (frame1 < max) frame1 = max;
    }
    if (frame2 < 0) frame2 = 0;
    if (frame3 < 0) frame3 = 0;
    if (frame2 >= n1) frame2 = n1-1;
    if (frame3 >= n2) frame3 = n2-1;

    if (!sf_getint("movie",&movie)) movie=0;
    /* 0: no movie, 1: movie over axis 1, 2: axis 2, 3: axis 3 */

    if (!sf_getfloat("dframe",&dframe)) dframe=1;
    /* frame increment in a movie */

    switch (movie) {
	case 0:
	    nframe = 1;
	    break;
	case 1:
	    nframe = 0.5+(max-frame1)/dframe;
	    break;
	case 2:
	    nframe = (n1-frame2)/dframe;
	    break;
	case 3:
	    nframe = (n2-frame3)/dframe;
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

    x = sf_floatalloc(n1);
    y = sf_floatalloc(n2);

    for (i1=0; i1 < n1; i1++) {
	x[i1] = o1 + i1*d1;
    }

    for (i2=0; i2 < n2; i2++) {
	y[i2] = o2 + i2*d2;
    }

    if (!sf_getbool("flat",&flat)) flat=true;
    /* if n, display perspective view */

    vp_cubeplot_init (n1pix, n2pix, n1front, n2front, flat, false); 
    vp_frame_init (in,"blt",false);
    vp_plot_init(n3);
    cnt = vp_contour_init(false,
			  n1,o1,d1,0.,
			  n2,o2,d2,0.); 

    for (iframe=0; iframe < nframe; iframe++) {
	if (iframe > 0) vp_erase (); 

	for (i3=0; i3 < n3; i3++) {
	    vp_plot_set (i3);

	    /* front face */

	    switch(orient) {
		case 1:
		default:
		    sf_seek(in,(off_t) (i3*n1*n2+frame3*n1)*esize,SEEK_SET);
		    sf_floatread(front[0],n1,in);
		    
		    vp_cubecoord(3,x[0],x[n1-1],min,max);
		    vp_umove(x[0],front[0][0]);
		    for (i1=1; i1 < n1; i1++) {
			vp_udraw(x[i1],front[0][i1]);
		    }

		    break;
		case 2:
		    sf_seek(in,(off_t) (i3*n1*n2+frame3*n1)*esize,SEEK_SET);
		    sf_floatread(front[0],n1,in);

		    if (yreverse) {
			vp_cubecoord(3,min,max,x[n1-1],x[0]);
			vp_umove(front[0][n1-1],x[n1-1]);
			for (i1=n1-2; i1 >= 0; i1--) {
			    vp_udraw(front[0][i1],x[i1]);
			}
		    } else {
			vp_cubecoord(3,min,max,x[0],x[n1-1]);
			vp_umove(front[0][0],x[0]);
			for (i1=1; i1 < n1; i1++) {
			    vp_udraw(front[0][i1],x[i1]);
			}
		    }
		    
		    break;
		case 3:
		    sf_seek(in,(off_t) (i3*n1*n2*esize),SEEK_SET);
		    sf_floatread(top[0],n1*n2,in);

		    vp_cubecoord(3,x[0],x[n1-1],y[0],y[n2-1]);
		    vp_contour_draw(cnt,false,top,frame1);

		    break;
	    }
		    

	    /* side face */

	    switch(orient) {
		case 1:
		default:	    
		    for (i2=0; i2 < n2; i2++) {
			sf_seek(in,(off_t) (i3*n1*n2+i2*n1+frame2)*esize,SEEK_SET);
			sf_floatread(side[i2],1,in);
		    }
		    
		    vp_cubecoord(2,y[0],y[n2-1],min,max);
		    vp_umove(y[0],side[0][0]);
		    for (i2=1; i2 < n2; i2++) {
			vp_udraw(y[i2],side[i2][0]);
		    }

		    break;
		case 2:
		    sf_seek(in,(off_t) (i3*n1*n2*esize),SEEK_SET);
		    sf_floatread(top[0],n1*n2,in);

		    /* transpose */
		    for (i2=0; i2 < n2; i2++) {
			for (i1=0; i1 < n1; i1++) {
			    topt[i1][i2] = top[i2][i1];
			}
		    }
		    
		    if (yreverse) {
			vp_cubecoord(2,y[0],y[n2-1],x[n1-1],x[0]);
		    } else {
			vp_cubecoord(2,y[0],y[n2-1],x[0],x[n1-1]);
		    }
		    vp_contour_draw(cnt,false,topt,frame1);
		    
		    break;
		case 3:
		    for (i2=0; i2 < n2; i2++) {
			sf_seek(in,(off_t) (i3*n1*n2+i2*n1+frame2)*esize,SEEK_SET);
			sf_floatread(side[i2],1,in);
		    }
		    
		    vp_cubecoord(2,min,max,y[0],y[n2-1]);
		    vp_umove(side[0][0],y[0]);
		    for (i2=1; i2 < n2; i2++) {
			vp_udraw(side[i2][0],y[i2]);
		    }

		    break;
	    }

	    /* top face */

	    switch(orient) {
		case 1:
		default:
		    sf_seek(in,(off_t) (i3*n1*n2*esize),SEEK_SET);
		    sf_floatread(top[0],n1*n2,in);
		    
		    vp_cubecoord(1,x[0],x[n1-1],y[0],y[n2-1]);
		    vp_contour_draw(cnt,false,top,frame1);

		    break;
		case 2:
		    for (i2=0; i2 < n2; i2++) {
			sf_seek(in,(off_t) (i3*n1*n2+i2*n1+frame2)*esize,SEEK_SET);
			sf_floatread(side[i2],1,in);
		    }
		    
		    vp_cubecoord(1,min,max,y[0],y[n2-1]);
		    vp_umove(side[0][0],y[0]);
		    for (i2=1; i2 < n2; i2++) {
			vp_udraw(side[i2][0],y[i2]);
		    }
		    
		    break;
		case 3:
		    sf_seek(in,(off_t) (i3*n1*n2+frame3*n1)*esize,SEEK_SET);
		    sf_floatread(front[0],n1,in);
		    
		    vp_cubecoord(1,x[0],x[n1-1],min,max);
		    vp_umove(x[0],front[0][0]);
		    for (i1=1; i1 < n1; i1++) {
			vp_udraw(x[i1],front[0][i1]);
		    }

		    break;
	    }
	}
	
	vp_plot_unset();
	vp_coordinates();

	switch(orient) {
	    case 1:
	    default:
		vp_cubeframe((frame1-oo1)/dd1,frame2,frame3);
		break;
	    case 2:
		vp_cubeframe(frame2,(frame1-oo1)/dd1,frame3);
		break;
	    case 3:
		vp_cubeframe(frame2,frame3,(frame1-oo1)/dd1);
		break;
	}
	
	switch (movie) {
	    case 2:
		frame2 += (int) dframe;
		break;
	    case 3:
		frame3 += (int) dframe;
		break;
	    case 1:
		frame1 += dframe;
		break;
	    default:
		break;
	}

	vp_purge(); 
    } /* frame loop */


    exit(0);
}
