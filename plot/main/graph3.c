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
    int n1,n2,n3, frame2,frame3, i1,i2,i3, iframe, np=3;
    int n1pix,n2pix, m1pix,m2pix, n1front,n2front, movie, nframe=1, dframe; 
    float point1, point2, *front, **top, *x, *y, *side, o1, d1, o2, d2;    
    float min, max, f, frame1, oo1, dd1;
    bool nomin, nomax;
    char *label1, *label2, *unit1, *unit2;
    off_t esize;
    bool flat;
    sf_file in;

    sf_init(argc,argv);
    in = sf_input("in");
    vp_init();

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    n3 = sf_leftsize(in,2);

    esize = sf_esize(in);
    sf_unpipe(in,n1*n2*n3*esize);

    front = sf_floatalloc(n1);
    side = sf_floatalloc(n2);
    top = sf_floatalloc2(n1,n2);

    nomin = !sf_getfloat("min",&min);
    /* minimum function value */
    nomax = !sf_getfloat("max",&max);
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

    /* for proper frame */
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

    if (!sf_getint("dframe",&dframe)) dframe=1;
    /* frame increment in a movie */

    switch (movie) {
	case 0:
	    nframe = 1;
	    break;
	case 1:
	    sf_error("movie=1 is not implemented yet");
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

    vp_cubeplot_init (n1pix, n2pix, n1front, n2front, flat); 
    vp_frame_init (in,"blt",false);
    vp_plot_init(n3);

    for (iframe=0; iframe < nframe; iframe++) {
	for (i3=0; i3 < n3; i3++) {
	    vp_plot_set (i3);

	    if (0 == iframe || 3 == movie) {		
		vp_cubecoord(true,x[0],x[n1-1],min,max);

		sf_seek(in,(off_t) (i3*n1*n2+frame3*n1)*esize,SEEK_SET);
		sf_floatread(front,n1,in);
		
		vp_umove(x[0],front[0]);
		for (i1=1; i1 < n1; i1++) {
		    vp_udraw(x[i1],front[i1]);
		}
	    }
	    
	    if (0 == iframe || 2 == movie) {
		vp_cubecoord(false,y[0],y[n2-1],min,max);

		for (i2=0; i2 < n2; i2++) {
		    sf_seek(in,(off_t) (i3*n1*n2+i2*n1+frame2)*esize,SEEK_SET);
		    sf_floatread(&f,1,in);
		    if (i2==0) {
			vp_umove(y[0],f);
		    } else {
			vp_udraw(y[i2],f);
		    }
		}
	    }

	    vp_coordinates();

	    if (0 == iframe || 1 == movie) {
		sf_seek(in,(off_t) (i3*n1*n2*esize),SEEK_SET);
		sf_floatread(top[0],n1*n2,in);		
	    }

	    /*
	    if (flat) {
		for (i=0; i < n2front; i++) {
		    i2 = n2*i/(float) n2front;
		    for (j=n1front; j < n1pix; j++) {
			i3 = n3*(j-n1front)/(float) (n1pix-n1front);
			buf[i][j] = top[i2][i3];
		    }
		}
	    } 
	    */
	}

	/*
	if ((0 == iframe || 1 == movie || 2 == movie) && !flat) {
	    for (j=n1front; j < n1pix; j++) {
		i3 = n3*(j-n1front)/(float) (n1pix-n1front);
		i0 = (j-n1front)*(n2pix-n2front)/(float) (n1pix-n1front);
		for (i=i0; i < i0+n2front; i++) {
		    i2 = n2*(i-i0)/(float) n2front;
		    if (i2 < n2)
			buf[i][j] = top[i2][i3];
		}
	    }
	}
	*/

	if (iframe > 0) vp_erase (); 
	
	vp_plot_unset();
	vp_cubeframe((frame1-oo1)/dd1,frame2,frame3);
	
	switch (movie) {
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
    
    sf_close();
    exit(0);
}

/* 	$Id: grey3.c 1055 2005-03-16 21:18:43Z savap $	 */
