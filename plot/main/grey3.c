/* Generate 3-D cube plot.
   
   Takes: > plot.vpl

   Requires an "unsigned char" input (the output of sfbyte).

April 2012 program of the month:
http://ahay.org/rsflog/index.php?/archives/287-Program-of-the-month-sfgrey3.html
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
    int n1,n2,n3, nreserve, frame1,frame2,frame3, i1,i2,i3, i,j, i0,j0, iframe;
    int n1pix,n2pix, m1pix,m2pix, n1front,n2front, movie, nframe=1, dframe; 
    float point1, point2, barmin, barmax, minmax[2];
    bool flat, scalebar, nomin=true, nomax=true, barreverse;
    char *color, *barfile;
    unsigned char **front, **top, **side, **buf, b, *barbuf[1];
    sf_file in, bar=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    vp_init();

    if (SF_UCHAR != sf_gettype(in)) sf_error("Need uchar input");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    n3 = sf_leftsize(in,2);

    if (!sf_getfloat("point1",&point1)) point1=0.5;
    /* fraction of the vertical axis for front face */
    if (!sf_getfloat("point2",&point2)) point2=0.5;
    /* fraction of the horizontal axis for front face */
    
    if (!sf_getint("frame1",&frame1)) frame1=0;
    /* top frame number */
    if (!sf_getint("frame2",&frame2)) frame2=n2-1;
    /* side frame number */
    if (!sf_getint("frame3",&frame3)) frame3=0;
    /* front frame number */

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

    front = sf_ucharalloc2(n1,n2);
    top   = sf_ucharalloc2(n3,n2);
    side  = sf_ucharalloc2(n1,n3);
    buf   = sf_ucharalloc2(n1pix,n2pix);

    if (!sf_getbool("flat",&flat)) flat=true;
    /* if n, display perspective view */

    if (!sf_getbool ("wantscalebar",&scalebar) && 
	!sf_getbool ("scalebar",&scalebar)) scalebar = false;
    /* if y, draw scalebar */
    if (scalebar) {
	nomin = (bool) !sf_getfloat("minval",&barmin);
	/* minimum value for scalebar (default is the data minimum) */
	nomax = (bool) !sf_getfloat("maxval",&barmax);
	/* maximum value for scalebar (default is the data maximum) */
	barbuf[0] = (unsigned char*) sf_alloc(VP_BSIZE,sizeof(unsigned char));
	if (!sf_getbool("barreverse",&barreverse)) barreverse=false;
	/* if y, go from small to large on the bar scale */

	barfile = sf_getstring("bar");
	/* file for scalebar data */

	if (NULL == barfile) {
	    barfile=sf_histstring(in,"bar");
	    if (NULL == barfile) sf_error("Need bar=");
	}

	bar = sf_input(barfile);
	if (SF_UCHAR != sf_gettype(bar)) sf_error("Need uchar in bar");
	
	if (nomin) nomin = (bool) !sf_histfloat(bar,"minval",&barmin);
	/* minimum value for scalebar */
	if (nomax) nomax = (bool) !sf_histfloat(bar,"maxval",&barmax);
	/* maximum value for scalebar */
    }

    sf_unpipe(in,(off_t) n1*n2*n3);

    vp_cubeplot_init (n1pix, n2pix, n1front, n2front, flat, false); 
    vp_frame_init (in,"blt",false);
    if (scalebar && !nomin && !nomax) vp_barframe_init (in,barmin,barmax);

    /* initialize color table */
    if (NULL == (color = sf_getstring("color"))) color="i";
    /* color scheme (default is i) */
    if (!sf_getint ("nreserve",&nreserve)) nreserve = 8;
    /* reserved colors */
    vp_rascoltab (nreserve, color);

    /* fill empty areas */
    b = '\0';

    if (flat) {
	/* rectangle */
	for (i=n2front; i < n2pix; i++) {
	    for (j=n1front; j < n1pix; j++) {
		buf[i][j] = b;
	    }
	}
    } else {
	/* two triangles */
	for (i=n2front; i < n2pix; i++) {
	    j0 = (i-n2front)*(n1pix-n1front)/(float) (n2pix-n2front);
	    for (j=0; j < j0; j++) {
		buf[i][j] = b;
	    }
	}
	for (j=n1front; j < n1pix; j++) {
	    i0 = (j-n1front)*(n2pix-n2front)/(float) (n1pix-n1front);
	    for (i=0; i < i0; i++) {
		buf[i][j] = b;
	    }
	}
    }

    for (iframe=0; iframe < nframe; iframe++) {
	/* read data and fill display buffer */

	if (0 == iframe || 3 == movie) { 
	    sf_seek(in,(off_t) frame3*n1*n2,SEEK_SET);
	    sf_ucharread(front[0],n1*n2,in);
	    
	    for (i=0; i < n2front; i++) {
		i2 = n2*i/(float) n2front;
		for (j=0; j < n1front; j++) {
		    i1 = n1*(n1front-j-0.5)/(float) n1front;
		    buf[i][j] = front[i2][i1];
		}
	    }
	    
	    if (scalebar) {
		sf_seek(bar,(off_t) frame3*(VP_BSIZE+2*sizeof(float)),
			SEEK_SET);

		sf_floatread(minmax,2,bar);
		sf_ucharread(barbuf[0],VP_BSIZE,bar);

		if (nomin) barmin=minmax[0];
		if (nomax) barmax=minmax[1];
	    }
	}

	if (0 == iframe || 2 == movie) {
	    for (i3=0; i3 < n3; i3++) {
		sf_seek(in,(off_t) i3*n1*n2+frame2*n1,SEEK_SET);
		sf_ucharread(side[i3],n1,in);
	    }

	    for (i=n2front; i < n2pix; i++) {
		i3 = n3*(i-n2front)/(float) (n2pix-n2front);
		if (flat) {
		    for (j=0; j < n1front; j++) {
			i1 = n1*(n1front-j-0.5)/(float) n1front;
			buf[i][j] = side[i3][i1];
		    }
		} else {
		    j0 = (i-n2front)*(n1pix-n1front)/
			(float) (n2pix-n2front);
		    for (j=j0; j < n1pix; j++) {
			i1 = n1*(n1front+j0-j-0.5)/(float) n1front;
			if (i1 >= 0)
			    buf[i][j] = side[i3][i1];
		    }
		}
	    }
	}

	if (0 == iframe || 1 == movie) {
	    for (i3=0; i3 < n3; i3++) {
		for (i2=0; i2 < n2; i2++) {
		    sf_seek(in,(off_t) i3*n1*n2+i2*n1+frame1,SEEK_SET);
		    sf_ucharread(&top[i2][i3],1,in);
		}
	    }

	    if (flat) {
		for (i=0; i < n2front; i++) {
		    i2 = n2*i/(float) n2front;
		    for (j=n1front; j < n1pix; j++) {
			i3 = n3*(j-n1front)/(float) (n1pix-n1front);
			buf[i][j] = top[i2][i3];
		    }
		}
	    } 
	}

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

	if (iframe > 0) vp_erase (); 

	vp_cuberaster(n1pix,n2pix,buf,frame1,frame2,frame3);
	
	if (scalebar) {
	    if (barreverse) {
		vp_barframe_init (in,barmax,barmin);
	    } else {
		vp_barframe_init (in,barmin,barmax);
	    }
	    vp_barraster(VP_BSIZE, barbuf);
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


    exit (0);
}

/* 	$Id$	 */
