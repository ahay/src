/* Contour plot.

Takes: > plot.vpl

Run "sfdoc stdplot" for more parameters.

December 2011 program of the month:
http://ahay.org/rsflog/index.php?/archives/277-Programs-of-the-month-sfcontour.html
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

#include <math.h>
#include <rsf.h>
#include <rsfplot.h>

int main (int argc, char* argv[])
{
    int n1, n2, n3, i3, nc0, nc, ic, n12, i1, i2, maxstr;
    char *cfilename;
    float **z, zi, dc, c0, zmin=0., zmax=0., *c;
    float  min1, min2, max1, max2, bmin, bmax, o1, o2, d1, d2;
    bool hasc, hasdc, hasc0, scalebar, nomin=false, nomax=false, pos, transp;
    vp_contour cnt;
    sf_file in=NULL, cfile=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    vp_init();

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    n3 = sf_leftsize(in,2);
    n12 = n1*n2;

    if (!sf_histfloat(in,"o1",&o1)) o1=0.;
    if (!sf_histfloat(in,"o2",&o2)) o2=0.;
    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;

    if (!sf_getfloat("min1",&min1)) min1=o1;           /* minimum on 1st axis */
    if (!sf_getfloat("min2",&min2)) min2=o2;           /* minimum on 2nd axis */
    if (!sf_getfloat("max1",&max1)) max1=o1+(n1-1)*d1; /* maximum on 1st axis */
    if (!sf_getfloat("max2",&max2)) max2=o2+(n2-1)*d2; /* maximum on 2nd axis */

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
    vp_plot_init(nc);

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

    if (!sf_getbool ("transp",&transp)) transp=true;
    /* if y, transpose the axes */

    z = sf_floatalloc2(n1,n2);

    if ((!sf_getbool ("wantscalebar",&scalebar) &&
	 !sf_getbool ("scalebar",&scalebar)) ||
	NULL == sf_getstring("barlabel")) scalebar = false;
    /* scale bar label */
    if (scalebar) {
	nomin = (bool) (!sf_getfloat("minval",&bmin));
	/* minimum value for scalebar (default is the data minimum) */
	nomax = (bool) (!sf_getfloat("maxval",&bmax));
	/* maximum value for scalebar (default is the data maximum) */
    }

    if (!sf_getbool("allpos",&pos)) pos=true;
    /* contour positive values only */

    cnt = vp_contour_init(transp,
			  n1,o1,d1,0.,
			  n2,o2,d2,0.);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(z[0],n12,in);
	
	if (!hasc) {
	    if (!hasdc || !hasc0) {
		zmin = z[0][0];
		zmax = z[0][0];
		for (i2=0; i2 < n2; i2++) {
		    for (i1=0; i1 < n1; i1++) {
			zi= z[i2][i1];
			if      (zi < zmin) zmin=zi;
			else if (zi > zmax) zmax=zi;
		    }
		}
		if (hasdc) {
		    for (c0 = floorf(zmin/dc) * dc - dc; c0 < zmin; c0 += dc) ;
		} else if (hasc0) {		
		    nc = vp_optimal_scale(nc0, true, false, "%g",
					  zmin-c0, zmax-c0, 
					  &zi, &dc, &maxstr);
		} else {
		    nc = vp_optimal_scale(nc0, true, false, "%g",
					  zmin,    zmax,    
					  &c0, &dc, &maxstr);
		}
	    }
	    for (ic=0; ic < nc; ic++) {
		c[ic] = c0 + dc*ic;
	    }
	}

	vp_stdplot_init (min1, max1, min2, max2,
			 transp,false,true,false);
	vp_frame_init(in,"tlb",false);

	if (i3 > 0) vp_erase();
	vp_frame();

	for (ic = 0; ic < nc; ic++) {
	    vp_plot_set (ic);
	    vp_contour_draw (cnt, pos, z, c[ic]);
	} 

	if (scalebar) {
	    if (nomin || nomax) {
		if (!hasc && (!hasdc || !hasc0)) {
		    if (nomin) bmin=zmin;
		    if (nomax) bmax=zmax;
		} else {
		    if (nomin) bmin = z[0][0];
		    if (nomax) bmax = z[0][0];
		    for (i2=0; i2 < n2; i2++) {
			for (i1=0; i1 < n1; i1++) {
			    zi= z[i2][i1];
			    if      (nomin && zi < bmin) bmin=zi;
			    else if (nomax && zi > bmax) bmax=zi;
			}
		    }
		}
	    }

	    vp_barframe_init (in,bmin,bmax);
	    vp_barline(nc,c,bmin,bmax);
	}

    } /* i3 */


    exit(0);
}
