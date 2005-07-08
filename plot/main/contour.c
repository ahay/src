/* Contour plot.

Takes: > plot.vpl

Run "sfdoc stdplot" for more parameters.
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

static void contour (bool pos, float **z, int n1, int n2, float c);
static float delta (float a, float b, float  c);
static bool cconnect (bool pos, bool **south, bool** west, float** z,
		      int n1, int n2, float  c, int *ix, int *iy);
static void draw (bool mask, float x, float y);

static bool transp;
static float o1, o2, d1, d2;

int main (int argc, char* argv[])
{
    int n1, n2, n3, i3, nc0, nc, ic, n12, i1, i2;
    float **z, zi, dc, c0, zmin=0., zmax=0., *c;
    float  min1, min2, max1, max2, bmin, bmax;
    bool hasc, hasdc, hasc0, scalebar, nomin=false, nomax=false, pos;
    sf_file in;
    
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

    if (!sf_getfloat("min1",&min1)) min1=o1;
    if (!sf_getfloat("min2",&min2)) min2=o2;
    if (!sf_getfloat("max1",&max1)) max1=o1+(n1-1)*d1;
    if (!sf_getfloat("max2",&max2)) max2=o2+(n2-1)*d2;
    /* data window to plot */

    if (!sf_getint("nc",&nc0)) nc0=50;
    /* number of contours */
    nc=nc0;

    c = sf_floatalloc(nc);
    vp_plot_init(nc);

    hasc = sf_getfloats("c",c,nc);
    hasdc = sf_getfloat("dc",&dc);
    /* contour increment */
    hasc0 = sf_getfloat("c0",&c0);
    /* first contour */

    if (!sf_getbool ("transp",&transp)) transp=true;
    /* if y, transpose the axes */

    z = sf_floatalloc2(n1,n2);

    if ((!sf_getbool ("wantscalebar",&scalebar) &&
	 !sf_getbool ("scalebar",&scalebar)) ||
	NULL == sf_getstring("barlabel")) scalebar = false;
    if (scalebar) {
	nomin = !sf_getfloat("minval",&bmin);
	/* minimum value for scalebar (default is the data minimum) */
	nomax = !sf_getfloat("maxval",&bmax);
	/* maximum value for scalebar (default is the data maximum) */
    }

    if (!sf_getbool("allpos",&pos)) pos=true;
    /* contour positive values only */

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
		    nc = vp_optimal_scale(nc0, zmin-c0, zmax-c0, &zi, &dc);
		} else {
		    nc = vp_optimal_scale(nc0, zmin, zmax, &c0, &dc);
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
	    contour (pos, z,n1,n2,c[ic]);
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
	    
	    vp_barframe_init (bmin,bmax);
	    vp_barline(nc,c,bmin,bmax);
	}

    } /* i3 */

    exit(0);
}

/*
            north (0)
  (ix,iy+1) --------- (ix+1,iy+1)
            | cell  |
   west (3) | ix,iy | east (1)
            |       |
  (ix,iy)   --------- (ix+1,iy)
            south (2)
*/
static void contour (bool pos, float **z, int n1, int n2, float c)
{
    int ix, iy, non, jx, jy;
    float  zxy, ze, zn, x, y;
    bool **west, **south;

    west = sf_boolalloc2(n1,n2);
    south = sf_boolalloc2(n1,n2);

    /* find all the intersections */
    non = 0;      /* clear intersection counter */
    for (iy = 0; iy < n2 - 1; iy++) {
	for (ix = 0; ix < n1 - 1; ix++) {
	    zxy = z[iy][ix] - c;   /* z(x,y) - c */
	    ze = z[iy][ix+1] - c;  /* (z to the east) - c */
	    zn = z[iy+1][ix] - c;  /* (z to the north) - c */

	    if ((0. < zxy) != (0. < zn)) { /* if west edge intersected */
		west[iy][ix] = true;
		non++; 
	    } else { 
		west[iy][ix] = false;
	    }

	    if ((0. < zxy) != (0. < ze)) {  /* if south edge intersected */
		south[iy][ix] = true;
		non++; 
	    } else { 
		south[iy][ix] = false;
	    }
	}
    }

    for (ix = 0, iy = n2 - 1; ix < n1 - 1; ix++) { /* northern boundary */
	zxy = z[iy][ix]  - c;  /* z(x,y) - c */
	ze = z[iy][ix+1] - c;  /* (z to the east) - c */
	
	if ((0. < zxy) != (0. < ze)) { /* if south edge intersected */
	    south[iy][ix] = true;
	    non++; 
	} else { 
	    south[iy][ix] = false;
	}
	west[iy][ix] = false;
    }
    
    for (iy = 0, ix = n1 - 1; iy < n2 - 1; iy++) { /* eastern boundary */
	zxy = z[iy][ix]  - c;  /* z(x,y) - c */
	zn = z[iy+1][ix] - c;  /* (z to the north) - c */
	
	if ((0. < zxy) != (0. < zn)) { /* if west edge intersected */
	    west[iy][ix] = true;
	    non++; 
	} else { 
	    west[iy][ix] = false;
	}
	south[iy][ix] = false;
    }
    
    /* draw contours intersecting a boundary */
    for (ix = 0, iy = n2 - 1; ix < n1 - 1 && non > 0; ix++) { /* north */
	if (south[iy][ix]) {
	    south[iy][ix] = false;
	    non--;

	    x = ix + delta (c, z[iy][ix], z[iy][ix + 1]);
	    y = iy;

	    draw (false, x, y);

	    jx = ix;
	    jy = iy - 1;
	    while (cconnect (pos, south, west, z, n1, n2, c, &jx, &jy))
		non--;
	}
    }
    
    for (ix = n1 - 1, iy = 0; iy < n2 - 1 && non > 0; iy++) { /* east */
	if (west[iy][ix]) {
	    west[iy][ix] = false;
	    non--;

	    x = ix;
	    y = iy + delta (c, z[iy][ix], z[iy + 1][ix]);
	    draw (false, x, y);

	    jx = ix - 1;
	    jy = iy;
	    while (cconnect (pos, south, west, z, n1, n2, c, &jx, &jy))
		non--;
	}
    }

    for (ix = 0, iy = 0; ix < n1 - 1 && non > 0; ix++) {  /* south */
	if (south[iy][ix]) {
	    south[iy][ix] = false;
	    non--;

	    x = ix + delta (c, z[iy][ix], z[iy][ix + 1]);
	    y = iy;
	    draw (false, x, y);

	    jx = ix;
	    jy = iy;
	    while (cconnect (pos, south, west, z, n1, n2, c, &jx, &jy))
		non--;
	}
    }

    for (ix = 0, iy = 0; iy < n2 - 1 && non > 0; iy++) { /* west */
	if (west[iy][ix]) {
	    west[iy][ix] = false;
	    non--;

	    x = ix;
	    y = iy + delta (c, z[iy][ix], z[iy + 1][ix]);
	    draw (false, x, y);

	    jx = ix;
	    jy = iy;
	    while (cconnect (pos, south, west, z, n1, n2, c, &jx, &jy))
		non--;
	}
    }

    /* draw interior contours */ 
    for (iy = 0; iy < n2 - 1 && non > 0; iy++) {
	for (ix = 0; ix < n1 - 1 && non > 0; ix++) {
	    if (south[iy][ix]) {  /* check south edge of cell */
		south[iy][ix] = false;
		non--;

		x = ix + delta (c, z[iy][ix], z[iy][ix + 1]);
		y = iy;
		draw (false, x, y);

		jx = ix;
		jy = iy;
		south[iy][ix] = 
		    cconnect (pos, south, west, z, n1, n2, c, &jx, &jy);
		while (cconnect (pos, south, west, z, n1, n2, c, &jx, &jy))
		    non--;
	    }
	}
    }
}

/* 
   cconnect draws a line from one intersection of the cell (ix,iy)
   to another intersection of the cell, 
   provided the latter intersection exists,
   and then clears the latter intersection and updates ix and iy.
   cconnect returns false if the latter intersection does not exist or if the 
   latter intersection is a grid boundary; otherwise returns true.
*/
static bool cconnect (bool pos, bool **south, bool** west, float** z,
		      int n1, int n2, float  c, int *ix, int *iy)
{
    bool mask; 
    int jx, jy;
    float x, y;

    jx = (*ix);
    jy = (*iy);

    /**** mask ****/
    mask = !pos || (z[jy][jx] >= 0.   && z[jy][jx+1] >= 0. && 
		    z[jy+1][jx] >= 0. && z[jy+1][jx+1] >= 0.);

    if (south[jy + 1][jx]) {  /* if exiting north */
	jy++;
	south[jy][jx] = false;

	x = jx + delta (c, z[jy][jx], z[jy][jx + 1]);
	y = jy;
	draw (mask, x, y);

	if (++(*iy) >= n2 - 1) return false;
    } else if (west[jy][jx + 1]) { /* if exiting east */
	jx++;
	west[jy][jx] = false;

	x = jx;
	y = jy + delta (c, z[jy][jx], z[jy + 1][jx]);
	draw (mask, x, y);

	if (++(*ix) >= n1 - 1) return false;
    } else if (south[jy][jx]) {  /* if exiting south */
	south[jy][jx] = false;

	x = jx + delta (c, z[jy][jx], z[jy][jx + 1]);
	y = jy;
	draw (mask, x, y);

	if (--(*iy) < 0) return false;
    } else if (west[jy][jx]) { /* if exiting west */
	west[jy][jx] = false;
	
	x = jx;
	y = jy + delta (c, z[jy][jx], z[jy + 1][jx]);
	draw (mask, x, y);

	if (--(*ix) < 0) return false;
    } else {
	return false; /* no exit found */
    }
    
    return true;
}

/* compute (a-b)/(c-b) for use in linear interpolation */
static float delta (float a, float b, float  c)
{
    float t;

    t = c - b;      /* avoids division by zero */
    return (t != 0.0)? ((a - b) / t): 0.5;
}

static void draw (bool mask, float x, float y) {
    if (transp) {    
	if (mask) {
	    vp_udraw (o2+y*d2, o1+x*d1);
	} else {
	    vp_umove (o2+y*d2, o1+x*d1);
	}
    } else {
	if (mask) {
	    vp_udraw (o1+x*d1, o2+y*d2);
	} else {
	    vp_umove (o1+x*d1, o2+y*d2);
	}
    }
}

/* 	$Id$	 */
