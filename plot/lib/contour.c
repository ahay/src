/* Contour plot. */
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
/*^*/

#include "vplot.h"

static float delta (float a, float b, float  c);
static bool cconnect (bool pos, float** z, float  c, int *ix, int *iy);
static void draw (bool mask, float x, float y);

static bool transp, **west, **south;
static int n1, n2;
static float o1, o2, d1, d2;

void vp_contour_init(bool transp_in,                      /* transpose flag */
		     int n1_in, float o1_in, float d1_in, /* first axis */
		     int n2_in, float o2_in, float d2_in) /* second axis */
/*< initialize >*/
{
    transp = transp_in;
    n1 = n1_in; o1 = o1_in; d1 = d1_in;
    n2 = n2_in; o2 = o2_in; d2 = d2_in;

    west = sf_boolalloc2(n1,n2);
    south = sf_boolalloc2(n1,n2);
}

void vp_contour_close(void)
/*< free allocated storage >*/
{
    free(*west); free(west);
    free(*south); free(south);
}

void vp_contour (bool pos  /* positive only flag */, 
		 float **z /* data to contour */, 
		 float c   /* contour level */)
/*< draw contours >*/
{
    int ix, iy, non, jx, jy;
    float  zxy, ze, zn, x, y;

/*
            north (0)
  (ix,iy+1) --------- (ix+1,iy+1)
            | cell  |
   west (3) | ix,iy | east (1)
            |       |
  (ix,iy)   --------- (ix+1,iy)
            south (2)
*/

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
	    while (cconnect (pos, z, c, &jx, &jy)) non--;
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
	    while (cconnect (pos, z, c, &jx, &jy)) non--;
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
	    while (cconnect (pos, z, c, &jx, &jy)) non--;
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
	    while (cconnect (pos, z, c, &jx, &jy)) non--;
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
		south[iy][ix] = cconnect (pos, z, c, &jx, &jy);
		while (cconnect (pos, z, c, &jx, &jy)) non--;
	    }
	}
    }
}


static bool cconnect (bool pos, float** z, float  c, int *ix, int *iy)
/* 
   draws a line from one intersection of the cell (ix,iy)
   to another intersection of the cell, 
   provided the latter intersection exists,
   and then clears the latter intersection and updates ix and iy.
   cconnect returns false if the latter intersection does not exist or if the 
   latter intersection is a grid boundary; otherwise returns true.
*/
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

/* 	$Id: contour.c 1778 2006-03-23 13:52:59Z fomels $	 */
