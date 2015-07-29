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

#ifndef _vp_contour_h

typedef struct sf_Contour *vp_contour;
/* abstract data type */
/*^*/

#endif

struct sf_Contour {
    int n1, n2;
    float o1, o2, d1, d2, s1, s2;
    bool transp, **west, **south;
};
/* concrete data type */

static float delta (float a, float b, float  c);
static bool cconnect (vp_contour cnt, 
		      bool pos, float** z, float  c, int *ix, int *iy);
static void draw (vp_contour cnt, bool mask, float x, float y);

vp_contour vp_contour_init(bool transp,                          /* transpose flag */
			   int n1, float o1, float d1, float s1, /* first axis */
			   int n2, float o2, float d2, float s2) /* second axis */
/*< initialize >*/
{
    vp_contour cnt;

    cnt = (vp_contour) sf_alloc(1,sizeof(*cnt));
    
    cnt->transp = transp;
    cnt->n1 = n1; cnt->o1 = o1; cnt->d1 = d1; cnt->s1 = s1;
    cnt->n2 = n2; cnt->o2 = o2; cnt->d2 = d2; cnt->s2 = s2;

    cnt->west = sf_boolalloc2(n1,n2);
    cnt->south = sf_boolalloc2(n1,n2);

    return cnt;
}

void vp_contour_close(vp_contour cnt)
/*< free allocated storage >*/
{
    free(*(cnt->west)); free(cnt->west);
    free(*(cnt->south)); free(cnt->south);
    free(cnt);
}

void vp_contour_draw (vp_contour cnt,
		      bool pos  /* positive only flag */, 
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
    for (iy = 0; iy < cnt->n2 - 1; iy++) {
	for (ix = 0; ix < cnt->n1 - 1; ix++) {
	    zxy = z[iy][ix] - c;   /* z(x,y) - c */
	    ze = z[iy][ix+1] - c;  /* (z to the east) - c */
	    zn = z[iy+1][ix] - c;  /* (z to the north) - c */

	    if ((0. < zxy) != (0. < zn)) { /* if west edge intersected */
		cnt->west[iy][ix] = true;
		non++; 
	    } else { 
		cnt->west[iy][ix] = false;
	    }

	    if ((0. < zxy) != (0. < ze)) {  /* if south edge intersected */
		cnt->south[iy][ix] = true;
		non++; 
	    } else { 
		cnt->south[iy][ix] = false;
	    }
	}
    }

    for (ix = 0, iy = cnt->n2 - 1; ix < cnt->n1 - 1; ix++) { /* northern boundary */
	zxy = z[iy][ix]  - c;  /* z(x,y) - c */
	ze = z[iy][ix+1] - c;  /* (z to the east) - c */
	
	if ((0. < zxy) != (0. < ze)) { /* if south edge intersected */
	    cnt->south[iy][ix] = true;
	    non++; 
	} else { 
	    cnt->south[iy][ix] = false;
	}
	cnt->west[iy][ix] = false;
    }
    
    for (iy = 0, ix = cnt->n1 - 1; iy < cnt->n2 - 1; iy++) { /* eastern boundary */
	zxy = z[iy][ix]  - c;  /* z(x,y) - c */
	zn = z[iy+1][ix] - c;  /* (z to the north) - c */
	
	if ((0. < zxy) != (0. < zn)) { /* if west edge intersected */
	    cnt->west[iy][ix] = true;
	    non++; 
	} else { 
	    cnt->west[iy][ix] = false;
	}
	cnt->south[iy][ix] = false;
    }
    
    /* draw contours intersecting a boundary */
    for (ix = 0, iy = cnt->n2 - 1; ix < cnt->n1 - 1 && non > 0; ix++) { /* north */
	if (cnt->south[iy][ix]) {
	    cnt->south[iy][ix] = false;
	    non--;

	    x = ix + delta (c, z[iy][ix], z[iy][ix + 1]);
	    y = iy;

	    draw (cnt, false, x, y);

	    jx = ix;
	    jy = iy - 1;
	    while (cconnect (cnt, pos, z, c, &jx, &jy)) non--;
	}
    }
    
    for (ix = cnt->n1 - 1, iy = 0; iy < cnt->n2 - 1 && non > 0; iy++) { /* east */
	if (cnt->west[iy][ix]) {
	    cnt->west[iy][ix] = false;
	    non--;

	    x = ix;
	    y = iy + delta (c, z[iy][ix], z[iy + 1][ix]);
	    draw (cnt, false, x, y);

	    jx = ix - 1;
	    jy = iy;
	    while (cconnect (cnt, pos, z, c, &jx, &jy)) non--;
	}
    }

    for (ix = 0, iy = 0; ix < cnt->n1 - 1 && non > 0; ix++) {  /* south */
	if (cnt->south[iy][ix]) {
	    cnt->south[iy][ix] = false;
	    non--;

	    x = ix + delta (c, z[iy][ix], z[iy][ix + 1]);
	    y = iy;
	    draw (cnt, false, x, y);

	    jx = ix;
	    jy = iy;
	    while (cconnect (cnt, pos, z, c, &jx, &jy)) non--;
	}
    }

    for (ix = 0, iy = 0; iy < cnt->n2 - 1 && non > 0; iy++) { /* west */
	if (cnt->west[iy][ix]) {
	    cnt->west[iy][ix] = false;
	    non--;

	    x = ix;
	    y = iy + delta (c, z[iy][ix], z[iy + 1][ix]);
	    draw (cnt, false, x, y);

	    jx = ix;
	    jy = iy;
	    while (cconnect (cnt, pos, z, c, &jx, &jy)) non--;
	}
    }

    /* draw interior contours */ 
    for (iy = 0; iy < cnt->n2 - 1 && non > 0; iy++) {
	for (ix = 0; ix < cnt->n1 - 1 && non > 0; ix++) {
	    if (cnt->south[iy][ix]) {  /* check south edge of cell */
		cnt->south[iy][ix] = false;
		non--;

		x = ix + delta (c, z[iy][ix], z[iy][ix + 1]);
		y = iy;
		draw (cnt, false, x, y);

		jx = ix;
		jy = iy;
		cnt->south[iy][ix] = cconnect (cnt, pos, z, c, &jx, &jy);
		while (cconnect (cnt, pos, z, c, &jx, &jy)) non--;
	    }
	}
    }
}


static bool cconnect (vp_contour cnt, 
		      bool pos, float** z, float  c, int *ix, int *iy)
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
    mask = (bool) (!pos || (z[jy][jx] >= 0.   && z[jy][jx+1] >= 0. && 
			    z[jy+1][jx] >= 0. && z[jy+1][jx+1] >= 0.));

    if (cnt->south[jy + 1][jx]) {  /* if exiting north */
	jy++;
	cnt->south[jy][jx] = false;

	x = jx + delta (c, z[jy][jx], z[jy][jx + 1]);
	y = jy;
	draw (cnt, mask, x, y);

	if (++(*iy) >= cnt->n2 - 1) return false;
    } else if (cnt->west[jy][jx + 1]) { /* if exiting east */
	jx++;
	cnt->west[jy][jx] = false;

	x = jx;
	y = jy + delta (c, z[jy][jx], z[jy + 1][jx]);
	draw (cnt, mask, x, y);

	if (++(*ix) >= cnt->n1 - 1) return false;
    } else if (cnt->south[jy][jx]) {  /* if exiting south */
	cnt->south[jy][jx] = false;

	x = jx + delta (c, z[jy][jx], z[jy][jx + 1]);
	y = jy;
	draw (cnt, mask, x, y);

	if (--(*iy) < 0) return false;
    } else if (cnt->west[jy][jx]) { /* if exiting west */
	cnt->west[jy][jx] = false;
	
	x = jx;
	y = jy + delta (c, z[jy][jx], z[jy + 1][jx]);
	draw (cnt, mask, x, y);

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

static void draw (vp_contour cnt, bool mask, float x, float y) 
{
    float x2, y2;

    x2 = cnt->o1+(x+y*cnt->s1)*cnt->d1;
    y2 = cnt->o2+(y+x*cnt->s2)*cnt->d2;

    if (cnt->transp) {    
	if (mask) {
	    vp_udraw (y2,x2);
	} else {
	    vp_umove (y2,x2);
	}
    } else {
	if (mask) {
	    vp_udraw (x2,y2);
	} else {
	    vp_umove (x2,y2);
	}
    }
}

/* 	$Id: contour.c 1778 2006-03-23 13:52:59Z fomels $	 */
