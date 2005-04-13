/* Kirchhoff integral modeling */
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
#include <stdlib.h>

#include <rsf.h>

#include "kirmod.h"

#ifndef _kirmod_h

typedef enum {CONST, S2LIN, VLIN} maptype;
/* type of velocity distribution */
/*^*/

typedef struct Surface *surface;
/* abstract data type */
/*^*/

#endif

struct Surface {
    int is, ih;
    float x, **ta;
};


static int ny, **map;

static int surface_comp(const void *a, const void *b)
/* compare by the surface coordinate */
{
    float aa, bb;

    aa = ((surface) a)->x;
    bb = ((surface) b)->x;

    if (aa <  bb) return (-1);
    if (aa == bb) return 0;
    return 1;
}

surface kirmod_init(int ns, float s0, float ds /* source axis */,
		    int nh, float h0, float dh /* offset axis */)
/*< Initialize surface locations >*/ 
{
    int is, ih, iy;
    float s;
    surface yi, y;

    ny = ns*(nh+1);
    y = (surface) sf_alloc(ny,sizeof(*y));
    map = sf_intalloc2(nh+1,ns);

    yi = y;
    for (is=0; is < ns; is++, yi++) {
	s = s0 + is*ds;
	for (ih=0; ih < nh; ih++, yi++) {
	    yi->x = s + h0 + ih*dh;
	    yi->is = is;
	    yi->ih = ih;
	}
	yi->x = s;
	yi->is = is;
	yi->ih = nh;	
    }

    qsort(y,ny,sizeof(*y),surface_comp);

    for (iy=0; iy < ny; iy++) {
	yi = y+iy;
	map[yi->is][yi->ih] = iy;
    }

    return y;
}

void kirmod_close(surface y) 
/*< Free allocated storage >*/
{
    int iy;
    float **ta, **ta2;

    ta2 = y[0].ta;;
    for (iy=0; iy < ny; iy++) {
	ta = y[iy].ta;
	if (ta != ta2) {
	    free(*ta2);
	    free(ta2);
	    ta2 = NULL;
	}
	ta2 = ta;
    }
    
    if (ta2 != NULL) {
	free(*ta2);
	free(ta2);
	ta2 = NULL;
    }

    free (y);
    free (map[0]);
    free (map);
}

void kirmod_table (surface y                  /* surface structure */,
		   maptype type               /* velocity distribution */, 
		   int nx, float x0, float dx /* reflector axis */,
		   float *curve               /* reflector */,
		   float *dip                 /* reflector dip */,
		   float *veloc               /* velocity attributes */) 
/*< Compute traveltime/amplitude map >*/
{
    int ix, iy;
    float x, z, x1, x2=0., **ta=NULL, t, a, p;

    for (iy=0; iy < ny; iy++) {	
	x1 = y[iy].x;
	if (0==iy || x1 != x2) {
	    ta = sf_floatalloc2(3,nx);

	    for (ix=0; ix < nx; ix++) {
		x = x0-x1 + ix*dx;
		z = curve[ix];
		switch (type) {
		    case CONST:	
			a = hypotf(x,z);
			t = a/veloc[0];
			p = (x+z*dip[ix])/(a*veloc[0]);
			break;		    
		    default:
			a = 0.;
			t = 0.;
			p = 0.;
			sf_error("__FILE__: case %d is not implemented",type);
			break;
		}
		ta[ix][0] = t; /* traveltime */
		ta[ix][1] = a; /* amplitude */
		ta[ix][2] = p; /* slope */
	    }
	}
	y[iy].ta = ta;
	x2 = x1;
    }
}

float *kirmod_map(surface y, int is, int ih, int ix) 
/*< Extract from traveltime/amplitude map >*/
{
    return y[map[is][ih]].ta[ix];
}
