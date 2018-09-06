/* 2-D Kirchhoff integral modeling */
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
#include <float.h>

#include <rsf.h>

#include "mpikirmod.h"
/*^*/

#ifndef _mpikirmod2_h

typedef struct Surface *surface;
/* abstract data type */
/*^*/

typedef struct Velocity {
    float v0, gx, gz, x0, z0, vz, n;
} *velocity;
/*^*/

#endif

struct Surface {
    int is, ih;
    float x;
    ktable **ta;
};

static int ny, nc, nx, **map;
static float x0, dx;

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

surface kirmod2_init(int ns,  float s0,  float ds  /* source/midpoint axis */,
		     int nh,  float h0,  float dh  /* offset axis */,
		     int nx1, float x01, float dx1 /* reflector axis */,
		     int nc1                       /* number of reflectors */,
                     bool cmp                      /* if CMP instead of shot gather */,
                     bool absoff                   /* use absolute offset */)
/*< Initialize surface locations >*/ 
{
    int is, ih, iy;
    float s, h;
    surface yi, y;

    nx = nx1;
    dx = dx1;
    x0 = x01;
    nc = nc1;

    if (cmp) {
	ny = 2*ns*nh;
	map = sf_intalloc2(2*nh,ns);
    } else {
	ny = ns*(nh+1);
	map = sf_intalloc2(nh+1,ns);
    }	

    y = (surface) sf_alloc(ny,sizeof(*y));

    yi = y;

    for (is=0; is < ns; is++) {
	s = s0 + is*ds;
	if (cmp) {
	    for (ih=0; ih < nh; ih++, yi++) {
		h = 0.5*(h0 + ih*dh);
		yi->x = s - h;
		yi->is = is;
		yi->ih = 2*ih;
		yi++;
		yi->x = s + h;
		yi->is = is;
		yi->ih = 2*ih+1;
	    }
	} else {
	    for (ih=0; ih < nh; ih++, yi++) {
		yi->x = absoff ? h0 + ih*dh : s + h0 + ih*dh;
		yi->is = is;
		yi->ih = ih;
	    }
	    yi->x = s;
	    yi->is = is;
	    yi->ih = nh;
	    yi++;
	}
    }

    qsort(y,ny,sizeof(*y),surface_comp);

    for (iy=0; iy < ny; iy++) {
	yi = y+iy;
	map[yi->is][yi->ih] = iy;
    }

    return y;
}

void kirmod2_close(surface y) 
/*< Free allocated storage >*/
{
    int iy, ix, ic;
    ktable **ta, **ta2;

    ta2 = y[0].ta;;
    for (iy=0; iy < ny; iy++) {
	ta = y[iy].ta;
	if (ta != ta2) {
	    for (ix=0; ix < nx; ix++) {
		for (ic=0; ic < nc; ic++) { 
		    free(ta2[ix][ic]);
		}
		free(ta2[ix]);
	    }
	    free(ta2);
	    ta2 = NULL;
	}
	ta2 = ta;
    }
    
    if (ta2 != NULL) {
	for (ix=0; ix < nx; ix++) {
	    for (ic=0; ic < nc; ic++) { 
		free(ta2[ix][ic]);
	    }
	    free(ta2[ix]);
	}
	free(ta2);
	ta2 = NULL;
    }

    free (y);
    free (map[0]);
    free (map);
}

void kirmod2_table (surface y                  /* surface structure */,
		    velocity v                 /* velocity attributes */, 
		    char type                  /* velocity distribution */,
		    bool twod                  /* 2-D or 2.5-D */, 
		    float **curve              /* reflectors */,
		    float **dip                /* reflector dip */)
/*< Compute traveltime/amplitude map >*/
{
    int ix, iy, ic;
    float x, z, x2, zx, x1, xp=0., px, pz, v1, v2, g, gx, gz, dz;
    ktable **ta=NULL;

    for (iy=0; iy < ny; iy++) {	
	x1 = y[iy].x; /* x1 is on the surface */
	if (0==iy || x1 != xp) { /* new point */
	    v1 = (v->v0)+(v->gx)*(x1-(v->x0))-(v->gz)*(v->z0);
	    ta = (ktable**) sf_alloc(nx,sizeof(ktable*));
	    
	    for (ix=0; ix < nx; ix++) {
		ta[ix] = (ktable*) sf_alloc(nc,sizeof(ktable));

		x2 = x0 + ix*dx; /* x2 is on the reflector */
		x = x2 - x1;

		for (ic=0; ic < nc; ic++) { 
		    ta[ix][ic] = (ktable) sf_alloc(1,sizeof(ta[ix][ic][0]));

		    z = curve[ic][ix];
		    zx = dip[ic][ix];
		    dz = hypotf(1.0,zx);
		    v2 = (v->v0)+(v->gx)*(x2-(v->x0))+(v->gz)*(z-(v->z0));
		    
		    g = hypotf(v->gx,v->gz);        /* gradient */
		    gx = (v->gx)+(v->gz)*zx;        /* dw/dx */
		    gz = (v->gz)-(v->gx)*zx;
		    px = x+z*zx;                    /* r*dr/dx */
		    pz = z-x*zx;
		    kirmod_table(type,twod,
				 z,x,0.,g,gx,0.,gz,v1,v2,v->vz,v->n,px,0.,pz,dz,ta[ix][ic]);
		} 
	    } 
	} 
	y[iy].ta = ta;
	xp = x1;
    }
}

ktable kirmod2_map(surface y, int is, int ih, int ix, int ic) 
/*< Extract from traveltime/amplitude map >*/
{
    return y[map[is][ih]].ta[ix][ic];
}
