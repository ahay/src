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
#include <float.h>

#include <rsf.h>

#include "kirmod.h"

#ifndef _kirmod_h

typedef struct Surface *surface;
/* abstract data type */
/*^*/

#endif

struct Surface {
    int is, ih;
    float x, ***ta;
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
    float ***ta, ***ta2;

    ta2 = y[0].ta;;
    for (iy=0; iy < ny; iy++) {
	ta = y[iy].ta;
	if (ta != ta2) {
	    free(**ta2);
	    free(*ta2);
	    free(ta2);
	    ta2 = NULL;
	}
	ta2 = ta;
    }
    
    if (ta2 != NULL) {
	free(**ta2);
	free(*ta2);
	free(ta2);
	ta2 = NULL;
    }

    free (y);
    free (map[0]);
    free (map);
}

void kirmod_table (surface y                  /* surface structure */,
		   char type                  /* velocity distribution */, 
		   int nx, float x0, float dx /* reflector axis */,
		   int nc                     /* number of reflectors */,
		   float **curve              /* reflectors */,
		   float **dip                /* reflector dip */,
		   float *veloc               /* velocity attributes */) 
/*< Compute traveltime/amplitude map >*/
{
    int ix, iy, ic;
    float x, z, x2, zx, x1, xp=0., ***ta=NULL, t, a, p, q, r, v0, v1, v2;
    float sigma, rad, g, gx, gy, dz;

    for (iy=0; iy < ny; iy++) {	
	x1 = y[iy].x; /* x1 is on the surface */
	if (0==iy || x1 != xp) { /* new point */
	    ta = sf_floatalloc3(4,nc,nx);

	    for (ix=0; ix < nx; ix++) {
		x2 = x0 + ix*dx; /* x2 is on the reflector */
		x = x2 - x1;

		for (ic=0; ic < nc; ic++) { 
		    z = curve[ic][ix];
		    zx = dip[ic][ix];
		    dz = hypotf(1.0,zx);

		    r = hypotf(x,z)+FLT_EPSILON*dx; /* distance */
		    g = hypotf(veloc[1],veloc[2]);  /* gradient */
		    gx = veloc[2]+veloc[1]*zx;      /* dw/dx */
		    gy = veloc[1]-veloc[2]*zx;
		    p = x+z*zx;                     /* r*dr/dx */
		    q = z-x*zx;
		    
		    switch (type) {
			case 'c': /* constant velocity */
			    v0 = veloc[0];
			    a = r;
			    t = r/v0;
			    p /= (r*v0);
			    q = fabsf(q)/(r*dz);
			    break;		    
			case 's': /* linear sloth */
			    v1 = veloc[0]+2*(veloc[2]*(x1-veloc[4])-veloc[1]*veloc[3]);
			    v2 = veloc[0]+2*(veloc[2]*(x2-veloc[4])+veloc[1]*(z-veloc[3]));
			    v0 = 0.5*(v1+v2);
			    rad = v0*v0-r*r*g*g;           /* => v0^2=1/v^4 */
			    if (rad < 0.) 
				sf_error("%s: shadow zone, fix later",__FILE__);
			    rad = sqrtf(rad);              /* => v0=1/v^2 */
			    sigma = r*sqrtf(2./(v0+rad));  /* => r*v */
			    a = sigma*sqrtf(rad);          /* => r */
			    t = sigma*(2*v0+rad)/3.;       /* => r/v */
			    q = fabsf(q+0.5*sigma*sigma*gy)/(sqrtf(v2)*sigma*dz);
			    p = (p/sigma+0.5*gx*sigma);    /* => p/(rv) */
			    break;
			case 'v': /* linear velocity */
			    v1 = veloc[0]+veloc[2]*(x1-veloc[4])-veloc[1]*veloc[3];
			    v2 = veloc[0]+veloc[2]*(x2-veloc[4])+veloc[1]*(z-veloc[3]);
			    v0 = sqrtf(v1*v2);                   /* => v */
			    a = r*hypotf(1.,0.5*r*g/v0);         /* => r */
			    t = acoshf(1.+0.5*r*r*g*g/(v0*v0))/g;    /* => r/v */
			    q = fabsf(q*v2-0.5*r*r*gy)/(v0*a*dz);
			    p = (p-0.5*r*r*gx/v2)/(a*v0);       /* => p/(r*v) */
			    break;
			default:
			    a = 0.;
			    t = 0.;
			    p = 0.;
			    q = 0.;
			    sf_error("%s: case %c is not implemented",__FILE__,type);
			    break;
		    } /* type */
		    ta[ix][ic][0] = t;                                   /* traveltime */
		    ta[ix][ic][1] = a;                                   /* geometrical spreading */
		    ta[ix][ic][2] = p;                                   /* traveltime slope */
		    ta[ix][ic][3] = (q >= 1.0)? 0.: -acosf(q)*SF_SIG(p); /* angle from the normal */
		} /* ic */
	    } /* ix */
	} /* if new */
	y[iy].ta = ta;
	xp = x1;
    }
}

float **kirmod_map(surface y, int is, int ih, int ix) 
/*< Extract from traveltime/amplitude map >*/
{
    return y[map[is][ih]].ta[ix];
}
