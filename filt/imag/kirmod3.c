/* 3-D Kirchhoff integral modeling */
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

#include "kirmod3.h"

#include "kirmod.h"
/*^*/

#ifndef _kirmod3_h

typedef struct Surface3 *surface3;
/* abstract data type */
/*^*/

typedef struct Velocity3 {
    float v0, gx, gy, gz, x0, y0, z0;
} *velocity3;
/*^*/

#endif

struct Surface3 {
    int isx, isy, ihx, ihy;
    float x, y;
    ktable ***ta;
};

static int nm, nc, nx, ny, nh, ***map;
static float fx, dx, fy, dy;

static int surface3_comp(const void *a, const void *b)
/* compare by the surface coordinate */
{
    float ax, bx, ay, by;

    ax = ((surface3) a)->x;
    bx = ((surface3) b)->x;
    ay = ((surface3) a)->y;
    by = ((surface3) b)->y;

    if (ax <  bx) return (-1);
    if (ax >  bx) return 1;
    if (ay <  by) return (-1);
    if (ay == by) return 0;
    return 1;
}

surface3 kirmod3_init(int nsx, float s0x, float dsx  /* source inline */,
		      int nsy, float s0y, float dsy  /* source crossline */,
		      int nhx, float h0x, float dhx  /* offset inline */,
		      int nhy, float h0y, float dhy  /* offset crossline */,
		      int nx1, float fx1, float dx1  /* reflector inline */,
		      int ny1, float fy1, float dy1  /* reflector crossline */,
		      int nc1                        /* number of reflectors */)
/*< Initialize surface locations >*/ 
{
    int isx, isy, ihx, ihy, im;
    float sx, sy;
    surface3 yi, y;

    nx = nx1; dx = dx1; fx = fx1;
    ny = ny1; dy = dy1; fy = fy1;
    nc = nc1;
    nh = nhx;

    nm = nsx*nsy*(nhx*nhy+1);
    y = (surface3) sf_alloc(nm,sizeof(*y));
    map = sf_intalloc3(nhx*nhy+1,nsx,nsy);

    yi = y;
    for (isy=0; isy < nsy; isy++) {
	sy = s0y + isy*dsy;
	for (isx=0; isx < nsx; isx++, yi++) {
	    sx = s0x + isx*dsx;
	    ihx=0;
	    for (ihy=0; ihy < nhy; ihy++) {
		for (ihx=0; ihx < nhx; ihx++, yi++) {
		    yi->x = sx + h0x + ihx*dhx;
		    yi->y = sy + h0y + ihy*dhy;
		    yi->isx = isx;
		    yi->isy = isy;
		    yi->ihx = ihx;
		    yi->ihy = ihy;
		}
	    }
	    yi->x = sx;
	    yi->y = sy;
	    yi->isx = isx;
	    yi->isy = isy;
	    yi->ihx = ihx;
	    yi->ihy = ihy;
	}
    }

    qsort(y,nm,sizeof(*y),surface3_comp);

    for (im=0; im < nm; im++) {
	yi = y+im;
	map[yi->isy][yi->isx][yi->ihy*nhx+yi->ihx] = im;
    }

    return y;
}

void kirmod3_close(surface3 y) 
/*< Free allocated storage >*/
{
    int im, ix, iy, ic;
    ktable ***ta, ***ta2;

    ta2 = y[0].ta;;
    for (im=0; im < nm; im++) {
	ta = y[im].ta;
	if (ta != ta2) {
	    for (iy=0; iy < ny; iy++) {
		for (ix=0; ix < nx; ix++) {
		    for (ic=0; ic < nc; ic++) { 
			free(ta2[iy][ix][ic]);
		    }
		    free(ta2[iy][ix]);
		}
		free(ta2[iy]);
	    }
	    free(ta2);
	    ta2 = NULL;
	}
	ta2 = ta;
    }
    
    if (ta2 != NULL) {
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		for (ic=0; ic < nc; ic++) { 
		    free(ta2[iy][ix][ic]);
		}
		free(ta2[iy][ix]);
	    }
	    free(ta2[iy]);
	}
	free(ta2);
	ta2 = NULL;
    }

    free (y);
    free (map[0]);
    free (map);
}

void kirmod3_table (surface3 s                 /* surface structure */,
		    velocity3 v                /* velocity attributes */,
		    char type                  /* velocity distribution */,
		    float ***curve             /* reflectors */,
		    float ***dipx              /* reflector inline dip */,
		    float ***dipy              /* reflector crossline dip */)
/*< Compute traveltime/amplitude map >*/
{
    int ix, iy, im, ic;
    float x, y, z, x2, y2, zx, zy, x1, y1, px, py, pz, r, v1, g, gy, gx, gz, dz;
    float xp=0., yp=0.;
    ktable ***ta=NULL;

    for (im=0; im < nm; im++) {	
	x1 = s[im].x; /* x1 is on the surface */
	y1 = s[im].y; /* x1 is on the surface */
	if (0==im || x1 != xp || y1 != yp) { /* new point */
	    v1 = 
		(v->v0)+
		(v->gy)*(y1-(v->y0))+
		(v->gx)*(x1-(v->x0))-
		(v->gz)*(v->z0);
	    ta = (ktable***) sf_alloc(ny,sizeof(ktable**));
	    
	    for (iy=0; iy < ny; iy++) {
		ta[iy] = (ktable**) sf_alloc(nx,sizeof(ktable*));
		
		y2 = fy + iy*dy; /* y2 is on the reflector */
		y = y2 - y1;

		for (ix=0; ix < nx; ix++) {
		    ta[iy][ix] = (ktable*) sf_alloc(nc,sizeof(ktable));

		    x2 = fx + ix*dx; /* x2 is on the reflector */
		    x = x2 - x1;

		    for (ic=0; ic < nc; ic++) { 
			ta[iy][ix][ic] = (ktable) 
			    sf_alloc(1,sizeof(ta[iy][ix][0]));
			
			z = curve[ic][iy][ix];
			zx = dipx[ic][iy][ix];
			zy = dipy[ic][iy][ix];
			dz = sqrtf(1.0+zx*zx+zy*zy);
		    
			r = sqrtf(x*x+y*y+z*z)+
			    FLT_EPSILON*hypotf(dx,dy); /* distance */
			g = sqrtf((v->gz)*(v->gz)+
				  (v->gx)*(v->gx)+
				  (v->gy)*(v->gy));
			gx = v->gx+v->gz*zx;            /* dw/dx */
			gy = v->gy+v->gz*zy;
			gz = v->gz-v->gx*zx-v->gy*zy;
			px = x+z*zx;                    /* r*dr/dx */
			py = y+z*zy;
			pz = z-x*zx-y*zy;
			kirmod_table(type,false,r,g,gx,gy,gz,v1,v1,px,py,pz,dz,
				     ta[iy][ix][ic]);
		    } 
		}
	    } 
	} 
	s[im].ta = ta;
	xp = x1;
	yp = y1;
    }
}

ktable kirmod3_map(surface3 s, int isx, int isy, int ihx, int ihy, 
		   int ix, int iy, int ic) 
/*< Extract from traveltime/amplitude map >*/
{
    return s[map[isy][isx][ihy*nh+ihx]].ta[iy][ix][ic];
}
