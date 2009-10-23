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

typedef struct Velocity3 {
    float v0, gx, gy, gz, x0, y0, z0, vz, n;
} *velocity3;
/*^*/

#endif

static float fx, dx, fy, dy;
static velocity3 v;
static char type;
static float ***curve, ***dipx, ***dipy; 

void kirmod3_init(float fx1, float dx1,
		  float fy1, float dy1             /* reflector sampling */,
		  velocity3 v1                     /* velocity attributes */,
		  char type1                       /* velocity distribution */,
		  float ***curve1                  /* reflectors */,
		  float ***dipx1                   /* reflector inline dip */,
		  float ***dipy1                   /* reflector cross dip */)
/*< Initialize >*/ 
{
    dx = dx1; fx = fx1;
    dy = dy1; fy = fy1;

    v = v1;
    type = type1;
    curve = curve1;
    dipx = dipx1;
    dipy = dipy1;
}

void kirmod3_map(ktable ta          /* ray attribute object */,
		 const float *xy    /* surface position */,
		 int ix, int iy     /* point on the surface */,
		 int ic             /* surface number */) 
/*< Compute traveltimes and amplitudes >*/
{
    float x, y, z, x2, y2, zx, zy, x1, y1;
    float px, py, pz, v1, g, gy, gx, gz, dz;
    
    x1 = xy[0];
    y1 = xy[1];

    v1 = (v->v0)+
	(v->gy)*(y1-(v->y0))+
	(v->gx)*(x1-(v->x0))-
	(v->gz)*(v->z0);

    y2 = fy + iy*dy; /* y2 is on the reflector */
    y = y2 - y1;

    x2 = fx + ix*dx; /* x2 is on the reflector */
    x = x2 - x1;

    z = curve[ic][iy][ix];
    zx = dipx[ic][iy][ix];
    zy = dipy[ic][iy][ix];
    dz = sqrtf(1.0+zx*zx+zy*zy);
		    
    g = sqrtf((v->gz)*(v->gz)+
	      (v->gx)*(v->gx)+
	      (v->gy)*(v->gy));
    gx = v->gx+v->gz*zx;            /* dw/dx */
    gy = v->gy+v->gz*zy;
    gz = v->gz-v->gx*zx-v->gy*zy;
    px = x+z*zx;                    /* r*dr/dx */
    py = y+z*zy;
    pz = z-x*zx-y*zy;
    kirmod_table(type,false,z,x,y,g,gx,gy,gz,v1,v1,v->vz,v->n,px,py,pz,dz,ta);
}

