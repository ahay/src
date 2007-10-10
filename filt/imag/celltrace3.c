/* Cell ray tracing in 3-D. */
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

#include "celltrace3.h"

#include "eno3.h"

#ifndef _celltrace3_h

typedef struct Celltrace3 *celltrace3;
/* abstract data type */
/*^*/

#endif

struct Celltrace3 {
    int nt, nx, ny, nz;
    float dx, dy, dz, x0, y0, z0;
    eno3 pnt;
};  


celltrace3 celltrace3_init (int order   /* interpolation accuracy */, 
			    int nt      /* maximum time steps */,
			    int nz      /* depth samples */, 
			    int ny      /* inline samples */,
			    int nx      /* crossline samples */, 
			    float dz    /* depth sampling */,
			    float dy    /* inline sampling */,
			    float dx    /* crossline sampling */, 
			    float z0    /* depth origin */,
			    float y0    /* inline origin */,
			    float x0    /* crossline origin */, 
			    float* slow /* slowness [nz*nx] */)
/*< Initialize ray tracing object >*/
{
    celltrace3 ct;

    ct = (celltrace3) sf_alloc (1,sizeof(*ct));

    ct->nt = nt;
    ct->dx = dx; ct->dy = dy; ct->dz = dz; 
    ct->nx = nx; ct->ny = ny; ct->nz = nz;
    ct->x0 = x0; ct->y0 = y0; ct->z0 = z0;
    ct->pnt = eno3_init (order, nz, ny, nx); 
    eno3_set1 (ct->pnt, slow); 
    
    return ct;
} 

void celltrace3_close (celltrace3 ct)
/*< Free allocated storage >*/
{
    eno3_close (ct->pnt);
    free (ct);
}

float cell_trace3 (celltrace3 ct, 
		  float* xp    /* position */, 
		  float* p     /* ray parameter */, 
		  int* it      /* steps till boundary */, 
		  float** traj /* trajectory */)
/*< ray trace >*/
{
    const float eps = 1.e-5;
    float t, v, x, y, z, s, sx, sy, sz, g[3];
    int i, ix, iy, iz, jx, jy, jz;
    bool onx, ony, onz;

    x = (xp[2]-ct->x0)/ct->dx; ix = floorf(x); x -= ix;
    y = (xp[1]-ct->y0)/ct->dy; iy = floorf(y); y -= iy;
    z = (xp[0]-ct->z0)/ct->dz; iz = floorf(z); z -= iz;
    onx = sf_cell_snap (&x,&ix,eps);
    ony = sf_cell_snap (&y,&iy,eps);
    onz = sf_cell_snap (&z,&iz,eps);

    eno3_apply(ct->pnt,iz,iy,ix,z,y,x,&v,g,BOTH);
    g[2] /= ct->dx;
    g[1] /= ct->dy;
    g[0] /= ct->dz;

    t = sf_cell_update2 (3, 0., v, p, g);
    /* p is normal vector now ||p|| = 1 */
  
    if (traj != NULL) {
	traj[0][0] = xp[0];
	traj[0][1] = xp[1];
	traj[0][2] = xp[2];
    }
  
    for (i=0; i < ct->nt; i++) {
	/* decide if we are out already */
	if (iz < 0 || iz > ct->nz-1 || 
	    (onz && iz == 0 && p[0] < 0.) ||
	    (iz == ct->nz-1 && (!onz || p[0] > 0.)) ||
	    iy < 0 || iy > ct->ny-1 || 
	    (ony && iy == 0 && p[1] < 0.) ||
	    (iy == ct->ny-1 && (!ony || p[1] > 0.)) ||
	    ix < 0 || ix > ct->nx-1 || 
	    (onx && ix == 0 && p[1] < 0.) ||
	    (ix == ct->nx-1 && (!onx || p[2] > 0.))) break;

	/* interesection with cell walls */
	sf_cell_intersect (g[2],x,ct->dx/v,p[2],&sx,&jx);
	sf_cell_intersect (g[1],y,ct->dy/v,p[1],&sy,&jy);
	sf_cell_intersect (g[0],z,ct->dz/v,p[0],&sz,&jz);

	/* choose first intersection */
	s = SF_MIN(SF_MIN(sx,sy),sz);

	t += sf_cell_update1 (3, s, v, p, g);
	/* p is slowness vector now ||p||=v */

	if (s == sz) {
	    z = 0.; onz = true; iz += jz;
	    x += p[2]*s/ct->dx;
	    onx = sf_cell_snap (&x,&ix,eps);
	} else if (s == sy) {
	    y = 0.; ony = true; iy += jy;
	    y += p[1]*s/ct->dy;
	    ony = sf_cell_snap (&y,&iy,eps);  
	} else {
	    x = 0.; onx = true; ix += jx;
	    z += p[0]*s/ct->dz;
	    onz = sf_cell_snap (&z,&iz,eps);
	}

	if (traj != NULL) { /* save trajectory */
	    traj[i+1][0] = ct->z0 + (z+iz)*ct->dz;
	    traj[i+1][1] = ct->y0 + (y+iy)*ct->dy;
	    traj[i+1][2] = ct->x0 + (x+ix)*ct->dx;
	}    
	
	/* interpolate velocity and gradient from the grid */
	eno3_apply(ct->pnt,iz,iy,ix,z,y,x,&v,g,BOTH);
	g[1] /= ct->dx;
	g[2] /= ct->dy;
	g[0] /= ct->dz;

	t += sf_cell_update2 (3, s, v, p, g);
	/* p is normal vector now ||p||=1 */
    }

    xp[2] = ct->x0 + (x+ix)*ct->dx;
    xp[1] = ct->y0 + (y+iy)*ct->dy;
    xp[0] = ct->z0 + (z+iz)*ct->dz;

    *it = (i > ct->nt)? -i:i;
    
    return t;
}

/* 	$Id$	 */
