/* Cell ray tracing. */
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

#include "celltrace.h"
#include "eno2.h"
#include "lsint2.h"
#include "alloc.h"
#include "cell.h"

#include "_bool.h"
/*^*/

#ifndef _sf_celltrace_h

typedef struct sf_CellTrace *sf_celltrace;
/* abstract data type */
/*^*/

#endif

struct sf_CellTrace {
    int nt, nx, nz;
    float dx, dz, x0, z0;
    sf_eno2 pnt;
    sf_lsint2 lnt;
};  


sf_celltrace sf_celltrace_init (bool lsint  /* use local ls interpolation */,
				int order   /* interpolation accuracy */, 
				int nt      /* maximum time steps */,
				int nz      /* depth samples */, 
				int nx      /* lateral samples */, 
				float dz    /* depth sampling */, 
				float dx    /* lateral sampling */, 
				float z0    /* depth origin */, 
				float x0    /* lateral origin */, 
				float* slow /* slowness [nz*nx] */)
/*< Initialize ray tracing object >*/
{
    sf_celltrace ct;

    ct = (sf_celltrace) sf_alloc (1,sizeof(*ct));

    ct->nt = nt;
    ct->dx = dx; ct->dz = dz;
    ct->nx = nx; ct->nz = nz;
    ct->x0 = x0; ct->z0 = z0;

    if (lsint) {
	ct->lnt = sf_lsint2_init (nz, nx); 
	sf_lsint2_set1 (ct->lnt, slow);
	ct->pnt = NULL;
    } else {
	ct->pnt = sf_eno2_init (order, nz, nx); 
	sf_eno2_set1 (ct->pnt, slow); 
	ct->lnt = NULL;
    }
    
    return ct;
} 

void sf_celltrace_close (sf_celltrace ct)
/*< Free allocated storage >*/
{
    if (NULL != ct->pnt) sf_eno2_close (ct->pnt);
    if (NULL != ct->lnt) sf_lsint2_close (ct->lnt);
    free (ct);
}

float sf_cell_trace (sf_celltrace ct, 
		     float* xp    /* position */, 
		     float* p     /* ray parameter */, 
		     int* it      /* steps till boundary */, 
		     float** traj /* trajectory */)
/*< ray trace >*/
{
    const float eps = 1.e-5;
    float t, v, x, z, s, sx, sz, g[2];
    int i, ix, iz, jx, jz;
    bool onx, onz;

    x = (xp[1]-ct->x0)/ct->dx; ix = floor(x); x -= ix;
    z = (xp[0]-ct->z0)/ct->dz; iz = floor(z); z -= iz;
    onx = sf_cell_snap (&x,&ix,eps);
    onz = sf_cell_snap (&z,&iz,eps);

    if (NULL != ct->lnt) {
	sf_lsint2_apply(ct->lnt,iz,ix,z,x,&v,g,BOTH);
    } else {
	sf_eno2_apply(ct->pnt,iz,ix,z,x,&v,g,BOTH);
    }
    g[1] /= ct->dx;
    g[0] /= ct->dz;

    t = sf_cell_update2 (2, 0., v, p, g);
    /* p is normal vector now ||p|| = 1 */
  
    if (traj != NULL) {
	traj[0][0] = xp[0];
	traj[0][1] = xp[1];
    }
  
    for (i=0; i < ct->nt; i++) {
	/* decide if we are out already */
	if (iz < 0 || iz > ct->nz-1 || 
	    (onz && iz == 0 && p[0] < 0.) ||
	    (iz == ct->nz-1 && (!onz || p[0] > 0.)) ||
	    ix < 0 || ix > ct->nx-1 || 
	    (onx && ix == 0 && p[1] < 0.) ||
	    (ix == ct->nx-1 && (!onx || p[1] > 0.))) break;

	sf_cell_intersect (g[1],x,ct->dx/v,p[1],&sx,&jx);
	sf_cell_intersect (g[0],z,ct->dz/v,p[0],&sz,&jz);

	s = SF_MIN(sx,sz);

	t += sf_cell_update1 (2, s, v, p, g);
	/* p is slowness vector now ||p||=v */

	if (s == sz) {
	    z = 0.; onz = true; iz += jz;
	    x += p[1]*s/ct->dx;
	    onx = sf_cell_snap (&x,&ix,eps);
	} else {
	    x = 0.; onx = true; ix += jx;
	    z += p[0]*s/ct->dz;
	    onz = sf_cell_snap (&z,&iz,eps);
	}

	if (traj != NULL) {
	    traj[i+1][0] = ct->z0 + (z+iz)*ct->dz;
	    traj[i+1][1] = ct->x0 + (x+ix)*ct->dx;
	}    
	
	if (NULL != ct->lnt) {
	    sf_lsint2_apply(ct->lnt,iz,ix,z,x,&v,g,BOTH);
	} else {
	    sf_eno2_apply(ct->pnt,iz,ix,z,x,&v,g,BOTH);
	}
	g[1] /= ct->dx;
	g[0] /= ct->dz;

	t += sf_cell_update2 (2, s, v, p, g);
	/* p is normal vector now ||p||=1 */
    }

    xp[1] = ct->x0 + (x+ix)*ct->dx;
    xp[0] = ct->z0 + (z+iz)*ct->dz;

    *it = (i > ct->nt)? -i:i;
    
    return t;
}

/* 	$Id: celltrace.c 9901 2013-02-14 22:49:26Z sfomel $	 */
