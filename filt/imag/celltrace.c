#include <math.h>

#include <rsf.h>

#include "celltrace.h"
#include "eno2.h"
#include "cell.h"

#ifndef MIN
#define MIN(a,b) ((a)<(b))?(a):(b)
#endif

struct CellTrace {
    int nt, nx, nz;
    float dx, dz, x0, z0;
    eno2 pnt;
};  

/* slow - slowness */
celltrace celltrace_init (int order, int nt,
			  int nz, int nx, 
			  float dz, float dx, 
			  float z0, float x0, float* slow)
{
    celltrace ct;

    ct = (celltrace) sf_alloc (1,sizeof(*ct));

    ct->nt = nt;
    ct->dx = dx; ct->dz = dz;
    ct->nx = nx; ct->nz = nz;
    ct->x0 = x0; ct->z0 = z0;
    ct->pnt = eno2_init (order, nz, nx); 
    eno2_set1 (ct->pnt, slow); 
    
    return ct;
} 

void celltrace_close (celltrace ct)
{
    eno2_close (ct->pnt);
    free (ct);
}

float cell_trace (celltrace ct, float* xp, float* p, int* it, float** traj)
{
    const float eps = 1.e-5;
    float t, v, x, z, s, sx, sz, g[2];
    int i, ix, iz, jx, jz;
    bool onx, onz;

    x = (xp[1]-ct->x0)/ct->dx; ix = floor(x); x -= ix;
    z = (xp[0]-ct->z0)/ct->dz; iz = floor(z); z -= iz;
    onx = cell_snap (&x,&ix,eps);
    onz = cell_snap (&z,&iz,eps);

    eno2_apply(ct->pnt,iz,ix,z,x,&v,g,BOTH);
    g[1] /= ct->dx;
    g[0] /= ct->dz;

    t = cell_update2 (2, 0., v, p, g);
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

	cell_intersect (g[1],x,ct->dx/v,p[1],&sx,&jx);
	cell_intersect (g[0],z,ct->dz/v,p[0],&sz,&jz);

	s = MIN(sx,sz);

	t += cell_update1 (2, s, v, p, g);
	/* p is slowness vector now ||p||=v */

	if (s == sz) {
	    z = 0.; onz = true; iz += jz;
	    x += p[1]*s/ct->dx;
	    onx = cell_snap (&x,&ix,eps);
	} else {
	    x = 0.; onx = true; ix += jx;
	    z += p[0]*s/ct->dz;
	    onz = cell_snap (&z,&iz,eps);
	}

	if (traj != NULL) {
	    traj[i+1][0] = ct->z0 + (z+iz)*ct->dz;
	    traj[i+1][1] = ct->x0 + (x+ix)*ct->dx;
	}    
	
	eno2_apply(ct->pnt,iz,ix,z,x,&v,g,BOTH);
	g[1] /= ct->dx;
	g[0] /= ct->dz;

	t += cell_update2 (2, s, v, p, g);
	/* p is normal vector now ||p||=1 */
    }

    xp[1] = ct->x0 + (x+ix)*ct->dx;
    xp[0] = ct->z0 + (z+iz)*ct->dz;

    *it = (i > ct->nt)? -i:i;
    
    return t;
}

/* 	$Id: celltrace.c,v 1.2 2003/09/30 14:30:52 fomels Exp $	 */
