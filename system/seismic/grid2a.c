/* 2-D velocity grid for ray tracing in VTI. */
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

#include "grid2a.h"

#ifndef _grid2a_h

typedef struct Grid2a* grid2a;
/* abstract data type */
/*^*/

#endif

struct Grid2a {
    sf_eno2 px, pz, pq;
    int n1, n2;
    float o1, d1, o2, d2;
};
/* concrete data type */

grid2a grid2a_init (int n1, float o1, float d1 /* first axis */, 
		    int n2, float o2, float d2 /* second axis */,
		    float *vz2                 /* vertical velocity squared [n1*n2], gets corrupted */, 
		    float *vx2                 /* horizontal velocity squared [n1*n2], vz2 (isotropic) if NULL */, 
		    float *q                   /* nonellipticity [n1*n2], 1 (elliptic) if NULL */,
		    int order                  /* interpolation order */)
/*< Initialize grid object >*/
{
    grid2a grd;
    int i;
    
    grd = (grid2a) sf_alloc(1,sizeof(*grd));

    grd->n1 = n1; grd->o1 = o1; grd->d1 = d1; 
    grd->n2 = n2; grd->o2 = o2; grd->d2 = d2;
    
    grd->pz = sf_eno2_init (order, n1, n2);
    grd->px = sf_eno2_init (order, n1, n2);
    grd->pq = sf_eno2_init (order, n1, n2);

    sf_eno2_set1 (grd->pz, vz2);
    sf_eno2_set1 (grd->px, (NULL==vx2)? vz2:vx2);

    if (NULL == q) {
	for (i=0; i < n1*n2; i++) {
	    vz2[i] = 1.;
	}
	sf_eno2_set1 (grd->pq, vz2);
    } else {
	sf_eno2_set1 (grd->pq, q);
    }

    return grd;
}

void grid2a_rhs(void* par /* grid */, 
		float* xy /* coordinate [4] */,
		float* g)
/*< right-hand side for the ray tracing system >*/
{
    grid2a grd;
    float x, y, vz, vx, q, q2, a, ax, az, e, d, r2, r, v, v2, num, den, one;
    float vz1[2], vx1[2], q1[2];
    int i, j;
    
    grd = (grid2a) par;
    
    x = (xy[0]-grd->o1)/grd->d1; i = x; x -= i;
    y = (xy[1]-grd->o2)/grd->d2; j = y; y -= j;
    
    sf_eno2_apply(grd->pz, i, j, x, y, &vz, vz1, BOTH);
    sf_eno2_apply(grd->px, i, j, x, y, &vx, vx1, BOTH);
    sf_eno2_apply(grd->pq, i, j, x, y, &q, q1, BOTH);

    one = hypotf(xy[2],xy[3]);
    a = xy[2]/one;
    /* cosine of phase angle from the symmetry axis */
    az = a*a;  /* cosine squared */
    ax = 1-az; /* sine squared */

    e = vz*az+vx*ax; /* elliptic part */
    d = vz*az-vx*ax;
    q2 = q*vx*vz;

    r2 = d*d + 4.*ax*az*q2;
    r = sqrtf(r2);

    /* phase velocity */
    v2 = 0.5*(e+r);
    v = sqrtf(v2);

    g[2] = -(q1[0]*2.*vx*vz*ax*az +
	     vx1[0]*ax*(r-d+2*az*vz*q) +
	     vz1[0]*az*(r+d+2*ax*vx*q))/(4*r*v2*grd->d1);
    g[3] = -(q1[1]*2.*vx*vz*ax*az +
	     vx1[1]*ax*(r-d+2*az*vz*q) +
	     vz1[1]*az*(r+d+2*ax*vx*q))/(4*r*v2*grd->d2);
    
    xy[2] /= v*one;
    xy[3] /= v*one;

    num = d*d*r2+ax*az*q2*q2;
    den = d*d+2*(az-ax)*d*r+r2;

    /* group velocity */
    v2 *= 4*num/(r2*den);
    v = sqrtf(v2);
    
    /* group angle */
    if (ax > az) {
	az = (d+r)*(d+r)*ax/(d*d+2*(ax-az)*d*r+r2+SF_EPS); 
    } else {
	az = (d+r)*(d+r)*den/(16.*az*num);
    }    
    if (az > 1.) az=1.;
    else if (az < 0.) az=0.;

    g[0] = SF_SIG(xy[2])*v*sqrtf(az);
    g[1] = SF_SIG(xy[3])*v*sqrtf(1.-az);
}

int grid2a_term (void* par /* grid */, 
		 float* xy /* location [2] */)
/*< Termination criterion. returns 0 if xy (data coordinates)
  are inside the grid >*/
{
    grid2a grd;
    
    grd = (grid2a) par;
    return (xy[0] < grd->o1 || xy[0] > grd->o1 + (grd->n1-1)*grd->d1 || 
	    xy[1] < grd->o2 || xy[1] > grd->o2 + (grd->n2-1)*grd->d2);
}

void grid2a_close(grid2a grd)
/*< Free internal storage >*/
{
    sf_eno2_close (grd->pz);
    sf_eno2_close (grd->px);
    sf_eno2_close (grd->pq);
    free (grd);
}

/* 	$Id$	 */
