/* 3-D velocity grid for ray tracing in VTI. */
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

#include "grid3a.h"
#include "eno3.h"

#ifndef _grid3a_h

typedef struct Grid3a* grid3a;
/* abstract data type */
/*^*/

#endif

struct Grid3a {
    eno3 px, pz, pq;
    int n1, n2, n3;
    float o1, d1, o2, d2, o3, d3;
};
/* concrete data type */

grid3a grid3a_init (int n1, float o1, float d1 /* first axis */, 
		    int n2, float o2, float d2 /* second axis */,
		    int n3, float o3, float d3 /* second axis */,
		    float *vz2                 /* vertical velocity squared [n1*n2*n3], gets corrupted */, 
		    float *vx2                 /* horizontal velocity squared [n1*n2*n3], vz2 (isotropic) if NULL */, 
		    float *q                   /* nonellipticity [n1*n2*n3], 1 (elliptic) if NULL */,
		    int order                  /* interpolation order */)
/*< Initialize grid object >*/
{
    grid3a grd;
    int i;
    
    grd = (grid3a) sf_alloc(1,sizeof(*grd));

    grd->n1 = n1; grd->o1 = o1; grd->d1 = d1; 
    grd->n2 = n2; grd->o2 = o2; grd->d2 = d2;
    grd->n3 = n3; grd->o3 = o3; grd->d3 = d3;

    grd->pz = eno3_init (order, n1, n2, n3);
    grd->px = eno3_init (order, n1, n2, n3);
    grd->pq = eno3_init (order, n1, n2, n3);

    eno3_set1 (grd->pz, vz2);
    eno3_set1 (grd->px, (NULL==vx2)? vz2:vx2);

    if (NULL == q) {
	for (i=0; i < n1*n2; i++) {
	    vz2[i] = 1.;
	}
	eno3_set1 (grd->pq, vz2);
    } else {
	eno3_set1 (grd->pq, q);
    }

    return grd;
}

void grid3a_rhs(void* par /* grid */, 
		float* xy /* coordinate [6] */,
		float* g  /* right-hand side [6] */)
/*< right-hand side for the ray tracing system >*/
{
    grid3a grd;
    float x, y, z, vz, vx, q, q2, a, ax, az, e, d, r2, r, v, v2, num, den, one;
    float vz1[3], vx1[3], q1[3];
    int i, j, k;
    
    grd = (grid3a) par;
    
    x = (xy[0]-grd->o1)/grd->d1; i = x; x -= i;
    y = (xy[1]-grd->o2)/grd->d2; j = y; y -= j;
    z = (xy[2]-grd->o3)/grd->d3; k = z; z -= k;

    eno3_apply(grd->pz, i, j, k, x, y, z, &vz, vz1, BOTH);
    eno3_apply(grd->px, i, j, k, x, y, z, &vx, vx1, BOTH);
    eno3_apply(grd->pq, i, j, k, x, y, z, &q, q1, BOTH);

    one = sqrtf(xy[3]*xy[3]+xy[4]*xy[4]+xy[5]*xy[5]);
    a = xy[3]/one;
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
    
    g[3] = -(q1[0]*2.*vx*vz*ax*az +
	     vx1[0]*ax*(r-d+2*az*vz*q) +
	     vz1[0]*az*(r+d+2*ax*vx*q))/(4*r*v2*grd->d1);
    g[4] = -(q1[1]*2.*vx*vz*ax*az +
	     vx1[1]*ax*(r-d+2*az*vz*q) +
	     vz1[1]*az*(r+d+2*ax*vx*q))/(4*r*v2*grd->d2);
    g[5] = -(q1[2]*2.*vx*vz*ax*az +
	     vx1[2]*ax*(r-d+2*az*vz*q) +
	     vz1[2]*az*(r+d+2*ax*vx*q))/(4*r*v2*grd->d3);
    
    xy[3] /= v*one;
    xy[4] /= v*one;
    xy[5] /= v*one;

    num = d*d*r2+ax*az*q2*q2;
    den = d*d+2*(az-ax)*d*r+r2;

    /* group velocity */
    v *= 4*num/(r2*den);
    v = sqrtf(v);
    
    /* group angle */
    if (ax > az) {
	az = (d+r)*(d+r)*ax/(d*d+2*(ax-az)*d*r+r2); 
    } else {
	az = (d+r)*(d+r)*den/(16.*az*num);
    }    
    if (az > 1.) az=1.;
    else if (az < 0.) az=0.;

    g[0] = SF_SIG(xy[3])*v*sqrtf(az);
    g[1] = v*sqrtf(1.-az)*xy[4]/hypotf(xy[4],xy[5]);
    g[2] = v*sqrtf(1.-az)*xy[5]/hypotf(xy[4],xy[5]);
}

int grid3a_term (void* par /* grid */, 
		 float* xy /* location [3] */)
/*< Termination criterion. returns 0 if xy (data coordinates)
  are inside the grid >*/
{
    grid3a grd;
    
    grd = (grid3a) par;
    return (xy[0] < grd->o1 || xy[0] > grd->o1 + (grd->n1-1)*grd->d1 || 
	    xy[1] < grd->o2 || xy[1] > grd->o2 + (grd->n2-1)*grd->d2 ||
	    xy[2] < grd->o3 || xy[2] > grd->o3 + (grd->n3-1)*grd->d3);
}

void grid3a_close(grid3a grd)
/*< Free internal storage >*/
{
    eno3_close (grd->pz);
    eno3_close (grd->px);
    eno3_close (grd->pq);
    free (grd);
}

/* 	$Id$	 */
