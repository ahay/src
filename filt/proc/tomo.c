#include <math.h>

#include <rsf.h>

#include "tomo.h"

/* velocity grid dimensions */
static int nz, nx, np;
static float dz, dx, dp, p0;

/* initializer (constructor) 
np1 - number of slopes
n1, n2 - grid dimensions
d1, d2 - cell size
*/
void tomo_init(int np1, int n1, int n2, float d1, float d2)
{
  nz = n1;
  nx = n2;
  dz = d1;
  dx = d2;

  np = np1;
  dp = 2./(np-1);
  p0 = -1.;
}

/* linear operator: forward and ajoint
--------------------------------------
adj: adjoint flag (either t = L s (adj=false) or s = L' t (adj=true) 
add: addition flag (if true, t = t + L v or s = s + L' t)
ns: total size of s (ns = nz*nx) nz is the fast axis
nt: size of t (nt = np*nx) np is the fast axis
s: slowness
t: traveltime
*/
void tomo_lop(bool adj, bool add, int ns, int nt, float* s, float* t)
{
  int is, ix, iz, ip, iy, it;
  float p, x, deltax, distance;

  if (ns != nx*nz || nt != nx*np) sf_error("%s: wrong size",__FILE__);

  /* initialize for lop
     if !add && adj, zero s
     if !add && !adj, zero t
     if add, do nothing */
  sf_adjnull(adj,add,ns,nt,s,t);

  for (iy=0; iy < nx; iy++) { /* loop over sources */
    for (ip=0; ip < np; ip++) { /* loop over initial directions */
      p = p0 + ip*dp; /* initial slope */
      x = iy*dx; /* initial position (origin iz zero) */

      is = iy*nz + nz-1; /* where in the grid is [nz-1,is] */
      it = iy*np + ip;

      deltax = dz*p; /* shift in x */
      distance = hypotf(dz,deltax); /* Pythagor */

      for (iz=nz-1; iz > 0; iz--) { /* loop up in depth */
	if (adj) {
	  s[is] += t[it]*distance;
	} else {
	  t[it] += s[is]*distance; /* time is slowness * distance */
	}

	x += deltax;
	ix = 0.5 + x/dx; /* nearest point on the grid */

	if (ix < 0 || ix >= nx) break; /* outside the grid */

	is = ix*nz + iz;       
      } /* iz */
    } /* ip */
  } /* is */
}
      
