#include <math.h>

#include <rsf.h>

#include "ztrace.h"
#include "cell.h"
#include "eno2.h"
#include "grad2fill.h"

#ifndef MIN
#define MIN(a,b) ((a)<(b))?(a):(b)
#endif

#ifndef MAX
#define MAX(a,b) ((a)>(b))?(a):(b)
#endif

static int nx, nz, na, nax, nt;
static float dx,dz,da, x0,z0,a0;
static eno2 cvel, fslc[NS];
static float **slice;
static bool *known;
static const float eps = 1.e-5;

static void psnap (float* p, float* q, int* iq);

void ztrace_init (int order, int iorder,
		  int nx1, int nz1, int na1, int nt1,
		  float dx1, float dz1, float da1,
		  float x01, float z01, float a01,
		  float** vel, float** slice_in)
{
    int is;

    slice = slice_in;
    nx = nx1; nz = nz1; na = na1;
    dx = dx1; dz = dz1; da = da1;
    x0 = x01; z0 = z01; a0 = a01;
    nax = na*nx; nt = nt1;

    cvel = eno2_init (order, nz, nx);
    eno2_set (cvel, vel); /* vel is slowness */

    known = sf_boolalloc (nax);

    grad2fill_init (na, nx);

    for (is=0; is < NS; is++) {
	fslc[is] = eno2_init (iorder, na, nx);
    }
}

void ztrace_close (void)
{
    int is;

    free (known);

    eno2_close (cvel);
    
    for (is = 0; is < NS; is++) {
	eno2_close (fslc[is]);
    }

    grad2fill_close ();
}

void ztrace_step (int kz) {
    float v, v0, g[2], g0[2], t, z, x, a, p[2], q, s, sx, sz;
    int is, it, nk, kx, kp, k, ix, iz, ip, jx, jz;
    bool onx, onz;
    
    for (is=0; is < NS; is++) {
	eno2_set1 (fslc[is], slice[is]);
    }

    nk = 0; /* number of known */
    
    /* First pass */
    for (kx=0; kx < nx; kx++) {
	eno2_apply(cvel,kz,kx,0.,0.,&v0,g0,BOTH);
	g0[1] /= dx;
	g0[0] /= dz;
	
	for (kp=0; kp < na; kp++) {
	    k = kp + kx*na;

	    a = a0+kp*da;
	    p[1] = sin(a);
	    p[0] = -cos(a);

	    v = v0;
	    g[0] = g0[0];
	    g[1] = g0[1];

	    ix = kx; x = 0.; onx = true;
	    iz = kz; z = 0.; onz = true;

	    t = cell_update2 (2, 0., v, p, g);
	    /* p is normal vector now ||p|| = 1 */

	    /* decide if we are out already */
	    if ((iz == 0 && p[0] < 0.) ||
		(iz == nz-1 && p[0] > 0.) ||
		(ix == 0 && p[1] < 0.) ||
		(ix == nx-1 && p[1] > 0.)) {
		slice[0][k] = 0.;
		slice[1][k] = x0+ix*dx;
		slice[2][k] = z0+iz*dz;
		slice[3][k] = cell_p2a(p)*180./SF_PI;
		known[k] = true;
		nk++;
		continue;
	    } else {
		known[k] = false;
	    }

	    for (it=0; it < nt; it++) {
		cell_intersect (g[1],x,dx/v,p[1],&sx,&jx);
		cell_intersect (g[0],z,dz/v,p[0],&sz,&jz);

		s = MIN(sx,sz);

		t += cell_update1 (2, s, v, p, g);
		/* p is slowness vector now ||p||=v */

		if (s == sz) {
		    z = 0.; onz = true; iz += jz;
		    x += p[1]*s/dx;
		    onx = cell_snap (&x,&ix,eps);
		} else {
		    x = 0.; onx = true; ix += jx;
		    z += p[0]*s/dz;
		    onz = cell_snap (&z,&iz,eps);
		}
		
		eno2_apply(cvel,iz,ix,z,x,&v,g,BOTH);
		g[1] /= dx;
		g[0] /= dz;

		t += cell_update2 (2, s, v, p, g);
		/* p is normal vector now ||p||=1 */

		/* decide if we are out */
		if (iz < 0 || iz > nz-1 || 
		    (onz && iz == 0 && p[0] < 0.) ||
		    (iz == nz-1 && (!onz || p[0] > 0.)) ||
		    ix < 0 || ix > nx-1 || 
		    (onx && ix == 0 && p[1] < 0.) ||
		    (ix == nx-1 && (!onx || p[1] > 0.))) {
		    slice[0][k] = t;
		    slice[1][k] = x0+(x+ix)*dx;
		    slice[2][k] = z0+(z+iz)*dz;
		    slice[3][k] = cell_p2a(p)*180./SF_PI;
		    known[k] = true;
		    nk++;
		    break;
		} else if (onz && iz == kz-1) { /* went up */
		    psnap (p,&q,&ip);
		    for (is=0; is < NS; is++) {
			eno2_apply(fslc[is],ip,ix,q,x,&slice[is][k],NULL,FUNC);
		    }
		    slice[0][k] += t;
		    known[k] = true;
		    nk++;
		    break;
		}
	    } /* it */
	} /* kp */
    } /* kx */

    if (nk < nax) {
	fprintf(stderr,"known=%d (%d)\n",nk,nax);
	nk = MIN(MIN(na,nx),nax-nk);
	for (is=0; is < NS; is++) {
	    grad2fill (nk, slice[is], known);
	}
    }
}

static void psnap (float* p, float* q, int* iq) {
    int ia;
    float a2, a;

    a = (cell_p2a(p)-a0)/da;
    ia = floor (a); a2 = a-ia;
    cell_snap (&a2, &ia, eps);

    if (ia < 0) {
	ia=0; a2=0.; a=a0; 
    } else if (ia >= na || (ia == na-1 && a2 > 0.)) {
	ia=na-1; a2=0.; a=a0+(na-1)*da; 
    } 

    p[1] = sin(a);
    p[0] = -cos(a);

    *q = a2;
    *iq = ia;
}
