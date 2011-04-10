#include <float.h>
#include <math.h>

#include <rsf.h>

#include "coh.h"
#include "apfilt.h"

static int window (int j, int w, int n);

struct Coh {
    int nx, ny, nz, nw, nj;
    float*** pp;
};

coh coh_init(int nw, int nj, int nx, int ny, int nz, float ***pp)
{
    coh ap;

    ap = (coh) sf_alloc(1,sizeof(*ap));

    ap->nw = nw;
    ap->nj = nj;
    ap->nx = nx;
    ap->ny = ny;
    ap->nz = nz;
    ap->pp = pp;

    return ap;
}

static int window (int j, int w, int n) {
    j -= w/2; 
    if      (j < 0)     j = 0;
    else if (j + w > n) j = n - w;
    return j;
}

float coh1 (const coh ap, int jx, int jy, int jz, int nx, int ny, int nz,
	    float*** xx)
{
    int ix, iy, iz, iw, is, xmin, xmax, np;
    float f[7], ab, a2, b2, a, b, p;

    jy = window(jy,ny,ap->ny-1);
    jz = window(jz,nz,ap->nz);
        
    xmin = window(jx,nx,ap->nx);
    xmax = xmin+nx;

    if (xmin < ap->nw*ap->nj) xmin = ap->nw*ap->nj;
    if (xmax > ap->nx-ap->nw*ap->nj) xmax=ap->nx-ap->nw*ap->nj;

    /* find the average dip in window */
    p = 0.;
    np = 0;
    for (iz=jz; iz < jz+nz; iz++) {
	for (iy=jy; iy < jy+ny; iy++) {	    
	    for (ix = xmin; ix < xmax; ix++) {
		p += ap->pp[iz][iy][ix];
		np++;
	    }
	}
    }
    if (np) p /= np;

    ab = a2 = b2 = 0.;
    passfilter(ap->nw, p, f);

    for (iz=jz; iz < jz+nz; iz++) {
	for (iy=jy; iy < jy+ny; iy++) {	    
	    for (ix = xmin; ix < xmax; ix++) {
		a = b = 0.;
		for (iw = 0; iw <= 2*ap->nw; iw++) {
		    is = (iw-ap->nw)*ap->nj;
		    a += f[iw]*xx[iz][iy+1][ix+is];
		    b += f[iw]*xx[iz][iy  ][ix-is];
		}
		ab += a*b;
		a2 += a*a;
		b2 += b*b;
	    }
	}
    }
    
    return ab/(sqrtf(a2*b2)+FLT_EPSILON);
}

float coh2 (const coh ap, int jx, int jy, int jz, int nx, int ny, int nz,
	    float*** xx)
{
    int ix, iy, iz, iw, is, xmin, xmax, np;
    float f[7], ab, a2, b2, a, b, p;

    jy = window(jy,ny,ap->ny);
    jz = window(jz,nz,ap->nz-1);
        
    xmin = window(jx,nx,ap->nx);
    xmax = xmin+nx;

    if (xmin < ap->nw*ap->nj) xmin = ap->nw*ap->nj;
    if (xmax > ap->nx-ap->nw*ap->nj) xmax=ap->nx-ap->nw*ap->nj;

    /* find the average dip in window */
    p = 0.;
    np = 0;
    for (iz=jz; iz < jz+nz; iz++) {
	for (iy=jy; iy < jy+ny; iy++) {	    
	    for (ix = xmin; ix < xmax; ix++) {
		p += ap->pp[iz][iy][ix];
		np++;
	    }
	}
    }
    if (np) p /= np;

    ab = a2 = b2 = 0.;
    passfilter(ap->nw, p, f);

    for (iz=jz; iz < jz+nz; iz++) {
	for (iy=jy; iy < jy+ny; iy++) {	    
	    for (ix = xmin; ix < xmax; ix++) {
		a = b = 0.;
		for (iw = 0; iw <= 2*ap->nw; iw++) {
		    is = (iw-ap->nw)*ap->nj;
		    a += f[iw]*xx[iz+1][iy][ix+is];
		    b += f[iw]*xx[iz  ][iy][ix-is];
		}
		ab += a*b;
		a2 += a*a;
		b2 += b*b;
	    }
	}
    }
    
    return ab/(sqrt(a2*b2)+FLT_EPSILON);
}

/* 	$Id$	 */
