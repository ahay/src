/* 3-D dip estimation (including the non-stationary version) */
/*
  Copyright (C) 2020 University of Texas at Austin
  
  On Mar, 4, 2020, Austin
  
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
#include <rsf.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "divnn.h"
#include "diputil.h"


static int Nf; /*size of allpass filter*/
static double *b;

void apfilt_init(int nw /* filter order */)
/*< initialize >*/
{
    int j, k;
    double bk;

    Nf = nw*2;
    b = (double*) sf_alloc(Nf+1,sizeof(double));

    for (k=0; k <= Nf; k++) {
	bk = 1.0;
	for (j=0; j < Nf; j++) {
	    if (j < Nf-k) {
		bk *= (k+j+1.0)/(2*(2*j+1)*(j+1));
	    } else {
		bk *= 1.0/(2*(2*j+1));
	    }
	}
	b[k] = bk;
    }
}

void apfilt_close(void)
/*< free allocated storage >*/
{
    free(b);
}

void passfilter (float p  /* slope */, 
		 float* a /* output filter [n+1] */)
/*< find filter coefficients >*/
{
    int j, k;
    double ak;
    
    for (k=0; k <= Nf; k++) {
	ak = b[k];
	for (j=0; j < Nf; j++) {
	    if (j < Nf-k) {
		ak *= (Nf-j-p);
	    } else {
		ak *= (p+j+1);
	    }
	}
	a[k] = ak;
    }
}

void aderfilter (float p  /* slope */, 
		 float* a /* output filter [n+1] */)
/*< find coefficients for filter derivative >*/
{

    int i, j, k;
    double ak, ai;
    
    for (k=0; k <= Nf; k++) {
	ak = 0.;
	for (i=0; i < Nf; i++) {
	    ai = -1.0;
	    for (j=0; j < Nf; j++) {
		if (j != i) {			
		    if (j < Nf-k) {
			ai *= (Nf-j-p);
		    } else {
			ai *= (p+j+1);
		    }
		} else if (j < Nf-k) {
		    ai *= (-1);
		}
	    }
	    ak += ai;
	}
	a[k] = ak*b[k];
    }
}

#ifndef _diputil_h

typedef struct Allpass *allpass;
/* abstract data type */
/*^*/

#endif

struct Allpass {
    int nx, ny, nz, nw, nj;
    bool drift;
    float *flt, *pp;
};

static allpass ap1, ap2;

allpass allpass_init(int nw                 /* filter size */, 
		     int nj                 /* filter step */, 
		     int nx, int ny, int nz /* data size */, 
		     bool drift             /* if shift filter */,
		     float *pp              /* dip [nz*ny*nx] */)
/*< Initialize >*/
{
    allpass ap;

    ap = (allpass) sf_alloc(1,sizeof(*ap));

    ap->nw = nw;
    ap->nj = nj;
    ap->nx = nx;
    ap->ny = ny;
    ap->nz = nz;
    ap->drift = drift;
    ap->pp = pp;

    ap->flt = sf_floatalloc(2*nw+1);
    apfilt_init(nw);

    return ap;
}

void allpass_close(allpass ap)
/*< free allocated storage >*/
{
    apfilt_close();
    free(ap->flt);
    free(ap);
}

void allpass1 (bool left        /* left or right prediction */,
	       bool der         /* derivative flag */, 
	       const allpass ap /* PWD object */, 
	       float* xx        /* input */, 
	       float* yy        /* output */)
/*< in-line plane-wave destruction >*/
{
    int ix, iy, iz, iw, is, i, nx, ny, nz, i1, i2, ip, id;

    nx = ap->nx;
    ny = ap->ny;
    nz = ap->nz;

    if (left) {
	i1=1; i2=ny;   ip=-nx;
    } else {
	i1=0; i2=ny-1; ip=nx;
    }

    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		i = ix + nx * (iy + ny * iz);
		yy[i] = 0.;
	    }
	}
    }
  
    for (iz=0; iz < nz; iz++) {
	for (iy=i1; iy < i2; iy++) {
	    for (ix = ap->nw*ap->nj; ix < nx-ap->nw*ap->nj; ix++) {
		i = ix + nx * (iy + ny * iz);

		if (ap->drift) {
		    id = SF_NINT(ap->pp[i]);
		    if (ix-ap->nw*ap->nj-id < 0 || 
			ix+ap->nw*ap->nj-id >= nx) continue;

		    if (der) {
			aderfilter(ap->pp[i]-id, ap->flt);
		    } else {
			passfilter(ap->pp[i]-id, ap->flt);
		    }
		    
		    for (iw = 0; iw <= 2*ap->nw; iw++) {
			is = (iw-ap->nw)*ap->nj;
			
			yy[i] += (xx[i+is+ip] - xx[i-is-id]) * ap->flt[iw];
		    }		    
		} else {
		    if (der) {
			aderfilter(ap->pp[i], ap->flt);
		    } else {
			passfilter(ap->pp[i], ap->flt);
		    }
		    
		    for (iw = 0; iw <= 2*ap->nw; iw++) {
			is = (iw-ap->nw)*ap->nj;
			
			yy[i] += (xx[i+is+ip] - xx[i-is]) * ap->flt[iw];
		    }
		}
	    }
	}
    }
}

void allpass1_fb(bool left        /* left or right prediction */,
	       bool der         /* derivative flag */, 
	       const allpass ap /* PWD object */, 
	       float* xx        /* input */, 
	       float* yy        /* output */)
/*< in-line plane-wave destruction with forward and backward space derivative calculation >*/
{
    int ix, iy, iz, iw, is, i, nx, ny, nz, i1, i2, id;

    nx = ap->nx;
    ny = ap->ny;
    nz = ap->nz;

/*modified by YC on May, 12, 2020, for forward and backward space derivative calculation*/
/*    if (left) {
	i1=1; i2=ny;   ip=-nx;
    } else {
	i1=0; i2=ny-1; ip=nx;
    }*/
    i1=1;i2=ny-1;

    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		i = ix + nx * (iy + ny * iz);
		yy[i] = 0.;
	    }
	}
    }
  
    for (iz=0; iz < nz; iz++) {
	for (iy=i1; iy < i2; iy++) {
	    for (ix = ap->nw*ap->nj; ix < nx-ap->nw*ap->nj; ix++) {
		i = ix + nx * (iy + ny * iz);

		if (ap->drift) {
		    id = SF_NINT(ap->pp[i]);
		    if (ix-ap->nw*ap->nj-id < 0 || 
			ix+ap->nw*ap->nj-id >= nx) continue;

		    if (der) {
			aderfilter(ap->pp[i]-id, ap->flt);
		    } else {
			passfilter(ap->pp[i]-id, ap->flt);
		    }
		    
		    for (iw = 0; iw <= 2*ap->nw; iw++) {
			is = (iw-ap->nw)*ap->nj;
			yy[i] += 0.5*(xx[i+is+nx] - xx[i-is-id]) * ap->flt[iw]+0.5*(xx[i+is] - xx[i-is-id-nx]) * ap->flt[iw];
		    }		    
		} else {
		    if (der) {
			aderfilter(ap->pp[i], ap->flt);
		    } else {
			passfilter(ap->pp[i], ap->flt);
		    }
		    
		    for (iw = 0; iw <= 2*ap->nw; iw++) {
			is = (iw-ap->nw)*ap->nj;
/*modified by YC on May, 12, 2020, for forward and backward space derivative calculation*/
			yy[i] += 0.5*(xx[i+is+nx] - xx[i-is]) * ap->flt[iw]+0.5*(xx[i+is] - xx[i-is-nx]) * ap->flt[iw];
/* 			yy[i] += (xx[i+is+ip] - xx[i-is]) * ap->flt[iw]; */
		    }
		}
	    }
	}
    }
}



void allpass1t (bool left        /* left or right prediction */,
	       bool der         /* derivative flag */, 
	       const allpass ap /* PWD object */, 
	       float* xx        /* input */, 
	       float* yy        /* output */)
/*< adjoint of in-line plane-wave destruction >*/
{
    int ix, iy, iz, iw, is, i, nx, ny, nz, i1, i2, ip, id;

    nx = ap->nx;
    ny = ap->ny;
    nz = ap->nz;

    if (left) {
	i1=1; i2=ny;   ip=-nx;
    } else {
	i1=0; i2=ny-1; ip=nx;
    }

    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		i = ix + nx * (iy + ny * iz);
		xx[i] = 0.;
	    }
	}
    }
  
    for (iz=0; iz < nz; iz++) {
	for (iy=i1; iy < i2; iy++) {
	    for (ix = ap->nw*ap->nj; ix < nx-ap->nw*ap->nj; ix++) {
		i = ix + nx * (iy + ny * iz);

		if (ap->drift) {
		    id = SF_NINT(ap->pp[i]);
		    if (ix-ap->nw*ap->nj-id < 0 || 
			ix+ap->nw*ap->nj-id >= nx) continue;

		    if (der) {
			aderfilter(ap->pp[i]-id, ap->flt);
		    } else {
			passfilter(ap->pp[i]-id, ap->flt);
		    }
		    
		    for (iw = 0; iw <= 2*ap->nw; iw++) {
			is = (iw-ap->nw)*ap->nj;
			
			xx[i+is+ip] += yy[i] * ap->flt[iw];
			xx[i-is-id] -= yy[i] * ap->flt[iw];
		    }
		} else {
		    if (der) {
			aderfilter(ap->pp[i], ap->flt);
		    } else {
			passfilter(ap->pp[i], ap->flt);
		    }
		    
		    for (iw = 0; iw <= 2*ap->nw; iw++) {
			is = (iw-ap->nw)*ap->nj;
			
			xx[i+is+ip] += yy[i] * ap->flt[iw];
			xx[i-is]    -= yy[i] * ap->flt[iw];
		    }
		}
	    }
	}
    }
}

void left1 (bool left        /* left or right prediction */,
	       bool der         /* derivative flag */, 
	       const allpass ap /* PWD object */, 
	       float* xx        /* input */, 
	       float* yy        /* output */)
/*< left part of in-line plane-wave destruction >*/
{
    int ix, iy, iz, iw, is, i, nx, ny, nz, i1, i2, ip, id;

    nx = ap->nx;
    ny = ap->ny;
    nz = ap->nz;

    if (left) {
	i1=1; i2=ny;   ip=-nx;
    } else {
	i1=0; i2=ny-1; ip=nx;
    }

    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		i = ix + nx * (iy + ny * iz);
		yy[i] = 0.;
	    }
	}
    }
  
    for (iz=0; iz < nz; iz++) {
	for (iy=i1; iy < i2; iy++) {
	    for (ix = ap->nw*ap->nj; ix < nx-ap->nw*ap->nj; ix++) {
		i = ix + nx * (iy + ny * iz);

		if (ap->drift) {
		    id = SF_NINT(ap->pp[i]);
		    if (ix-ap->nw*ap->nj-id < 0 || 
			ix+ap->nw*ap->nj-id >= nx) continue;

		    if (der) {
			aderfilter(ap->pp[i]-id, ap->flt);
		    } else {
			passfilter(ap->pp[i]-id, ap->flt);
		    }
		} else {
		    if (der) {
			aderfilter(ap->pp[i], ap->flt);
		    } else {
			passfilter(ap->pp[i], ap->flt);
		    }
		}

		for (iw = 0; iw <= 2*ap->nw; iw++) {
		    is = (iw-ap->nw)*ap->nj;
		    
		    yy[i] += xx[i+is+ip] * ap->flt[iw];
		}
	    }
	}
    }
}

void right1 (bool left        /* left or right prediction */,
	       bool der         /* derivative flag */, 
	       const allpass ap /* PWD object */, 
	       float* xx        /* input */, 
	       float* yy        /* output */)
/*< right part of in-line plane-wave destruction >*/
{
    int ix, iy, iz, iw, is, i, nx, ny, nz, i1, i2, id;

    nx = ap->nx;
    ny = ap->ny;
    nz = ap->nz;

    if (left) {
	i1=1; i2=ny;   
    } else {
	i1=0; i2=ny-1;
    }

    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		i = ix + nx * (iy + ny * iz);
		yy[i] = 0.;
	    }
	}
    }
  
    for (iz=0; iz < nz; iz++) {
	for (iy=i1; iy < i2; iy++) {
	    for (ix = ap->nw*ap->nj; ix < nx-ap->nw*ap->nj; ix++) {
		i = ix + nx * (iy + ny * iz);

		if (ap->drift) {
		    id = SF_NINT(ap->pp[i]);
		    if (ix-ap->nw*ap->nj-id < 0 || 
			ix+ap->nw*ap->nj-id >= nx) continue;
		    
		    if (der) {
			aderfilter(ap->pp[i]-id, ap->flt);
		    } else {
			passfilter(ap->pp[i]-id, ap->flt);
		    }
		    
		    for (iw = 0; iw <= 2*ap->nw; iw++) {
			is = (iw-ap->nw)*ap->nj;
			
			yy[i] += xx[i-is-id] * ap->flt[iw];
		    }
		} else {
		    if (der) {
			aderfilter(ap->pp[i], ap->flt);
		    } else {
			passfilter(ap->pp[i], ap->flt);
		    }
		    
		    for (iw = 0; iw <= 2*ap->nw; iw++) {
			is = (iw-ap->nw)*ap->nj;
			
			yy[i] += xx[i-is] * ap->flt[iw];
		    }
		}
	    }
	}
    }
}

void allpass2 (bool left        /* left or right prediction */,
	       bool der         /* derivative flag */, 
	       const allpass ap /* PWD object */, 
	       float* xx        /* input */, 
	       float* yy        /* output */)
/*< cross-line plane-wave destruction >*/
{
    int ix, iy, iz, iw, is, i, nx, ny, nz, i1, i2, ip, id;

    nx = ap->nx;
    ny = ap->ny;
    nz = ap->nz;

    if (left) {
	i1=1; i2=nz;   ip=-nx*ny;
    } else {
	i1=0; i2=nz-1; ip=nx*ny;
    }
    
    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		i = ix + nx * (iy + ny * iz);
		yy[i] = 0.;
	    }
	}
    }
    
    for (iz=i1; iz < i2; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix = ap->nw*ap->nj; ix < nx-ap->nw*ap->nj; ix++) {
		i = ix + nx * (iy + ny * iz);
		
		if (ap->drift) {
		    id = SF_NINT(ap->pp[i]);
		    if (ix-ap->nw*ap->nj-id < 0 || 
			ix+ap->nw*ap->nj-id >= nx) continue;

		    if (der) {
			aderfilter(ap->pp[i]-id, ap->flt);
		    } else {
			passfilter(ap->pp[i]-id, ap->flt);
		    }
		    
		    for (iw = 0; iw <= 2*ap->nw; iw++) {
			is = (iw-ap->nw)*ap->nj;
			
			yy[i] += (xx[i+is+ip] - xx[i-is-id]) * ap->flt[iw];
		    }

		} else {
		    if (der) {
			aderfilter(ap->pp[i], ap->flt);
		    } else {
			passfilter(ap->pp[i], ap->flt);
		    }
		    
		    for (iw = 0; iw <= 2*ap->nw; iw++) {
			is = (iw-ap->nw)*ap->nj;
			
			yy[i] += (xx[i+is+ip] - xx[i-is]) * ap->flt[iw];
		    }
		}
	    }
	}
    }
}

void allpass2_fb (bool left        /* left or right prediction */,
	       bool der         /* derivative flag */, 
	       const allpass ap /* PWD object */, 
	       float* xx        /* input */, 
	       float* yy        /* output */)
/*< cross-line plane-wave destruction with forward and backward space derivative calculation>*/
{
    int ix, iy, iz, iw, is, i, nx, ny, nz, i1, i2, ip, id;

    nx = ap->nx;
    ny = ap->ny;
    nz = ap->nz;

    if (left) {
	i1=1; i2=nz;   ip=-nx*ny;
    } else {
	i1=0; i2=nz-1; ip=nx*ny;
    }
    
    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix=0; ix < nx; ix++) {
		i = ix + nx * (iy + ny * iz);
		yy[i] = 0.;
	    }
	}
    }
    
    for (iz=i1; iz < i2; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix = ap->nw*ap->nj; ix < nx-ap->nw*ap->nj; ix++) {
		i = ix + nx * (iy + ny * iz);
		
		if (ap->drift) {
		    id = SF_NINT(ap->pp[i]);
		    if (ix-ap->nw*ap->nj-id < 0 || 
			ix+ap->nw*ap->nj-id >= nx) continue;

		    if (der) {
			aderfilter(ap->pp[i]-id, ap->flt);
		    } else {
			passfilter(ap->pp[i]-id, ap->flt);
		    }
		    
		    for (iw = 0; iw <= 2*ap->nw; iw++) {
			is = (iw-ap->nw)*ap->nj;
			
			yy[i] += (xx[i+is+ip] - xx[i-is-id]) * ap->flt[iw];
		    }

		} else {
		    if (der) {
			aderfilter(ap->pp[i], ap->flt);
		    } else {
			passfilter(ap->pp[i], ap->flt);
		    }
		    
		    for (iw = 0; iw <= 2*ap->nw; iw++) {
			is = (iw-ap->nw)*ap->nj;
			
			yy[i] += (xx[i+is+ip] - xx[i-is]) * ap->flt[iw];
		    }
		}
	    }
	}
    }
}

void allpass3_init (allpass ap, allpass aq)
/*< Initialize linear operator >*/
{
    ap1 = ap;
    ap2 = aq;
}

void allpass3_lop (bool adj, bool add, int n1, int n2, float* xx, float* yy)
/*< PWD as linear operator >*/
{
    int i, ix, iy, iz, iw, is, nx, ny, nz, nw, nj, id;

    if (n2 != 2*n1) sf_error("%s: size mismatch: %d != 2*%d",__FILE__,n2,n1);

    sf_adjnull(adj, add, n1, n2, xx, yy);

    nx = ap1->nx;
    ny = ap1->ny;
    nz = ap1->nz;
    nw = ap1->nw;
    nj = ap1->nj;

    if (nx*ny*nz != n1) sf_error("%s: size mismatch",__FILE__);
    
    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny-1; iy++) {
	    for (ix = nw*nj; ix < nx-nw*nj; ix++) {
		i = ix + nx*(iy + ny*iz);

		if (ap1->drift) {
		    id = SF_NINT(ap1->pp[i]);
		    if (ix-nw*nj-id < 0 || 
			ix+nw*nj-id >= nx) continue;

		    passfilter(ap1->pp[i]-id, ap1->flt);
		    
		    for (iw = 0; iw <= 2*nw; iw++) {
			is = (iw-nw)*nj;
			
			if (adj) {
			    xx[i+nx+is] += yy[i] * ap1->flt[iw];
			    xx[i-is-id] -= yy[i] * ap1->flt[iw];
			} else {
			    yy[i] += (xx[i+nx+is] - xx[i-is-id]) * ap1->flt[iw];
			}
		    }
		} else {
		    passfilter(ap1->pp[i], ap1->flt);
		    
		    for (iw = 0; iw <= 2*nw; iw++) {
			is = (iw-nw)*nj;
			
			if (adj) {
			    xx[i+nx+is] += yy[i] * ap1->flt[iw];
			    xx[i-is]    -= yy[i] * ap1->flt[iw];
			} else {
			    yy[i] += (xx[i+nx+is] - xx[i-is]) * ap1->flt[iw];
			}
		    }
		}
	    }
	}
    }

    nx = ap2->nx;
    ny = ap2->ny;
    nz = ap2->nz;
    nw = ap2->nw;
    nj = ap2->nj;

    if (nx*ny*nz != n1) sf_error("%s: size mismatch",__FILE__);
    
    for (iz=0; iz < nz-1; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix = nw*nj; ix < nx-nw*nj; ix++) {
		i = ix + nx*(iy + ny*iz);

		if (ap2->drift) {
		    id = SF_NINT(ap2->pp[i]);
		    if (ix-nw*nj-id < 0 || 
			ix+nw*nj-id >= nx) continue;

		    passfilter(ap2->pp[i]-id, ap2->flt);
		    
		    for (iw = 0; iw <= 2*nw; iw++) {
			is = (iw-nw)*nj;
			
			if (adj) {
			    xx[i+nx*ny+is] += yy[i+n1] * ap2->flt[iw];
			    xx[i-is-id]    -= yy[i+n1] * ap2->flt[iw];
			} else {
			    yy[i+n1] += (xx[i+nx*ny+is] - xx[i-is-id]) * ap2->flt[iw];
			}
		    }
		} else {
		    passfilter(ap2->pp[i], ap2->flt);
		    
		    for (iw = 0; iw <= 2*nw; iw++) {
			is = (iw-nw)*nj;
			
			if (adj) {
			    xx[i+nx*ny+is] += yy[i+n1] * ap2->flt[iw];
			    xx[i-is]       -= yy[i+n1] * ap2->flt[iw];
			} else {
			    yy[i+n1] += (xx[i+nx*ny+is] - xx[i-is]) * ap2->flt[iw];
			}
		    }
		}
	    }
	}
    }
}





static float *u1, *u2, *dp, *p0, eps;
static int n, n1, n2, n3, nn[3];

void dip3_init(int m1, int m2, int m3 /* dimensions */, 
	       int* rect              /* smoothing radius [3] */, 	       
	       int niter              /* number of iterations */,
	       float eps1             /* regularization */,      
	       bool verb              /* verbosity flag */)
/*< initialize >*/
{
    n1=m1;
    n2=m2;
    n3=m3;
    n = n1*n2*n3;
    eps = eps1;

    u1 = sf_floatalloc(n);
    u2 = sf_floatalloc(n);
    dp = sf_floatalloc(n);
    p0 = sf_floatalloc(n);

    nn[0]=n1;
    nn[1]=n2;
    nn[2]=n3;

    sf_divn_init (3, n, nn, rect, niter, verb);
}

void dip3n_init(int m1, int m2, int m3 /* dimensions */, 
		   int *nbox 			  /* triangle radius [ndim] */, 
		   float **rct 			  /* triangle lengths [ndim][nd] */,
           int **sft 			  /* triangle shifts [ndim][nd] */,	       
	       int niter              /* number of iterations */,
	       float eps1             /* regularization */,      
	       bool verb              /* verbosity flag */)
/*< initialize (non-stationary) >*/
{
    n1=m1;
    n2=m2;
    n3=m3;
    n = n1*n2*n3;
    eps = eps1;

    u1 = sf_floatalloc(n);
    u2 = sf_floatalloc(n);
    dp = sf_floatalloc(n);
    p0 = sf_floatalloc(n);

    nn[0]=n1;
    nn[1]=n2;
    nn[2]=n3;

    divnn_init (3, n, nbox, nn, rct, sft, niter, verb);
}

void dip3_close(void)
/*< free allocated storage >*/
{
    free (u1);
    free (u2);
    free (dp);
    sf_divn_close();
}

void dip3(bool left               /* left or right prediction */,
	  int dip                 /* 1 - inline, 2 - crossline */, 
	  int niter               /* number of nonlinear iterations */, 
	  int nw                  /* filter size */, 
	  int nj                  /* filter stretch for aliasing */, 
	  bool drift              /* if shift filter */,
	  float *u                /* input data */, 
	  float* p                /* output dip */, 
	  bool* mask              /* input mask for known data */,
	  float pmin, float pmax  /* minimum and maximum dip */)
/*< estimate local dip >*/
{
    int i, iter, k;
    float usum, usum2, pi, lam;
    allpass ap;
 
    ap = allpass_init (nw,nj,n1,n2,n3,drift,p);

    if (dip == 1) {
	allpass1 (left, false, ap, u,u2);
    } else {
	allpass2 (left, false, ap, u,u2);
    }

    for (iter =0; iter < niter; iter++) {
	if (dip == 1) {
	    allpass1 (left, true,  ap, u,u1);
	} else {
	    allpass2 (left, true,  ap, u,u1);
	}

	usum = 0.0;
	for(i=0; i < n; i++) {
	    p0[i] = p[i];
	    usum += u2[i]*u2[i];
	}
	
	if (NULL != mask) {
	    for(i=0; i < n; i++) {
		if (mask[i]) {
		    u1[i] = 0.;
		    u2[i] = 0.;
		}
	    }
	}

	sf_divne (u2, u1, dp, eps);

	lam = 1.;
	for (k=0; k < 8; k++) {
	    for(i=0; i < n; i++) {
		pi = p0[i]+lam*dp[i];
		if (pi < pmin) pi=pmin;
		if (pi > pmax) pi=pmax;
		p[i] = pi;
	    }
	    if (dip == 1) {
		allpass1 (left, false, ap, u,u2);
	    } else {
		allpass2 (left, false, ap, u,u2);
	    }

	    usum2 = 0.;
	    for(i=0; i < n; i++) {
		usum2 += u2[i]*u2[i];
	    }
	    if (usum2 < usum) break;
	    lam *= 0.5;
	}
    } /* iter */

    allpass_close(ap);
}

void dip3n(bool left          /* left or right prediction */,
	  int dip                 /* 1 - inline, 2 - crossline */, 
	  int niter               /* number of nonlinear iterations */, 
	  int nw                  /* filter size */, 
	  int nj                  /* filter stretch for aliasing */, 
	  bool drift              /* if shift filter */,
	  float *u                /* input data */, 
	  float* p                /* output dip */, 
	  bool* mask              /* input mask for known data */,
	  float pmin, float pmax  /* minimum and maximum dip */)
/*< estimate local dip (non-stationary) >*/
{
    int i, iter, k;
    float usum, usum2, pi, lam;
    allpass ap;
 
    ap = allpass_init (nw,nj,n1,n2,n3,drift,p);

    if (dip == 1) {
	allpass1 (left, false, ap, u,u2);
    } else {
	allpass2 (left, false, ap, u,u2);
    }

    for (iter =0; iter < niter; iter++) {
	if (dip == 1) {
	    allpass1 (left, true,  ap, u,u1);
	} else {
	    allpass2 (left, true,  ap, u,u1);
	}

	usum = 0.0;
	for(i=0; i < n; i++) {
	    p0[i] = p[i];
	    usum += u2[i]*u2[i];
	}
	
	if (NULL != mask) {
	    for(i=0; i < n; i++) {
		if (mask[i]) {
		    u1[i] = 0.;
		    u2[i] = 0.;
		}
	    }
	}

	divnne (u2, u1, dp, eps);

	lam = 1.;
	for (k=0; k < 8; k++) {
	    for(i=0; i < n; i++) {
		pi = p0[i]+lam*dp[i];
		if (pi < pmin) pi=pmin;
		if (pi > pmax) pi=pmax;
		p[i] = pi;
	    }
	    if (dip == 1) {
		allpass1 (left, false, ap, u,u2);
	    } else {
		allpass2 (left, false, ap, u,u2);
	    }

	    usum2 = 0.;
	    for(i=0; i < n; i++) {
		usum2 += u2[i]*u2[i];
	    }
	    if (usum2 < usum) break;
	    lam *= 0.5;
	}
    } /* iter */

    allpass_close(ap);
}

void dip3n_fb(bool left          /* left or right prediction */,
	  int dip                 /* 1 - inline, 2 - crossline */, 
	  int niter               /* number of nonlinear iterations */, 
	  int nw                  /* filter size */, 
	  int nj                  /* filter stretch for aliasing */, 
	  bool drift              /* if shift filter */,
	  float *u                /* input data */, 
	  float* p                /* output dip */, 
	  bool* mask              /* input mask for known data */,
	  float pmin, float pmax  /* minimum and maximum dip */)
/*< estimate local dip (non-stationary) with forward and backward space derivative calculation>*/
{
    int i, iter, k;
    float usum, usum2, pi, lam;
    allpass ap;
 
    ap = allpass_init (nw,nj,n1,n2,n3,drift,p);

    if (dip == 1) {
	allpass1_fb (left, false, ap, u,u2);
    } else {
	allpass2_fb (left, false, ap, u,u2);
    }

    for (iter =0; iter < niter; iter++) {
	if (dip == 1) {
	    allpass1_fb (left, true,  ap, u,u1);
	} else {
	    allpass2_fb (left, true,  ap, u,u1);
	}

	usum = 0.0;
	for(i=0; i < n; i++) {
	    p0[i] = p[i];
	    usum += u2[i]*u2[i];
	}
	
	if (NULL != mask) {
	    for(i=0; i < n; i++) {
		if (mask[i]) {
		    u1[i] = 0.;
		    u2[i] = 0.;
		}
	    }
	}

	divnne (u2, u1, dp, eps);

	lam = 1.;
	for (k=0; k < 8; k++) {
	    for(i=0; i < n; i++) {
		pi = p0[i]+lam*dp[i];
		if (pi < pmin) pi=pmin;
		if (pi > pmax) pi=pmax;
		p[i] = pi;
	    }
	    if (dip == 1) {
		allpass1_fb (left, false, ap, u,u2);
	    } else {
		allpass2_fb (left, false, ap, u,u2);
	    }

	    usum2 = 0.;
	    for(i=0; i < n; i++) {
		usum2 += u2[i]*u2[i];
	    }
	    if (usum2 < usum) break;
	    lam *= 0.5;
	}
    } /* iter */

    allpass_close(ap);
}

void dip3_fb(bool left               /* left or right prediction */,
	  int dip                 /* 1 - inline, 2 - crossline */, 
	  int niter               /* number of nonlinear iterations */, 
	  int nw                  /* filter size */, 
	  int nj                  /* filter stretch for aliasing */, 
	  bool drift              /* if shift filter */,
	  float *u                /* input data */, 
	  float* p                /* output dip */, 
	  bool* mask              /* input mask for known data */,
	  float pmin, float pmax  /* minimum and maximum dip */)
/*< estimate local dip with forward and backward space derivative calculation >*/
{
    int i, iter, k;
    float usum, usum2, pi, lam;
    allpass ap;
 
    ap = allpass_init (nw,nj,n1,n2,n3,drift,p);

    if (dip == 1) {
	allpass1_fb (left, false, ap, u,u2);
    } else {
	allpass2_fb (left, false, ap, u,u2);
    }

    for (iter =0; iter < niter; iter++) {
	if (dip == 1) {
	    allpass1_fb (left, true,  ap, u,u1);
	} else {
	    allpass2_fb (left, true,  ap, u,u1);
	}

	usum = 0.0;
	for(i=0; i < n; i++) {
	    p0[i] = p[i];
	    usum += u2[i]*u2[i];
	}
	
	if (NULL != mask) {
	    for(i=0; i < n; i++) {
		if (mask[i]) {
		    u1[i] = 0.;
		    u2[i] = 0.;
		}
	    }
	}

	sf_divne (u2, u1, dp, eps);

	lam = 1.;
	for (k=0; k < 8; k++) {
	    for(i=0; i < n; i++) {
		pi = p0[i]+lam*dp[i];
		if (pi < pmin) pi=pmin;
		if (pi > pmax) pi=pmax;
		p[i] = pi;
	    }
	    if (dip == 1) {
		allpass1_fb (left, false, ap, u,u2);
	    } else {
		allpass2_fb (left, false, ap, u,u2);
	    }

	    usum2 = 0.;
	    for(i=0; i < n; i++) {
		usum2 += u2[i]*u2[i];
	    }
	    if (usum2 < usum) break;
	    lam *= 0.5;
	}
    } /* iter */

    allpass_close(ap);
}


void mask32 (bool both              /* left and right predictions */,
	     int nw                 /* filter size */, 
	     int nj1, int nj2       /* dealiasing stretch */, 
	     int nx, int ny, int nz /* data size */, 
	     float *yy              /* data [nz*ny*nx] */, 
	     bool **m               /* dip mask [both? 4:2][nz*ny*nx] */) 
/*< two-dip masks in 3-D >*/
{
    int ix, iy, iz, iw, is, i, n;
    bool *xx;

    n = nx*ny*nz;

    xx = sf_boolalloc(n);

    for (i=0; i < n; i++) {
	xx[i] = (bool) (yy[i] == 0.);
	m[0][i] = false;
	m[1][i] = false;
	if (both) {
	    m[2][i] = false;
	    m[3][i] = false;
	}
    }

    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny-1; iy++) {
	    for (ix = nw*nj1; ix < nx-nw*nj1; ix++) {
		i = ix + nx * (iy + ny * iz);

		for (iw = 0; iw <= 2*nw; iw++) {
		    is = (iw-nw)*nj1;		  
		    m[0][i] = (bool) (m[0][i] || xx[i-is] || xx[i+nx+is]);
		}
	    }
	}
    }

    for (iz=0; iz < nz-1; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix = nw*nj2; ix < nx-nw*nj2; ix++) {
		i = ix + nx * (iy + ny * iz);

		for (iw = 0; iw <= 2*nw; iw++) {
		    is = (iw-nw)*nj2;		  
		    m[1][i] = (bool) (m[1][i] || xx[i-is] || xx[i+ny*nx+is]);
		}
	    }
	}
    }

    if (!both) {
	free(xx);
	return;
    }

    for (iz=0; iz < nz; iz++) {
	for (iy=1; iy < ny; iy++) {
	    for (ix = nw*nj1; ix < nx-nw*nj1; ix++) {
		i = ix + nx * (iy + ny * iz);

		for (iw = 0; iw <= 2*nw; iw++) {
		    is = (iw-nw)*nj1;		  
		    m[2][i] = (bool) (m[2][i] || xx[i-is] || xx[i-nx+is]);
		}
	    }
	}
    }

    for (iz=1; iz < nz; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix = nw*nj2; ix < nx-nw*nj2; ix++) {
		i = ix + nx * (iy + ny * iz);

		for (iw = 0; iw <= 2*nw; iw++) {
		    is = (iw-nw)*nj2;		  
		    m[3][i] = (bool) (m[3][i] || xx[i-is] || xx[i-ny*nx+is]);
		}
	    }
	}
    }
	
    free(xx); 
}

void mask3 (int nw         /* filter size */, 
	    int nj         /* dealiasing stretch */, 
	    int nx, int ny /* data size */, 
	    float **yy     /* data */, 
	    bool **mm      /* mask */) 
/*< one-dip mask in 2-D >*/
{
    int ix, iy, iw, is;
    bool **xx;

    xx = sf_boolalloc2(nx,ny);
    
    for (iy=0; iy < ny; iy++) {
	for (ix=0; ix < nx; ix++) {
	    xx[iy][ix] = (bool) (yy[iy][ix] == 0.);
	    mm[iy][ix] = false;
	}
    }

    for (iy=0; iy < ny-1; iy++) {
	for (ix = nw*nj; ix < nx-nw*nj; ix++) {
	    for (iw = 0; iw <= 2*nw; iw++) {
		is = (iw-nw)*nj;
		mm[iy][ix] = (bool) (mm[iy][ix] || xx[iy+1][ix+is] || xx[iy][ix-is]);
	    }
	}
    }
    
    free(xx[0]);
    free(xx);
}

void mask6 (int nw           /* filter size */, 
	    int nj1, int nj2 /* dealiasing stretch */, 
	    int nx, int ny   /* data size */, 
	    float *yy       /* data [ny][nx] */, 
	    bool *mm        /* mask [ny][nx] */) 
/*< two-dip mask in 2-D >*/
{
    int ix, iy, iw, is, n, i;
    bool *xx;

    n = nx*ny;

    xx = sf_boolalloc(n);
    
    for (i=0; i < n; i++) {
	mm[i] = (bool) (yy[i] == 0.);
	xx[i] = false;
    }

    for (iy=0; iy < ny-1; iy++) {
	for (ix = nw*nj1; ix < nx-nw*nj1; ix++) {
	    i = ix + nx*iy;

	    for (iw = 0; iw <= 2*nw; iw++) {
		is = (iw-nw)*nj1;
		xx[i] = (bool) (xx[i] || mm[i+nx+is] || mm[i-is]);
	    }
	}
    }
    
    for (i=0; i < n; i++) {
	mm[i] = false;
    }
    
    for (iy=0; iy < ny-1; iy++) {
	for (ix = nw*nj2; ix < nx-nw*nj2; ix++) {
	    i = ix + nx*iy;
	    
	    for (iw = 0; iw <= 2*nw; iw++) {
		is = (iw-nw)*nj2;
		mm[i] = (bool) (mm[i] || xx[i+nx+is] || xx[i-is]);
	    }
	}
    }
    
    free(xx);
}
