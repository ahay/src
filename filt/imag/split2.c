/* 3-D extended split-step migration/modeling */
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
/*^*/

#include "split2.h"
#include "fft2.h"
#include "taper.h"
#include "slowref.h"

#include "slice.h"
/*^*/

#define LOOPxy(a) for(ix=0;ix<nx;ix++){ for(iy=0;iy<ny;iy++){ {a} }}
#define LOOPxy2(a) for(ix=0;ix<nx2;ix++){ for(iy=0;iy<ny2;iy++){ {a} }}

static int nx,ny,nx2,ny2,nz,nrmax;
static float dz, ds, ds2;
static float         **ss; /* slowness */
static float         **sm; /* reference slowness squared */
static int            *nr; /* number of references */

static float         **qq; /* image */
static float         **kk; /* wavenumber */

static float complex **pp; /* wavefield */
static float complex **wx; /* x wavefield */
static float complex **wk; /* k wavefield */

static float         **wt; /* interpolation weight */
static float         **ph; /* phase */

void split2_init(int nz1, float dz1  /* depth */,
		 int ny1, float dy1  /* in-line midpoint */,
		 int nx1, float dx1  /* cross-line midpoint */,
		 int ntx, int nty    /* taper size */,
		 int padx, int pady  /* padding */,
		 int nr1             /* maximum number of references */,
		 float dt            /* time error */) 
/*< initialize >*/
{
    int ix, iy, jx, jy;
    float x0, y0, dx, dy, kx, ky;

    nz = nz1;
    dz = dz1;

    ds = dt/dz;
    ds2 = ds*ds;
    ds2 *= ds2;

    nx = nx1;
    nx2 = nx+padx;
    dx =        2.0*SF_PI/(nx2*dx1);
    x0 = (1==nx2)?0:-SF_PI/dx1;

    ny = ny1;
    ny2 = ny+pady;
    dy =        2.0*SF_PI/(ny2*dy1);
    y0 = (1==ny2)?0:-SF_PI/dy1;

    nrmax = nr1;

    fft2_init(ny2,nx2);

    /* allocate workspace */
    ss = sf_floatalloc2 (ny,nx);    /* slowness */
    sm = sf_floatalloc2 (nrmax,nz); /* reference slowness squared*/
    nr = sf_intalloc    (      nz); /* number of reference slownesses */

    qq = sf_floatalloc2   (ny,nx);  /* image */
    kk = sf_floatalloc2   (ny2,nx2);  /* wavenumber */
 
    pp = sf_complexalloc2 (ny2,nx2);  /* wavefield */
    wx = sf_complexalloc2 (ny,nx);    /* x wavefield */
    wk = sf_complexalloc2 (ny2,nx2);  /* k wavefield */

    wt = sf_floatalloc2   (ny,nx);
    ph = sf_floatalloc2   (ny2,nx2);

    /* precompute wavenumbers */
    for (ix=0; ix<nx2; ix++) {
	jx = (ix < nx2/2)? ix + nx2/2: ix - nx2/2;
	kx = x0 + jx*dx;
	kx *= kx;
	
	for (iy=0; iy<ny2; iy++) {
	    jy = (iy < ny2/2)? iy + ny2/2: iy - ny2/2;
	    ky = y0 + jy*dy;
	    kk[ix][iy] = ky*ky + kx;
	}
    }    

    /* tapering */
    taper2_init(nx,ny,ntx,nty,true,true);
}

void split2_close(void)
/*< free allocated storage >*/
{
    fft2_close();

    free(*pp); free(pp);
    free(*wx); free(wx);
    free(*wk); free(wk);

    free(*ss); free(ss);
    free(*sm); free(sm);
    free( nr);

    free(*qq); free(qq);
    free(*kk); free(kk);

    free(*wt); free(wt);

    taper2_close();
    fft2_close();
}

void split2(bool verb                   /* verbosity flag */, 
	    bool inv                    /* migration/modeling flag */, 
	    float eps                   /* stability factor */,  
	    int nw, float dw, float w0  /* frequency (radian) */,
	    sf_file data                /* data  [nw][nx][ny] */,
	    slice imag                  /* image  [nz][nx][ny] */,
	    slice slow                  /* slowness  [nz][nx][ny] */)
/*< Apply migration/modeling >*/
{
    int iz,iw,ix,iy,ir,nxy;
    float complex cshift, cref, w, w2;
    float d;

    if (!inv) { /* prepare image for migration */
	LOOPxy( qq[ix][iy] = 0.0; );
	for (iz=0; iz<nz; iz++) {
	    slice_put(imag,iz,qq[0]);
	}
    }

    /* compute reference slowness */
    nxy = nx*ny;
    for (iz=0; iz<nz; iz++) {
	slice_get(slow,iz,  ss[0]);
	nr[iz] = slowref(nrmax,ds,nxy,ss[0],sm[iz]);
	if (verb) sf_warning("nr[%d]=%d",iz,nr[iz]);
    }
    for (iz=0; iz<nz-1; iz++) {
	for (ir=0; ir<nr[iz]; ir++) {
	    sm[iz][ir] = 2*(sm[iz][ir]+sm[iz+1][ir]);
	}
    }

    /* loop over frequencies w */
    for (iw=0; iw<nw; iw++) {
	if (verb) sf_warning ("frequency %d of %d",iw+1,nw);

	w = eps*dw + I*(w0+iw*dw);
	w2 = w*w;

	if (inv) { /* MODELING */
	    /* start from bottom */
	    slice_get(slow,nz-1,ss[0]);
	    
	    /* imaging condition */
	    slice_get(imag,nz-1,qq[0]);
	    LOOPxy( wx[ix][iy] = qq[ix][iy]; );

	    /* loop over migrated depths z */
	    for (iz=nz-2; iz>=0; iz--) {
		
		/* w-x @ bottom */
		LOOPxy2( pp[ix][iy] = 0.;);
		taper2(wx);
		LOOPxy( cshift = cexpf(-w*ss[ix][iy]*dz);
			pp[ix][iy] = wx[ix][iy]*cshift; 
			wt[ix][iy] = 0.; 
			wx[ix][iy] = 0.; 
		    );
		
		fft2(false,pp);
		slice_get(slow,iz,ss[0]);

		for (ir=0; ir<nr[iz]; ir++) {
		    /* w-k phase shift */
		    cref = csqrtf(w2*sm[iz][ir]);
		    LOOPxy2( cshift = 
			     cexpf((cref- 
				    csqrtf(w2*sm[iz][ir]+kk[ix][iy]))*dz);
			     wk[ix][iy] = pp[ix][iy]*cshift; );
		    fft2(true,wk);

		    /* accumulate wavefield */
		    LOOPxy( d = fabsf(4.*ss[ix][iy]*ss[ix][iy]-sm[iz][ir]);
			    d = ds2/(d*d+ds2);
			    wx[ix][iy] += wk[ix][iy]*d;
			    wt[ix][iy] += d;
			);
		} /* ir loop */

		slice_get(imag,iz,qq[0]);
		
		/* w-x @ top */
		LOOPxy( cshift = cexpf(-w*ss[ix][iy]*dz)/wt[ix][iy];
			wx[ix][iy] = qq[ix][iy] + wx[ix][iy]*cshift; );
	    } /* iz */

	    taper2(wx);
	    sf_complexwrite(wx[0],nxy,data);

	} else { /* MIGRATION */

	    sf_complexread(wx[0],nxy,data);
	    taper2(wx);

	    slice_get(slow,0,ss[0]);

	    /* loop over migrated depths z */
	    for (iz=0; iz< nz-1; iz++) {

		/* imaging condition */
		slice_get(imag,iz,qq[0]);
		LOOPxy( qq[ix][iy] += crealf(wx[ix][iy]); );
		slice_put(imag,iz,qq[0]);

		/* w-x @ top */
		LOOPxy2( pp[ix][iy]=0.;);
		LOOPxy( cshift= conjf(cexpf(-w*ss[ix][iy]*dz));
			pp[ix][iy] = wx[ix][iy]*cshift; 
			wt[ix][iy] = 0.; 
			wx[ix][iy] = 0.; );
		
		fft2(false,pp);
		slice_get(slow,iz+1,ss[0]);

		for (ir=0; ir<nr[iz]; ir++) {
		    /* w-k phase shift */
		    cref = csqrtf(w2*sm[iz][ir]);
		    LOOPxy2( cshift = 
			     cexpf((cref-csqrtf(w2*sm[iz][ir]+kk[ix][iy]))*dz);
			     wk[ix][iy] = pp[ix][iy]*conjf(cshift); );
		    fft2(true,wk);

		    /* accumulate wavefield */
		    LOOPxy( d = fabsf(4.*ss[ix][iy]*ss[ix][iy]-sm[iz][ir]);
			    d = ds2/(d*d+ds2);
			    wx[ix][iy] += wk[ix][iy]*d;
			    wt[ix][iy] += d; );
		} /* ir loop */
		
		/* w-x @ bottom */
		LOOPxy( 
		    cshift = conjf(cexpf(-w*ss[ix][iy]*dz))/wt[ix][iy];
		    wx[ix][iy] *= cshift; 
		    );
		taper2(wx);
	    } /* iz */

	    /* imaging condition @ bottom */ 
	    slice_get(imag,nz-1,qq[0]);
	    LOOPxy( qq[ix][iy] += crealf(wx[ix][iy]); );
	    slice_put(imag,nz-1,qq[0]);
	    
	} /* else */
    } /* iw */
}

/*
static float complex muir(int niter,float complex s,float complex q,
			  float k2,float dz)
{
    float complex q2;
    int iter;
    
    for (iter=0; iter < niter; iter++) {
	q2 = k2/(q + 2.*s);
	if (cabsf(q2-q) < FLT_EPSILON) break;
	q = q2;
    }

    return cexpf(-q2*dz);
}
*/
