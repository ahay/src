/* 3-D split-step migration/modeling */
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

#define LOOPxy(a) for(ix=0;ix<nx;ix++){ for(iy=0;iy<ny;iy++){ {a} }}
/*^*/

#include "slice.h"
/*^*/

static int nx,ny,nz,nt,nr;
static float dz;
static float         **ss; /* slowness */
static float         **so; /* slowness */
static float         **qq; /* image */
static float         **kk; /* wavenumber */
static float         **tt; /* taper */
static float         **sz; /* reference slowness */

static int           **mm; /* multi-reference slowness map  */
static float         **ma; /* multi-reference slowness mask */

static float complex **pp; /* wavefield */
static float complex **wt; /* wavefield top */
static float complex **wb; /* wavefield bot */

void split2_init(int nz1, float dz1  /* depth */,
		 int ny1, float dy1  /* in-line midpoint */,
		 int nx1, float dx1  /* cross-line midpoint */,
		 int nt1             /* taper size */, 
		 int nr1             /* number of reference slownesses */ )
/*< initialize >*/
{
    int ix, iy, jx, jy;
    float x0, y0, dx, dy, kx, ky;

    nz = nz1;
    dz = dz1;

    nx = nx1;
    dx = 2.0*SF_PI/(nx*dx1);
    x0 = (1==nx)?0:-SF_PI/dx1;

    ny = ny1;
    dy = 2.0*SF_PI/(ny*dy1);
    y0 = (1==ny)?0:-SF_PI/dy1;

    nt = nt1;
    nr = nr1;

    fft2_init(ny,nx);

    /* allocate workspace */
    ss  = sf_floatalloc2 (ny,nx);  /* slowness */
    so  = sf_floatalloc2 (ny,nx);  /* slowness */
    sz  = sf_floatalloc2 (nr,nz);  /* reference slowness */

    qq  = sf_floatalloc2 (ny,nx);  /* image */
    kk  = sf_floatalloc2 (ny,nx);  /* wavenumber */
    tt  = sf_floatalloc2 (ny,nx);  /* taper */

    mm  = sf_intalloc2   (ny,nx);  /* MRS map */
    ma  = sf_floatalloc2 (ny,nx);  /* MRS mask */

    pp  = sf_complexalloc2 (ny,nx);/* wavefield */
    wt  = sf_complexalloc2 (ny,nx);/* wavefield top */
    wb  = sf_complexalloc2 (ny,nx);/* wavefield bot */

    /* precompute wavenumbers */
    for (ix=0; ix<nx; ix++) {
	jx = (ix < nx/2)? ix + nx/2: ix - nx/2;
	kx = x0 + jx*dx;
	kx *= kx;
	
	for (iy=0; iy<ny; iy++) {
	    jy = (iy < ny/2)? iy + ny/2: iy - ny/2;
	    ky = y0 + jy*dy;
	    kk[ix][iy] = ky*ky + kx;
	}
    }    

    /* precompute taper array */
    taper_init(nt);
}

void split2_close(void)
/*< free allocated storage >*/
{
    free(*pp); free( pp);
    free(*wt); free( wt);
    free(*wb); free( wb);

    free(*ss); free( ss);
    free(*so); free( so);
    free(*sz); free( sz);

    free(*qq); free( qq);
    free(*kk); free( kk);
    free(*tt); free( tt);

    free(*mm); free( mm);
    free(*ma); free( ma);
}

void split2(bool verb                   /* verbosity flag */, 
	    bool inv                    /* migration/modeling flag */, 
	    float eps                   /* stability factor */,  
	    int nw, float dw, float w0  /* frequency (radian) */,
	    float complex *** cp        /* data  [nw][nx][ny] */,
	    slice imag                  /* imag  [nz][nx][ny] */,
	    slice slow                  /* slow  [nz][nx][ny] */)
/*< Apply migration/modeling >*/
{
    int iz,iw,ix,iy,ir;
    float complex cshift, cref, w, w2, **pp;
    float qr;

    if (!inv) { /* prepare image for migration */
	LOOPxy( qq[ix][iy] = 0.0; );
	for (iz=0; iz<nz; iz++) {
	    slice_put(imag,iz,qq[0]);
	}
    }

    /* find reference from average slowness in layer */
    for (iz=0; iz<nz-1; iz++) {
	slice_get(slow,iz,  so[0]);
	slice_get(slow,iz+1,ss[0]);
	LOOPxy( ss[ix][iy] += so[ix][iy]; 
		ss[ix][iy] /= 2; );
	for (ir=0; ir<nr; ir++) {
	    qr = (ir+1.0)/nr - 0.5 * 1./nr;
	    sz[iz][ir]  = sf_quantile(qr*nx*ny,nx*ny,ss[0]);
	}
    }
    for (ir=0; ir<nr; ir++) {
	sz[nz-1][ir] = sz[nz-2][ir];
    }
    
    /* compute reference slowness squared */
    for (iz=0; iz<nz; iz++) {
	sz[iz][ir] *= sz[iz][ir];
    }

    /* loop over frequencies w */
    for (iw=0; iw<nw; iw++) {
	if (verb) sf_warning ("frequency %d of %d",iw+1,nw);

	w = eps*dw + I*(w0+iw*dw);
	w2 = w*w;

	pp = cp[iw];

	if (inv) { /* MODELING */
	    /* start from bottom */
	    slice_get(slow,nz-1,ss[0]);
	    
	    /* imaging condition */
	    slice_get(imag,nz-1,qq[0]);
	    LOOPxy( pp[ix][iy] = qq[ix][iy]; );

	    /* loop over migrated depths z */
	    for (iz=nz-2; iz>=0; iz--) {
		
		/* w-x @ bottom */
		LOOPxy( cshift = cexpf(-0.5*w*ss[ix][iy]*dz);
			pp[ix][iy] *= cshift; );
		
		/* taper and FFT */
		LOOPxy( pp[ix][iy] *= tt[ix][iy]; );
		fft2(false,pp);
		
		/* w-k phase shift */
		cref = csqrtf(w2*sz[iz][ir]);
		LOOPxy( cshift = cexpf((cref-csqrtf(w2*sz[iz][ir]+kk[ix][iy]))*dz);
			pp[ix][iy] *= cshift; );


		/* taper and IFT */
		fft2(true,pp);
		LOOPxy( pp[ix][iy] *= tt[ix][iy]; );
		
		/* w-x at top */
		slice_get(slow,iz,ss[0]);		
		LOOPxy( cshift = cexpf(-0.5*w*ss[ix][iy]*dz);
			pp[ix][iy] *= cshift; );
		
		/* imaging condition */
		slice_get(imag,iz,qq[0]);
		LOOPxy( pp[ix][iy] += qq[ix][iy]; );

	    } /* iz */
	} else { /* MIGRATION */
	    slice_get(slow,0,ss[0]);

	    /* loop over migrated depths z */
	    for (iz=0; iz< nz-1; iz++) {

		/* imaging condition */
		slice_get(imag,iz,qq[0]);
		LOOPxy( qq[ix][iy] += crealf(pp[ix][iy]); );
		slice_put(imag,iz,qq[0]);

		/* w-x @ top */
		LOOPxy( cshift= cexpf(-0.5*w*ss[ix][iy]*dz);
			pp[ix][iy] *= conjf(cshift); );

		/* taper and FFT */
		LOOPxy( pp[ix][iy] *= tt[ix][iy]; );
		fft2(false,pp);

/* START MRS */
		/* create MRS map */
		slice_get(slow,iz  ,so[0]);
		slice_get(slow,iz+1,ss[0]);
		LOOPxy( ss[ix][iy] += so[ix][iy]; 
			ss[ix][iy] /= 2; );

		LOOPxy( mm[ix][iy] = 0; );
		for (ir=0; ir<nr-1; ir++) {
		    LOOPxy( if(ss[ix][iy] > 0.5*(sz[iz][ir]+sz[iz][ir+1] )) mm[ix][iy] +=1; );
		}

		LOOPxy( wt[ix][iy] = pp[ix][iy]; ); /* save wavefield */
		LOOPxy( wb[ix][iy] = 0; );
		for (ir=0; ir<nr; ir++) {
		    LOOPxy( pp[ix][iy] = wt[ix][iy]; );
		    
		    /* w-k phase shift */
		    cref = csqrtf(w2*sz[iz][ir]);
		    LOOPxy( cshift = cexpf((cref-csqrtf(w2*sz[iz][ir]+kk[ix][iy]))*dz);
			    pp[ix][iy] *= conjf(cshift); );
		    
		    /* IFT and taper */
		    fft2(true,pp);
		    LOOPxy( pp[ix][iy] *= tt[ix][iy]; );

		    /* create MRS mask */
		    LOOPxy(                     ma[ix][iy]=0.; );
		    LOOPxy( if( mm[ix][iy]==ir) ma[ix][iy]=1.; );

		    /* accumulate wavefield */
		    LOOPxy( wb[ix][iy] += pp[ix][iy] * ma[ix][iy]; );
		    
		} /* ir loop */
		LOOPxy( pp[ix][iy] = wb[ix][iy]; ); /* collect wavefield */
/* END MRS */
		
		/* w-x @ bottom */
		slice_get(slow,iz+1,ss[0]);
		LOOPxy( cshift = cexpf(-0.5*w*ss[ix][iy]*dz);
			pp[ix][iy] *= conjf(cshift); );
	    } /* iz */

	    /* imaging condition @ bottom */ 
	    slice_get(imag,nz-1,qq[0]);
	    LOOPxy( qq[ix][iy] += crealf(pp[ix][iy]); );
	    slice_put(imag,nz-1,qq[0]);

	} /* else */
    } /* iw */
}

void taper_init(int nt)
/*< Initialize boundary taper >*/
{
    int it,ix,iy;

    LOOPxy( tt[ix][iy] = 1; );
    
    if(nt>=1) {
	
	    if(nx>=nt) {
		for(it=0; it<nt; it++) {
		    for(iy=0; iy<ny; iy++) {
			tt[   it  ][iy] = cos(SF_PI/2* (float)SF_ABS(nt-it)/nt);
			tt[nx-it-1][iy] = cos(SF_PI/2* (float)SF_ABS(nt-it)/nt);
		    }
		}
	    }
	    
	    if(ny>=nt) {
		for(it=0; it<nt; it++) {
		    for(ix=0; ix<nx; ix++) {
			tt[ix][   it  ] = cos(SF_PI/2* (float)SF_ABS(nt-it)/nt);
			tt[ix][ny-it-1] = cos(SF_PI/2* (float)SF_ABS(nt-it)/nt);
		    }
		}
	    }

    } /* nt>1 */

}
