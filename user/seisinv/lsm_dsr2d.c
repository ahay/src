/* 2-D prestack migration/modeling by split-step DSR */
/*
  Copyright (C) 2012 China University of Petroleum
  Copyright (C) 2006 Colorado School of Mines
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

#include "lsm_dsr2d.h"
#include "fft2.h"
#include "taper.h"
#include "slowref.h"

#define LOOPxh(a)  for(ix=0;ix<nx;ix++){ for(ih=0;ih<nh; ih++){ {a} }}
#define LOOPxh2(a) for(ix=0;ix<nx;ix++){ for(ih=0;ih<nh2;ih++){ {a} }}
#define LOOPuh(a)  for(iu=0;iu<nu;iu++){ for(ih=0;ih<nh; ih++){ {a} }}

static int nx,nh,nh2,nz,nu,nrmax;
static float dz;

static sf_complex         **qq;             /* image */
static int     **is, **ir, *ii;        /* indices */
static float   **ks, **kr;             /* wavenumber */

static float         **sz;             /* reference slowness */
static float         **sm;             /* reference slowness squared */
static int            *nr;             /* number of references */

static sf_complex **pk;   /* wavefield */
static sf_complex **wk;   /* wavefield k */
static sf_complex **wx;   /* wavefield x */

static int           **ms, **mr; /* multi-reference slowness map  */
static bool          ***skip;
static float         **ma; /* multi-reference slowness mask */
static sf_fslice mms, mmr;
static float dw,w0;
static int nw;
static float **slow;
static float dt;    /*time error */
static bool verb;   /* verbosiy falg */
static float eps;   /*stability flag */




void lsm_dsr2_init(int nz1, float dz1             /* depth */,
		   int nh1, float dh1, float h01  /* half-offset */,
		   int nx1, float dx1, float x01  /* midpoint */,
		   int nu1, float du,  float u0   /* slowness grid */,
		   int nw1, float dw1, float w01,
		   int ntx, int nth               /* taper size */,
		   int nr1                        /* number of references */,
		   int npad                       /* padding on nh */,
		   bool verb1, float eps1,
		   float **slow1,float dt1)
/*< initialize >*/
{
    int   ix,  ih,  iu;
    int   jx,  jh;
    float  x,   h;
    float  kx, kh, k;
    float  dx, dh;
    float   x0, h0;

    slow=slow1;
    verb=verb1;
    eps=eps1;
    dt=dt1;
    
    nz = nz1;
    dz = dz1;

    nx = nx1;
    dx = 2.0*SF_PI/(nx*dx1);
    x0 =    -SF_PI/    dx1 ;

    nh = nh1;
    nh2 = nh+npad;

    dh = 2.0*SF_PI/(nh2*dh1);
    h0 =    -SF_PI/     dh1 ;

    nu = nu1;

    nw=nw1;
    dw=dw1;
    w0=w01;

    nrmax = nr1;

    fft2_init(nh2,nx);

    /* allocate workspace */
    sz = sf_floatalloc2  (nrmax,nz);       /* reference slowness */
    sm = sf_floatalloc2  (nrmax,nz);       /* reference slowness squared*/
    nr = sf_intalloc     (      nz);       /* number of reference slownesses */

    qq = sf_complexalloc2     (nh,nu);       /* image */

    ks = sf_floatalloc2   (nh2,nx);        /* source wavenumber */
    kr = sf_floatalloc2   (nh2,nx);        /* receiver wavenumber */
    is = sf_intalloc2     (nh, nx);        /* source reference */
    ir = sf_intalloc2     (nh, nx);        /* receiver reference */
    ii = sf_intalloc          (nx);        /* midpoint reference */

    pk = sf_complexalloc2 (nh2,nx);        /* padded wavefield */ 
    wx = sf_complexalloc2 (nh, nx);        /* x wavefield */
    wk = sf_complexalloc2 (nh2,nx);        /* k wavefield */

    ms = sf_intalloc2     (nh, nx);        /* MRS map source */
    mr = sf_intalloc2     (nh, nx);        /* MRS map receiver */
    ma = sf_floatalloc2   (nh, nx);        /* MRS mask */

    skip = sf_boolalloc3 (nrmax,nrmax,nz); /* skip slowness combination */

    /* precompute wavenumbers */
    for (ix=0; ix<nx; ix++) {
	jx = (ix < nx/2)? ix + nx/2: ix - nx/2;
	kx = x0 + jx*dx;
	x = x01 + ix*dx1;
	
	iu = 0.5+(x-u0)/du;
	if      (iu <   0) iu=0;
	else if (iu >= nu) iu=nu-1;
	ii[ix] = iu;

	for (ih=0; ih<nh; ih++) {
	    h = h01 + ih*dh1;
	    
	    iu = 0.5+(x-h-u0)/du;
	    if      (iu <   0) iu=0;
	    else if (iu >= nu) iu=nu-1;
	    is[ix][ih] = iu;

	    iu = 0.5+(x+h-u0)/du;
	    if      (iu <   0) iu=0;
	    else if (iu >= nu) iu=nu-1;
	    ir[ix][ih] = iu;
	}

	for (ih=0; ih<nh2; ih++) {
	    jh = (ih < nh2/2)? ih + nh2/2: ih - nh2/2;
	    kh = h0 + jh*dh;

	    k = 0.5*(kx-kh);
	    ks[ix][ih] = k*k;

	    k = 0.5*(kx+kh);
	    kr[ix][ih] = k*k;
	}
    }

    /* precompute taper array */
    taper2_init(nx,nh,ntx,nth,true,false);

    mms = sf_fslice_init(nh*nx,nz,sizeof(int));
    mmr = sf_fslice_init(nh*nx,nz,sizeof(int));
}




void lsm_dsr2_close(void)
/*< free allocated storage >*/
{
    free( *pk); free( pk);
    free( *wk); free( wk);
    free( *wx); free( wx);

    free( *sz); free( sz);
    free( *sm); free( sm);
    free(  nr);

    free( *qq); free( qq);
    free( *ks); free( ks);
    free( *kr); free( kr);
    free( *is); free( is);
    free( *ir); free( ir);
    free(  ii);

    free( *ms); free( ms);
    free( *mr); free( mr);
    free( *ma); free( ma);
    
    free(**skip); 
    free( *skip); free(skip);
    sf_fslice_close(mms);
    sf_fslice_close(mmr);
    taper2_close();
    fft2_close();
}




void lsm_dsr2_lop(bool adj            /* if y,do migration */, 
		  bool add,
		  int m,int n,
		  sf_complex *img     /* model vector,number=[nz][nu][nh] */,
		  sf_complex *dat     /* data vector,number=[nw][nx][nh] */)
/*< Apply migration/modeling >*/
{
    int iz,iw,ix,ih,iu,j,k;
    float sy, *si;
    sf_complex cshift, w, w2, cs, cr, cref;

    sf_cadjnull(adj,false,m,n,img,dat);

    /* compute reference slowness */
    for (iz=0; iz<nz; iz++) {
	si = slow[iz];
	nr[iz] = slowref(nrmax,dt/dz,nu,si,sm[iz]);

/*	nr[iz] = slowref(nrmax,ds,alx.n*aly.n,ss[0],sm[iz]);*/

	if (verb) sf_warning("nr[%d]=%d",iz,nr[iz]);
	/* create MRS map */
	LOOPxh( ms[ix][ih] = 0; mr[ix][ih] = 0.;);

	for (j=0; j<nr[iz]; j++) {
            sz[iz][j] = sqrtf(sm[iz][j]);
        }
	for (j=0; j<nr[iz]-1; j++) {
	    sy = 0.5*(sz[iz][j]+sz[iz][j+1]);
	    LOOPxh( if(si[ is[ix][ih] ] > sy) ms[ix][ih]++;
		    if(si[ ir[ix][ih] ] > sy) mr[ix][ih]++; );
	}
	sf_fslice_put(mms,iz,ms[0]);
	sf_fslice_put(mmr,iz,mr[0]);
	for (j=0; j<nr[iz]; j++) {
	    for (k=0; k<nr[iz]; k++) {
		skip[iz][j][k] = true;
	    }
	}
	LOOPxh( skip[iz][ms[ix][ih]][mr[ix][ih]] = false; );
    }

    for (iz=0; iz<nz-1; iz++) {
	for (j=0; j<nr[iz]; j++) {
	    sm[iz][j] = 0.5*(sm[iz][j]+sm[iz+1][j]);
	}
    }
    
    /* loop over frequencies w */
    for (iw=0; iw<nw; iw++) {
	if (verb) sf_warning ("frequency %d of %d",iw+1,nw);

	w = sf_cmplx(eps*dw,w0+iw*dw);
#ifdef SF_HAS_COMPLEX_H
	w2 = w*w;
#else
	w2 = sf_cmul(w,w);
#endif

	if (!adj) { /* MODELING */
	    /* start from bottom */
	    si = slow[nz-1];

	    /* imaging condition */
	    LOOPuh(qq[iu][ih]=img[(nz-1)*nu*nh+iu*nh+ih];);
		
	    LOOPxh( wx[ix][ih] =qq[ ii[ix] ][ih];  );

	    /* loop over migrated depths z */
	    for (iz=nz-2; iz>=0; iz--) {
		/* w-x @ bottom */
		LOOPxh2( pk[ix][ih] = sf_cmplx(0.,0.););
		taper2(wx);
#ifdef SF_HAS_COMPLEX_H
		LOOPxh( sy = 0.5*(si[ is[ix][ih] ] + 
				  si[ ir[ix][ih] ]);
			cshift = cexpf(-w*sy*dz);
			pk[ix][ih] = 
			wx[ix][ih] * cshift; );
#else
		LOOPxh( sy = 0.5*(si[ is[ix][ih] ] + 
				  si[ ir[ix][ih] ]);
			cshift = cexpf(sf_crmul(w,-sy*dz));
			pk[ix][ih] = sf_cmul(wx[ix][ih],cshift); );
#endif
		
		/* FFT */
		fft2(false,(kiss_fft_cpx**) pk);		
		
		si = slow[iz];
		sf_fslice_get(mms,iz,ms[0]);
		sf_fslice_get(mmr,iz,mr[0]);

		LOOPxh( wx[ix][ih] = sf_cmplx(0,0); );
		for (j=0; j<nr[iz]; j++) {
		    for (k=0; k<nr[iz]; k++) {
			if (skip[iz][j][k]) continue; 
                        /* skip S-R reference combinations */

			/* w-k phase shift */
#ifdef SF_HAS_COMPLEX_H
			cref =       csqrtf(w2*sm[iz][j])
			    +        csqrtf(w2*sm[iz][k]);
			LOOPxh2(cs = csqrtf(w2*sm[iz][j] + ks[ix][ih]);
				cr = csqrtf(w2*sm[iz][k] + kr[ix][ih]);
				cshift = cexpf((cref-cs-cr)*dz); 
				wk[ix][ih] = 
				pk[ix][ih]*cshift; ); 
#else
			cref =       sf_cadd(csqrtf(sf_crmul(w2,sm[iz][j])),
					     csqrtf(sf_crmul(w2,sm[iz][k])));
			LOOPxh2(cs = csqrtf(sf_cadd(sf_crmul(w2,sm[iz][j]),
						    sf_cmplx(ks[ix][ih],0.)));
				cr = csqrtf(sf_cadd(sf_crmul(w2,sm[iz][k]), 
						    sf_cmplx(kr[ix][ih],0.)));
				cshift = cexpf(
				    sf_crmul(
					sf_csub(cref,sf_cadd(cs,cr)),dz)); 
				wk[ix][ih] = sf_cmul(pk[ix][ih],cshift); ); 
#endif
		
			/* IFT */
			fft2(true,(kiss_fft_cpx**) wk);

			/* create MRS mask */
			LOOPxh( ma[ix][ih]= (ms[ix][ih]==j && 
					     mr[ix][ih]==k)?1.:0.; );

			/* accumulate wavefield */
#ifdef SF_HAS_COMPLEX_H
			LOOPxh( wx[ix][ih] += wk[ix][ih] * ma[ix][ih]; );
#else
			LOOPxh( wx[ix][ih] = sf_cadd(wx[ix][ih],
						     sf_crmul(wk[ix][ih],
							      ma[ix][ih])); );
#endif	
		    }
		}
		
		LOOPuh( qq[iu][ih]=img[iz*nu*nh+iu*nh+ih];);

		/* w-x at top */
#ifdef SF_HAS_COMPLEX_H
		LOOPxh( sy = 0.5*(si[ is[ix][ih] ] + 
				  si[ ir[ix][ih] ]);
			cshift = cexpf(-w*sy*dz);
			wx   [ix] [ih] = 
			qq[ii[ix]][ih] + 
			wx   [ix] [ih] * cshift; );
#else
		LOOPxh( sy = 0.5*(si[ is[ix][ih] ] + 
				  si[ ir[ix][ih] ]);
			cshift = cexpf(sf_crmul(w,-sy*dz));
			wx   [ix] [ih] = sf_cadd(
			    qq[ii[ix]][ih], 
			    sf_cmul(wx[ix] [ih],cshift)); );
#endif
	    } /* iz */

	    taper2(wx);
	    LOOPxh(dat[iw*nx*nh+ix*nh+ih]=wx[ix][ih];);

	} else { /* MIGRATION */
	    si = slow[0];

	    LOOPxh(wx[ix][ih]=dat[iw*nx*nh+ix*nh+ih];);
	    taper2(wx);

	    /* loop over migrated depths z */
	    for (iz=0; iz<nz-1; iz++) {
		/* imaging condition */
		LOOPuh(qq[iu][ih]=img[iz*nu*nh+iu*nh+ih];);
#ifdef SF_HAS_COMPLEX_H		
		LOOPxh(        qq[ii[ix]][ih] += 
			       wx[ix][ih] ; );
#else
		LOOPxh(        qq[ii[ix]][ih] = 
			       sf_cadd(qq[ii[ix]][ih],wx[ix][ih]) ; );
#endif
		LOOPuh(img[iz*nu*nh+iu*nh+ih]=qq[iu][ih];);

		/* w-x @ top */
		LOOPxh2( pk[ix][ih] = sf_cmplx(0.,0.););
#ifdef SF_HAS_COMPLEX_H
		LOOPxh( sy = 0.5*(si[ is[ix][ih] ] + 
				  si[ ir[ix][ih] ]);
			cshift = conjf(cexpf(-w*sy*dz));
			pk[ix][ih] = 
			wx[ix][ih] * cshift; );
#else
		LOOPxh( sy = 0.5*(si[ is[ix][ih] ] + 
				  si[ ir[ix][ih] ]);
			cshift = conjf(cexpf(sf_crmul(w,-sy*dz)));
			pk[ix][ih] = sf_cmul(wx[ix][ih],cshift); );
#endif

		/* FFT */
		fft2(false,(kiss_fft_cpx**) pk);

		si = slow[iz+1];
		sf_fslice_get(mms,iz,ms[0]);
		sf_fslice_get(mmr,iz,mr[0]);
		
		LOOPxh( wx[ix][ih] = sf_cmplx(0,0); );
		for (j=0; j<nr[iz]; j++) {
		    for (k=0; k<nr[iz]; k++) {
			if (skip[iz][j][k]) continue;
		
			/* w-k phase shift */
#ifdef SF_HAS_COMPLEX_H
			cref = csqrtf(w2*sm[iz][j])+csqrtf(w2*sm[iz][k]);
			LOOPxh2( cs = csqrtf(w2*sm[iz][j] + ks[ix][ih]);
				 cr = csqrtf(w2*sm[iz][k] + kr[ix][ih]);
				 cshift = conjf(cexpf((cref-cs-cr)*dz)); 
				 wk[ix][ih] = 
				 pk[ix][ih] * cshift; ); 
#else
			cref = sf_cadd(
			    csqrtf(sf_crmul(w2,sm[iz][j])),
			    csqrtf(sf_crmul(w2,sm[iz][k])));
			LOOPxh2( cs = csqrtf(
				     sf_cadd(sf_crmul(w2,sm[iz][j]),
					     sf_cmplx(ks[ix][ih],0.)));
				 cr = csqrtf(
				     sf_cadd(sf_crmul(w2,sm[iz][k]),
					     sf_cmplx(kr[ix][ih],0.)));
				 cshift = conjf(
				     cexpf(
					 sf_crmul(
					     sf_csub(cref,sf_cadd(cs,cr)),
					     dz))); 
				 wk[ix][ih] = sf_cmul(pk[ix][ih],cshift); ); 
#endif
			
			/* IFT */
			fft2(true,(kiss_fft_cpx**) wk);

			/* create MRS mask */
			LOOPxh( ma[ix][ih]= (ms[ix][ih]==j && 
					     mr[ix][ih]==k)?1.:0.; );

			/* accumulate wavefield */
#ifdef SF_HAS_COMPLEX_H
			LOOPxh( wx[ix][ih] += wk[ix][ih] * ma[ix][ih]; );
#else
			LOOPxh( wx[ix][ih] = 
				sf_cadd(wx[ix][ih],
					sf_crmul(wk[ix][ih],ma[ix][ih])); );
#endif
		    }
		} /* j loop */


		/* w-x @ bottom */
#ifdef SF_HAS_COMPLEX_H
		LOOPxh( sy = 0.5*(si[ is[ix][ih] ] + 
				  si[ ir[ix][ih] ]);
			cshift = conjf(cexpf(-w*sy*dz));
			wx[ix][ih] *= cshift; );
#else
		LOOPxh( sy = 0.5*(si[ is[ix][ih] ] + 
				  si[ ir[ix][ih] ]);
			cshift = conjf(cexpf(sf_crmul(w,-sy*dz)));
			wx[ix][ih] = sf_cmul(wx[ix][ih],cshift); );
#endif
		taper2(wx);
	    } /* iz */
	    
	    /* imaging condition @ bottom */
	    LOOPuh(qq[iu][ih]=img[(nz-1)*nu*nh+iu*nh+ih];);
#ifdef SF_HAS_COMPLEX_H  
	    LOOPxh(        qq[ii[ix]][ih] += 
			   wx[ix][ih]; );
#else
	    LOOPxh(        qq[ii[ix]][ih] = 
			   sf_cadd(qq[ii[ix]][ih],wx[ix][ih]); );
#endif
	    LOOPuh(img[(nz-1)*nu*nh+iu*nh+ih]=qq[iu][ih];);

	} /* else */
    } /* iw */
}
