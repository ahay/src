/* Forward and backward lowrank FD modeling on a staggered grid */
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

#include <rsf.h>
#include "cfft2w.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef _lrosutil_h

#define SRCRANGE 10
#define SRCALPHA 0.5
#define SRCTRUNC 100
#define SRCDECAY false
/*^*/

typedef struct Geopar {
    int   nx;
    int   nz;
    int   nxb;
    int   nzb;
    float dx;
    float dz;
    float ox;
    float oz;
    int   spx;
    int   spz;
    int   gp;
    int   snpint;
/*absorbing boundary*/
    int top;
    int bot;
    int lft;
    int rht;
} * geopar; /*geometry parameters*/
/*^*/

typedef struct Srcpar {
    sf_complex *wavelet;
    int nt;
    float dt;
    int range;
    float alpha;
    bool decay;
    float trunc;
} * srcpar; /*source parameter*/
/*^*/

#endif

srcpar createsrc(void)
/*<Create source struct>*/
{
    srcpar srcp;

    srcp = (srcpar)sf_alloc(1, sizeof(*srcp));
    srcp->wavelet = NULL;
    srcp->nt = 0;
    srcp->dt = 0.0;
    srcp->range = SRCRANGE;
    srcp->alpha = SRCALPHA;
    srcp->decay = SRCDECAY;
    srcp->trunc = SRCTRUNC;
    return srcp;
}

void loadsrc(srcpar srcp, sf_file srcfile)
/*<allocate source wavelet>*/
{
    if (srcp->nt == 0) sf_error("Need nt in srcpar!");
    if (srcp->wavelet != NULL) sf_error("source has been loaded!");
    srcp->wavelet = sf_complexalloc(srcp->nt);
    sf_complexread(srcp->wavelet, srcp->nt, srcfile);
}

void freesrc(srcpar srcp)
/*< Free allocated storage >*/
{
    free(srcp->wavelet);
    free(srcp);
}


void explsourcet(sf_complex *curr/*@out@*/,
		 sf_complex *vwavlet, 
		 int vit, float vdt,
                 int vsx, int vsdepth, 
		 int nx2, int nz2,
		 srcpar vps/*decay parameters*/)
/*< explosive source >*/ 
{
    float phi = 0.0;
    int cent = vps->range/2;
    int ix, iz;

    if (vps->decay ==1){
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,phi)
#endif
	for (ix=0; ix<2*cent; ix++)
	    for (iz=0; iz<2*cent; iz++) {
		phi = exp( -1*vps->alpha*vps->alpha*((ix-cent)*(ix-cent)+(iz-cent)*(iz-cent)) );
#ifdef SF_HAS_COMPLEX_H
		curr[(vsx-cent+ix)*nz2+(vsdepth-cent+iz)] += vwavlet[vit]*phi;
#else
		curr[(vsx-cent+ix)*nz2+(vsdepth-cent+iz)] += sf_crmul(vwavlet[vit],phi);
#endif
	    }
    } else {
	curr[vsx*nz2+vsdepth] += vwavlet[vit];
    } 
}

geopar creategeo(void)
/*< Create geometry used in RTM >*/
{
    geopar geop;
    geop = (geopar) sf_alloc(1, sizeof(*geop));
    return geop;
}

 
int lrosfor2(float ***wavfld, sf_complex **rcd, bool verb,
	     sf_complex **lt, sf_complex **rt, int m2,
            geopar geop, srcpar srcp, int pad1)
/*< low-rank one-step forward modeling >*/
{
    int it,iz,im,ik,ix,i,j;     /* index variables */
    int nxb,nzb,dx,dz,spx,spz,gp,snpint,dt,nth,wfit;
    int nt,nz,nx, nk, nzx, nz2, nx2, nzx2;
    sf_complex c;
    sf_complex *cwave, *cwavem;
    sf_complex **wave, *curr;

    nx = geop->nx;
    nz = geop->nz;
    nxb = geop->nxb;
    nzb = geop->nzb;
    dx = geop->dx;
    dz = geop->dz;

    spx = geop->spx;
    spz = geop->spz;
    gp  = geop->gp;
    snpint = geop->snpint;
    
    nt = srcp->nt;
    dt = srcp->dt;

#ifdef _OPENMP
#pragma omp parallel  
{
    nth = omp_get_num_threads();
}
#endif
    sf_warning(">>>> Using %d threads <<<<<", nth);
    
    /*Matrix dimensions*/
    nk = cfft2_init(pad1,nzb,nxb,&nz2,&nx2);
    nzx = nzb*nxb;
    nzx2 = nz2*nx2;

    curr   = sf_complexalloc(nzx2);
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave   = sf_complexalloc2(nzx2,m2);

    icfft2_allocate(cwavem);

#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = sf_cmplx(0.,0.);
    }

    /*Main loop*/
    wfit = 0;
    for (it = 0; it < nt; it++) {
	if (verb) sf_warning("Forward it=%d/%d;", it, nt-1);
	
	/*matrix multiplication*/
	cfft2(curr,cwave);

	for (im = 0; im < m2; im++) {
#ifdef _OPENMP
#pragma omp parallel for private(ik)
#endif
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		cwavem[ik] = cwave[ik]*rt[ik][im];
#else
		cwavem[ik] = sf_cmul(cwave[ik],rt[ik][im]);
#endif
	    }
	    icfft2(wave[im],cwavem);
	}

#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,i,j,im,c) shared(curr,lt,wave)
#endif
	for (ix = 0; ix < nxb; ix++) {
	    for (iz=0; iz < nzb; iz++) {
		i = iz+ix*nzb;  /* original grid */
		j = iz+ix*nz2; /* padded grid */
		c = sf_cmplx(0.,0.); // source term originally, now adds later
		for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
		    c += lt[im][i]*wave[im][j];
#else
		    c += sf_cmul(lt[im][i], wave[im][j]);
#endif
		}
		curr[j] = c;
	    }
	    if ((it*dt)<=srcp->trunc) {
	      explsourcet(curr, srcp->wavelet, it, dt, spx+geop->lft, spz+geop->top, nx2, nz2, srcp); // dt redundant???
	    }
	}

#ifdef _OPENMP
#pragma omp parallel for private(ix,j)
#endif	 
	for ( ix =0 ; ix < nx; ix++) {
	    j = (gp+geop->top)+(ix+geop->lft)*nz2; /* padded grid */
	    rcd[ix][it] = curr[j];
	}
	
	if ( it%snpint == 0 ) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,j)
#endif
	    for ( ix = 0; ix < nx; ix++) {
		for ( iz = 0; iz<nz; iz++ ) { 
		    j = (iz+geop->top)+(ix+geop->lft)*nz2; /* padded grid */
		    wavfld[wfit][ix][iz] = crealf(curr[j]);
		}
	    }
	    wfit++;
	}
    } /*Main loop*/
    if (verb) sf_warning(".");
    cfft2_finalize();
    return wfit;
    
}


int lrosback2(float **img1, float **img2, float ***wavfld, sf_complex **rcd, 
	      bool verb, bool wantwf, sf_complex **lt, sf_complex **rt, int m2,
               geopar geop, srcpar srcp, int pad1, float ***wavfld2)  
/*< low-rank one-step backward propagation + imaging >*/
{
    int it,iz,im,ik,ix,i,j;     /* index variables */
    int nxb,nzb,dx,dz,gp,snpint,dt,wfit;
    int nt,nz,nx, nk, nzx, nz2, nx2, nzx2;
    sf_complex c;
    sf_complex *cwave, *cwavem;
    sf_complex **wave, *curr;
    float **sill, **ccr;

    nx = geop->nx;
    nz = geop->nz;
    nxb = geop->nxb;
    nzb = geop->nzb;
    dx = geop->dx;
    dz = geop->dz;
    
    gp  = geop->gp;
    snpint = geop->snpint;
    
    nt = srcp->nt;
    dt = srcp->dt;

    sill = sf_floatalloc2(nz, nx);
    ccr  = img1;

    nk = cfft2_init(pad1,nzb,nxb,&nz2,&nx2);
    nzx = nzb*nxb;
    nzx2 = nz2*nx2;

    curr = sf_complexalloc(nzx2);
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave = sf_complexalloc2(nzx2,m2);

    icfft2_allocate(cwavem);

#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = sf_cmplx(0.,0.);
    }

#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix = 0; ix < nx; ix++) {
	for (iz = 0; iz < nz; iz++) {
	    sill[ix][iz] = 0.0;
	}
    }
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix = 0; ix < nx; ix++) {
	for (iz = 0; iz < nz; iz++) {
	    ccr[ix][iz] = 0.0;
	}
    }
        
    /*Main loop*/
    wfit = (int)(nt-1)/snpint; // dblcheck
  
    for (it = nt-1; it>=0; it--) {
	if  (verb) sf_warning("Backward it=%d/%d;", it, nt-1);

#ifdef _OPENMP
#pragma omp parallel for private(ix,j)
#endif
        for (ix=0; ix<nx; ix++)  {
	    j = (gp+geop->top)+(ix+geop->lft)*nz2; /* padded grid */
	    curr[j]+=rcd[ix][it];
	}

	/*matrix multiplication*/
	cfft2(curr,cwave);

	for (im = 0; im < m2; im++) {
#ifdef _OPENMP
#pragma omp parallel for private(ik)
#endif
	    for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
		cwavem[ik] = cwave[ik]*rt[ik][im];
#else
		cwavem[ik] = sf_cmul(cwave[ik],rt[ik][im]);
#endif
	    }
	    icfft2(wave[im],cwavem);
	}

#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,i,j,im,c) shared(curr,lt,wave)
#endif
	for (ix = 0; ix < nxb; ix++) {
	    for (iz=0; iz < nzb; iz++) {
		i = iz+ix*nzb;  /* original grid */
		j = iz+ix*nz2; /* padded grid */
		c = sf_cmplx(0.,0.); // initialize
		for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
		    c += lt[im][i]*wave[im][j];
#else
		    c += sf_cmul(lt[im][i], wave[im][j]);
#endif
		}
		curr[j] = c;
	    }
	}

        if ( wantwf && it%snpint == 0 ) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,j)
#endif
	    for ( ix = 0; ix < nx; ix++) {
		for ( iz = 0; iz<nz; iz++ ) { 
		    j = (iz+geop->top)+(ix+geop->lft)*nz2; /* padded grid */
		    wavfld2[wfit][ix][iz] = crealf(curr[j]);
		}
	    }
        }
	/*cross-correlation imaging condition*/
	if (it%snpint == 0 ) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,j)
#endif
	    for (ix=0; ix<nx; ix++) {
		for (iz=0; iz<nz; iz++) {
		    j = (iz+geop->top)+(ix+geop->lft)*nz2; /* padded grid */
		    ccr[ix][iz] += wavfld[wfit][ix][iz]*crealf(curr[j]);
		}
	    }
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
	    for (ix=0; ix<nx; ix++) {
		for (iz=0; iz<nz; iz++) {
		    sill[ix][iz] += wavfld[wfit][ix][iz]*wavfld[wfit][ix][iz];
		}
	    }
	    wfit--;
	}
    } /*Main loop*/
    if (verb) sf_warning(".");
    cfft2_finalize();
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif    
    for (ix=0; ix<nx; ix++) {
	for (iz=0; iz<nz; iz++) {
	    img2[ix][iz] = ccr[ix][iz]/(sill[ix][iz]+SF_EPS);//
	}
    } 
    
    return 0;
    
}
