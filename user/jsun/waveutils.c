/* Wave propagation utilities for modeling and imaging */
/*
  Copyright (C) 2014 University of Texas at Austin
  
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
#ifdef _OPENMP
#include <omp.h>
#endif

/*automatically generated headers*/
#include "waveutils.h"
#include "cfft2nsps.h"

/*******************************************************/
/* wave propagation utils*/
/*******************************************************/

#ifndef _waveutils_h

typedef struct Geopar {
  /* acquisition geometry and mesh setup */
  int nx, nz; /* domain of interest */
  int nxb, nzb; /* domain of computation (including boundary) */
  int   spx, spz, gpz, gpx, gpl;
  float dx, dz, ox, oz;
  /* snapshot interval */
  int   snpint;
  /*absorbing boundary*/
  int top, bot, lft, rht;
  /*source parameters*/
  int nt;
  float dt;
  float trunc; /* truncation */
  /*misc*/
  bool adj, verb, illum;
  int m2, m2b, pad1;
  /* pointers */
  sf_complex **ltf, **rtf;
  sf_complex **ltb, **rtb;
  sf_complex *ww;
  float *rr;
  /*extras*/
  bool roll;
  int rectz,rectx,repeat; /*refl smoothing parameters*/
  int sht0,shtbgn,shtend,shtnum,shtnum0,shtint;
  /*mpi*/
  int cpuid, numprocs;
  /*switch*/
  int mode;
} * geopar; 
/*^*/

typedef struct Mpipar {
    int cpuid;
    int numprocs;
} * mpipar; 
/*^*/

#endif

void reflgen(int nzb, int nxb, int spz, int spx,
             int rectz, int rectx, int nrep, /*smoothing parameters*/
	     float *refl/*reflectivity map*/)
/*< Generate reflectivity map with smoothing >*/
{   
    int iz, i, j, i0, irep;
    int nzx=nzb*nxb; 
    sf_triangle tr;
    int n[2],s[2],rect[2];
    bool diff[2],box[2];

    n[0]=nzb; n[1]=nxb;
    s[0]=1;   s[1]=nzb;
    rect[0]=rectz; rect[1]=rectx;
    diff[0]=false; diff[1]=false;
    box[0]=false; box[1]=false;
    
#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
    for (iz=0; iz < nzx; iz++) {
      refl[iz]=0;
    } 
    j=spx*nzb+spz; /*point source position*/
    refl[j]=1;
    
    /* 2-d triangle smoothing */
    for (i=0;i<2;i++) {
      if (rect[i] <= 1) continue;
      tr = sf_triangle_init (rect[i],n[i]);
      for (j=0; j < nzx/n[i]; j++) {
	i0 = sf_first_index (i,j,2,n,s);
	for (irep=0; irep < nrep; irep++) {
	  sf_smooth2 (tr,i0,s[i],diff[i],box[i],refl); // why adjoint?
	}
      }
      sf_triangle_close(tr);
    }
}

int lrosfor2(sf_complex ***wvfld, float **sill, sf_complex **rcd, geopar geop)
/*< low-rank one-step forward modeling >*/
{
    /* acquisition geometry and mesh setup */
    int nx, nz; /* domain of interest */
    int nxb, nzb; /* domain of computation (including boundary) */
    float dx, dz, ox, oz;
    int spx, spz, gpz, gpx, gpl;
    /* snapshot interval */
    int   snpint;
    /*absorbing boundary*/
    int top, bot, lft, rht;
    /*source parameters*/
    int nt;
    float dt;
    float trunc; /* truncation */
    /*misc*/
    bool adj, verb, illum;
    int m2, pad1;
    /* pointers */
    sf_complex **lt, **rt;
    sf_complex *ww;
    float *rr;

    int it,iz,im,ik,ix,i,j,wfit; /* index variables */
    int nk, nz2, nx2, nzx2;
    sf_complex c;
    sf_complex *cwave, *cwavem;
    sf_complex **wave, *curr;

    nx = geop->nx; nz = geop->nz;
    nxb = geop->nxb; nzb = geop->nzb;
    dx = geop->dx; dz = geop->dz; ox = geop->ox; oz = geop->oz; /*not acutally used*/
    snpint = geop->snpint;
    spx = geop->spx; spz = geop->spz; /*not acutally used*/
    gpz = geop->gpz; gpx = geop->gpx; gpl = geop->gpl;
    top = geop->top; bot = geop->bot; lft = geop->lft; rht = geop->rht;
    nt = geop->nt;
    dt = geop->dt;
    trunc = geop->trunc;
    adj = geop->adj; verb = geop->verb; illum = geop->illum;
    pad1 = geop->pad1;
    if (geop->mode==0) {
      m2 = geop->m2;
      lt = geop->ltf; rt = geop->rtf;
    } else {
      m2 = geop->m2b;
      lt = geop->ltb; rt = geop->rtb;
    }
    ww = geop->ww; rr = geop->rr;

    /*Matrix dimensions*/
    nk = cfft2_init(pad1,nzb,nxb,&nz2,&nx2);
    nzx2 = nz2*nx2;

    curr   = sf_complexalloc(nzx2);
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave   = sf_complexalloc2(nzx2,m2);

    /*icfft2_allocate(cwavem);*/

#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = sf_cmplx(0.,0.);
    }

    if (illum) {
#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
      for (ix=0; ix < nx; ix++) {
	for (iz=0; iz < nz; iz++) {
	  sill[ix][iz] = 0.f;
	}
      }
    }

    /*Main loop*/
    wfit = 0;
    for (it = 0; it < nt; it++) {
	if (verb) sf_warning("Forward source it=%d/%d;", it, nt-1);
	
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
		if ((it*dt)<=trunc) {
#ifdef SF_HAS_COMPLEX_H
		  c = ww[it] * rr[i]; // source term
#else
		  c = sf_crmul(ww[it], rr[i]); // source term
#endif
		} else {
		  c = sf_cmplx(0.,0.);
		}
		for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
		    c += lt[im][i]*wave[im][j];
#else
		    c = sf_cadd(c,sf_cmul(lt[im][i], wave[im][j]));
#endif
		}
		curr[j] = c;
	    }
	}

#ifdef _OPENMP
#pragma omp parallel for private(ix,j)
#endif	 
	for ( ix =0 ; ix < gpl; ix++) {
	    j = (gpz+top)+(ix+gpx+lft)*nz2; /* padded grid */
	    rcd[ix][it] = curr[j];
	}
	
	if ( it%snpint == 0 ) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,j)
#endif
	    for ( ix = 0; ix < nx; ix++) {
		for ( iz = 0; iz<nz; iz++ ) { 
		    j = (iz+top)+(ix+lft)*nz2; /* padded grid */
		    wvfld[wfit][ix][iz] = curr[j];
		    if (illum) sill[ix][iz] += pow(hypotf(crealf(curr[j]),cimagf(curr[j])),2);
		}
	    }
	    wfit++;
	}
    } /*Main loop*/
    if (verb) sf_warning(".");
    cfft2_finalize();
    return wfit;
}

int lrosback2(sf_complex **img, sf_complex ***wvfld, float **sill, sf_complex **rcd, geopar geop)
/*< low-rank one-step backward propagation + imaging >*/
{
    /* acquisition geometry and mesh setup */
    int nx, nz; /* domain of interest */
    int nxb, nzb; /* domain of computation (including boundary) */
    float dx, dz, ox, oz;
    int spx, spz, gpz, gpx, gpl;
    /* snapshot interval */
    int snpint;
    /*absorbing boundary*/
    int top, bot, lft, rht;
    /*source parameters*/
    int nt;
    float dt;
    float trunc; /* truncation */
    /*misc*/
    bool adj, verb, illum;
    int m2, pad1;
    /* pointers */
    sf_complex **lt, **rt;

    int it,iz,im,ik,ix,i,j,wfit;     /* index variables */
    int nk, nz2, nx2, nzx2;
    sf_complex c;
    sf_complex *cwave, *cwavem, *currm;
    sf_complex **wave, *curr;
    sf_complex **ccr;

    nx = geop->nx; nz = geop->nz;
    nxb = geop->nxb; nzb = geop->nzb;
    dx = geop->dx; dz = geop->dz; ox = geop->ox; oz = geop->oz;
    snpint = geop->snpint;
    spx = geop->spx; spz = geop->spz;
    gpz = geop->gpz; gpx = geop->gpx; gpl = geop->gpl;
    top = geop->top; bot = geop->bot; lft = geop->lft; rht = geop->rht;
    nt = geop->nt;
    dt = geop->dt;
    trunc = geop->trunc;
    adj = geop->adj; verb = geop->verb; illum = geop->illum;
    pad1 = geop->pad1;
    if (geop->mode==0) {
      m2 = geop->m2;
      lt = geop->ltf; rt = geop->rtf;
    } else {
      m2 = geop->m2b;
      lt = geop->ltb; rt = geop->rtb;
    }

    /*Matrix dimensions*/
    nk = cfft2_init(pad1,nzb,nxb,&nz2,&nx2);
    nzx2 = nz2*nx2;

    curr = sf_complexalloc(nzx2);
    cwave  = sf_complexalloc(nk);
    wave = sf_complexalloc2(nzx2,m2);
    ccr = sf_complexalloc2(nz, nx);

    if (!adj) {
	currm  = sf_complexalloc(nzx2);
	/*icfft2_allocate(cwave);*/
    } else {
	cwavem = sf_complexalloc(nk);
	/*icfft2_allocate(cwavem);*/
    }

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
	  ccr[ix][iz] = sf_cmplx(0.,0.);
	}
    }

    if (adj) { /* migration */
      /* step backward in time (PSPI) */
      /*Main loop*/
      wfit = (int)(nt-1)/snpint;
      for (it = nt-1; it>=0; it--) {
	if  (verb) sf_warning("Backward receiver it=%d/%d;", it, nt-1);
#ifdef _OPENMP
#pragma omp parallel for private(ix,j)
#endif
        for (ix=0; ix<gpl; ix++)  {
	  j = (gpz+top)+(ix+gpx+lft)*nz2; /* padded grid */
#ifdef SF_HAS_COMPLEX_H
	  curr[j]+=rcd[ix][it]; /* data injection */
#else
	  curr[j]=sf_cadd(curr[j],rcd[ix][it]);
#endif
	}

	/*cross-correlation imaging condition*/
	if (it%snpint == 0 ) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,j)
#endif
	  for (ix=0; ix<nx; ix++) {
	    for (iz=0; iz<nz; iz++) {
	      j = (iz+top)+(ix+lft)*nz2; /* padded grid */
#ifdef SF_HAS_COMPLEX_H
	      ccr[ix][iz] += conjf(wvfld[wfit][ix][iz])*curr[j];
#else
	      ccr[ix][iz] = sf_cadd(ccr[ix][iz],sf_cmul(conjf(wvfld[wfit][ix][iz]),curr[j]));
#endif
	    }
	  }
	  wfit--;
	}

	/*matrix multiplication*/
	/*PSPI*/
	cfft2(curr,cwave);
	for (im = 0; im < m2; im++) {
#ifdef _OPENMP
#pragma omp parallel for private(ik)
#endif
	  for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	    cwavem[ik] = cwave[ik]*conjf(rt[ik][im]);
#else
	    cwavem[ik] = sf_cmul(cwave[ik],conjf(rt[ik][im]));
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
	      c += conjf(lt[im][i])*wave[im][j];
#else
	      c = sf_cadd(c,sf_cmul(conjf(lt[im][i]), wave[im][j]));
#endif
	    }
	    curr[j] = c;
	  }
	}

      } /*Main loop*/
      if (verb) sf_warning(".");
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif    
      for (ix=0; ix<nx; ix++) {
	for (iz=0; iz<nz; iz++) {
	  if (illum) {
#ifdef SF_HAS_COMPLEX_H
	    img[ix][iz] = ccr[ix][iz]/(sill[ix][iz]+SF_EPS);
#else
	    img[ix][iz] = sf_crmul(ccr[ix][iz],1./(sill[ix][iz]+SF_EPS));
#endif
	  } else img[ix][iz] = ccr[ix][iz];
	}
      } 
    } else { /* modeling */
      /* adjoint of source illumination */
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif    
      for (ix=0; ix<nx; ix++) {
	for (iz=0; iz<nz; iz++) {
	  if (illum) {
#ifdef SF_HAS_COMPLEX_H
	    ccr[ix][iz] = img[ix][iz]/(sill[ix][iz]+SF_EPS);
#else
	    ccr[ix][iz] = sf_crmul(img[ix][iz],1./(sill[ix][iz]+SF_EPS));
#endif
	  } else ccr[ix][iz] = img[ix][iz];
	}
      } 
      /* step forward in time (NSPS) */
      /*Main loop*/
      wfit=0;
      for (it=0; it<nt; it++) {
	if (verb) sf_warning("Forward receiver it=%d/%d;", it, nt-1);

	/*matrix multiplication*/
	/*NSPS*/
	for (im = 0; im < m2; im++) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,i,j) shared(currm,lt,curr)
#endif
	  for (ix = 0; ix < nxb; ix++) {
	    for (iz=0; iz < nzb; iz++) {
	      i = iz+ix*nzb;  /* original grid */
	      j = iz+ix*nz2; /* padded grid */
#ifdef SF_HAS_COMPLEX_H
	      currm[j] = lt[im][i]*curr[j];
#else
	      currm[j] = sf_cmul(lt[im][i], curr[j]);
#endif
	    }
	  }
	  cfft2(currm,wave[im]);
	}
#ifdef _OPENMP
#pragma omp parallel for private(ik,im,c) shared(wave,rt,cwave)
#endif	    
	for (ik = 0; ik < nk; ik++) {
	  c = sf_cmplx(0.,0.);
	  for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
	    c += wave[im][ik]*rt[ik][im];
#else
	    c = sf_cadd(c,sf_cmul(wave[im][ik],rt[ik][im])); 
#endif
	  }
	  cwave[ik] = c;
	}
	icfft2(curr,cwave);

	/*adjoint of cross-correlation imaging condition*/
	if (it%snpint == 0 ) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,j)
#endif
	  for (ix=0; ix<nx; ix++) {
	    for (iz=0; iz<nz; iz++) {
	      j = (iz+top)+(ix+lft)*nz2; /* padded grid */
#ifdef SF_HAS_COMPLEX_H
	      curr[j] += (wvfld[wfit][ix][iz])*ccr[ix][iz];/* adjoint of ccr[ix][iz] += conjf(wvfld[wfit][ix][iz])*curr[j]; ??? */
#else
	      curr[j] = sf_cadd(curr[j],sf_cmul((wvfld[wfit][ix][iz]),ccr[ix][iz]));
#endif
	    }
	  }
	  wfit++;
	}

#ifdef _OPENMP
#pragma omp parallel for private(ix,j)
#endif
        for (ix=0; ix<gpl; ix++)  {
	  j = (gpz+top)+(ix+gpx+lft)*nz2; /* padded grid */
	  rcd[ix][it]=curr[j];
	}
      } /*Main loop*/
    }
    cfft2_finalize();
    return 0;
}

int rcvill2(float **sill, sf_complex **rcd, geopar geop)
/*< low-rank one-step backward propagation for illumination calculation >*/
{
    /* acquisition geometry and mesh setup */
    int nx, nz; /* domain of interest */
    int nxb, nzb; /* domain of computation (including boundary) */
    float dx, dz, ox, oz;
    int spx, spz, gpz, gpx, gpl;
    /* snapshot interval */
    int snpint;
    /*absorbing boundary*/
    int top, bot, lft, rht;
    /*source parameters*/
    int nt;
    float dt;
    float trunc; /* truncation */
    /*misc*/
    bool adj, verb, illum;
    int m2, pad1;
    /* pointers */
    sf_complex **lt, **rt;

    int it,iz,im,ik,ix,i,j,wfit;     /* index variables */
    int nk, nz2, nx2, nzx2;
    sf_complex c;
    sf_complex *cwave, *cwavem;
    sf_complex **wave, *curr;

    nx = geop->nx; nz = geop->nz;
    nxb = geop->nxb; nzb = geop->nzb;
    dx = geop->dx; dz = geop->dz; ox = geop->ox; oz = geop->oz;
    snpint = geop->snpint;
    spx = geop->spx; spz = geop->spz;
    gpz = geop->gpz; gpx = geop->gpx; gpl = geop->gpl;
    top = geop->top; bot = geop->bot; lft = geop->lft; rht = geop->rht;
    nt = geop->nt;
    dt = geop->dt;
    trunc = geop->trunc;
    adj = geop->adj; verb = geop->verb; illum = geop->illum;
    pad1 = geop->pad1;
    if (geop->mode==0) {
      m2 = geop->m2;
      lt = geop->ltf; rt = geop->rtf;
    } else {
      m2 = geop->m2b;
      lt = geop->ltb; rt = geop->rtb;
    }

    /*Matrix dimensions*/
    nk = cfft2_init(pad1,nzb,nxb,&nz2,&nx2);
    nzx2 = nz2*nx2;

    curr = sf_complexalloc(nzx2);
    cwave  = sf_complexalloc(nk);
    wave = sf_complexalloc2(nzx2,m2);

    cwavem = sf_complexalloc(nk);
    /*icfft2_allocate(cwavem);*/

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
	  sill[ix][iz] = 0.f;
	}
    }

    /* step backward in time (PSPI) */
    /*Main loop*/
    wfit = (int)(nt-1)/snpint;
    for (it = nt-1; it>=0; it--) {
      if  (verb) sf_warning("Backward receiver it=%d/%d;", it, nt-1);
#ifdef _OPENMP
#pragma omp parallel for private(ix,j)
#endif
      for (ix=0; ix<gpl; ix++)  {
	j = (gpz+top)+(ix+gpx+lft)*nz2; /* padded grid */
#ifdef SF_HAS_COMPLEX_H
	curr[j]+=rcd[ix][it]; /* data injection */
#else
	curr[j]=sf_cadd(curr[j],rcd[ix][it]);
#endif
      }

      /*receiver illumination function*/
      if (it%snpint == 0 ) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,j)
#endif
	for (ix=0; ix<nx; ix++) {
	  for (iz=0; iz<nz; iz++) {
	    j = (iz+top)+(ix+lft)*nz2; /* padded grid */
	    sill[ix][iz] += pow(hypotf(crealf(curr[j]),cimagf(curr[j])),2);
	  }
	}
	wfit--;
      }

      /*matrix multiplication*/
      /*PSPI*/
      cfft2(curr,cwave);
      for (im = 0; im < m2; im++) {
#ifdef _OPENMP
#pragma omp parallel for private(ik)
#endif
	for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	  cwavem[ik] = cwave[ik]*conjf(rt[ik][im]);
#else
	  cwavem[ik] = sf_cmul(cwave[ik],conjf(rt[ik][im]));
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
	    c += conjf(lt[im][i])*wave[im][j];
#else
	    c = sf_cadd(c,sf_cmul(conjf(lt[im][i]), wave[im][j]));
#endif
	  }
	  curr[j] = c;
	}
      }

    } /*Main loop*/
    if (verb) sf_warning(".");

    cfft2_finalize();
    return 0;
}

int lrosfor2q(sf_complex ***wvfld, geopar geop)
/*< low-rank one-step forward modeling >*/
{
    /* acquisition geometry and mesh setup */
    int nx, nz; /* domain of interest */
    int nxb, nzb; /* domain of computation (including boundary) */
    float dx, dz, ox, oz;
    int spx, spz, gpz, gpx, gpl;
    /* snapshot interval */
    int   snpint;
    /*absorbing boundary*/
    int top, bot, lft, rht;
    /*source parameters*/
    int nt;
    float dt;
    float trunc; /* truncation */
    /*misc*/
    bool adj, verb, illum;
    int m2, pad1;
    /* pointers */
    sf_complex **lt, **rt;
    sf_complex *ww;
    float *rr;

    int it,iz,im,ik,ix,i,j,wfit; /* index variables */
    int nk, nz2, nx2, nzx2;
    sf_complex c;
    sf_complex *cwave, *cwavem;
    sf_complex **wave, *curr;

    nx = geop->nx; nz = geop->nz;
    nxb = geop->nxb; nzb = geop->nzb;
    dx = geop->dx; dz = geop->dz; ox = geop->ox; oz = geop->oz; /*not acutally used*/
    snpint = geop->snpint;
    spx = geop->spx; spz = geop->spz; /*not acutally used*/
    gpz = geop->gpz; gpx = geop->gpx; gpl = geop->gpl;
    top = geop->top; bot = geop->bot; lft = geop->lft; rht = geop->rht;
    nt = geop->nt;
    dt = geop->dt;
    trunc = geop->trunc;
    adj = geop->adj; verb = geop->verb; illum = geop->illum;
    pad1 = geop->pad1;
    if (geop->mode==0) {
      m2 = geop->m2;
      lt = geop->ltf; rt = geop->rtf;
    } else {
      m2 = geop->m2b;
      lt = geop->ltb; rt = geop->rtb;
    }
    ww = geop->ww; rr = geop->rr;

    /*Matrix dimensions*/
    nk = cfft2_init(pad1,nzb,nxb,&nz2,&nx2);
    nzx2 = nz2*nx2;

    curr   = sf_complexalloc(nzx2);
    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave   = sf_complexalloc2(nzx2,m2);

    /*icfft2_allocate(cwavem);*/

#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = sf_cmplx(0.,0.);
    }

    /*Main loop*/
    wfit = 0;
    for (it = 0; it < nt; it++) {
	if (verb) sf_warning("Forward source it=%d/%d;", it, nt-1);
	
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
		if ((it*dt)<=trunc) {
#ifdef SF_HAS_COMPLEX_H
		  c = ww[it] * rr[i]; // source term
#else
		  c = sf_crmul(ww[it], rr[i]); // source term
#endif
		} else {
		  c = sf_cmplx(0.,0.);
		}
		for (im = 0; im < m2; im++) {
#ifdef SF_HAS_COMPLEX_H
		    c += lt[im][i]*wave[im][j];
#else
		    c = sf_cadd(c,sf_cmul(lt[im][i], wave[im][j]));
#endif
		}
		curr[j] = c;
	    }
	}
	
	if ( it%snpint == 0 ) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,j)
#endif
	    for ( ix = 0; ix < nx; ix++) {
		for ( iz = 0; iz<nz; iz++ ) { 
		    j = (iz+top)+(ix+lft)*nz2; /* padded grid */
		    wvfld[wfit][ix][iz] = curr[j];
		}
	    }
	    wfit++;
	}
    } /*Main loop*/
    if (verb) sf_warning(".");
    cfft2_finalize();
    return wfit;
}

int lrosback2q(sf_complex ***wvfld, sf_complex **rcd, geopar geop)
/*< low-rank one-step backward propagation + imaging >*/
{
    /* acquisition geometry and mesh setup */
    int nx, nz; /* domain of interest */
    int nxb, nzb; /* domain of computation (including boundary) */
    float dx, dz, ox, oz;
    int spx, spz, gpz, gpx, gpl;
    /* snapshot interval */
    int snpint;
    /*absorbing boundary*/
    int top, bot, lft, rht;
    /*source parameters*/
    int nt;
    float dt;
    float trunc; /* truncation */
    /*misc*/
    bool adj, verb, illum;
    int m2, pad1;
    /* pointers */
    sf_complex **lt, **rt;

    int it,iz,im,ik,ix,i,j,wfit;     /* index variables */
    int nk, nz2, nx2, nzx2;
    sf_complex c;
    sf_complex *cwave, *cwavem, *currm;
    sf_complex **wave, *curr;
    sf_complex **ccr;

    nx = geop->nx; nz = geop->nz;
    nxb = geop->nxb; nzb = geop->nzb;
    dx = geop->dx; dz = geop->dz; ox = geop->ox; oz = geop->oz;
    snpint = geop->snpint;
    spx = geop->spx; spz = geop->spz;
    gpz = geop->gpz; gpx = geop->gpx; gpl = geop->gpl;
    top = geop->top; bot = geop->bot; lft = geop->lft; rht = geop->rht;
    nt = geop->nt;
    dt = geop->dt;
    trunc = geop->trunc;
    adj = geop->adj; verb = geop->verb; illum = geop->illum;
    pad1 = geop->pad1;
    if (geop->mode==0) {
      m2 = geop->m2;
      lt = geop->ltf; rt = geop->rtf;
    } else {
      m2 = geop->m2b;
      lt = geop->ltb; rt = geop->rtb;
    }

    /*Matrix dimensions*/
    nk = cfft2_init(pad1,nzb,nxb,&nz2,&nx2);
    nzx2 = nz2*nx2;

    curr = sf_complexalloc(nzx2);
    cwave  = sf_complexalloc(nk);
    wave = sf_complexalloc2(nzx2,m2);
    ccr = sf_complexalloc2(nz, nx);

    cwavem = sf_complexalloc(nk);
    /*icfft2_allocate(cwavem);*/

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
	  ccr[ix][iz] = sf_cmplx(0.,0.);
	}
    }

    /* step backward in time (PSPI) */
    /*Main loop*/
    wfit = (int)(nt-1)/snpint;
    for (it = nt-1; it>=0; it--) {
      if  (verb) sf_warning("Backward receiver it=%d/%d;", it, nt-1);
#ifdef _OPENMP
#pragma omp parallel for private(ix,j)
#endif
      for (ix=0; ix<gpl; ix++)  {
	j = (gpz+top)+(ix+gpx+lft)*nz2; /* padded grid */
#ifdef SF_HAS_COMPLEX_H
	curr[j]+=rcd[ix][it]; /* data injection */
#else
	curr[j]=sf_cadd(curr[j],rcd[ix][it]);
#endif
      }

      /*record the backward wavefield*/
      if (it%snpint == 0 ) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,j)
#endif
	for (ix=0; ix<nx; ix++) {
	  for (iz=0; iz<nz; iz++) {
	    j = (iz+top)+(ix+lft)*nz2; /* padded grid */
	    wvfld[wfit][ix][iz] = curr[j];
	  }
	}
	wfit--;
      }

      /*matrix multiplication*/
      /*PSPI*/
      cfft2(curr,cwave);
      for (im = 0; im < m2; im++) {
#ifdef _OPENMP
#pragma omp parallel for private(ik)
#endif
	for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	  cwavem[ik] = cwave[ik]*conjf(rt[ik][im]);
#else
	  cwavem[ik] = sf_cmul(cwave[ik],conjf(rt[ik][im]));
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
	    c += conjf(lt[im][i])*wave[im][j];
#else
	    c = sf_cadd(c,sf_cmul(conjf(lt[im][i]), wave[im][j]));
#endif
	  }
	  curr[j] = c;
	}
      }

    } /*Main loop*/
    if (verb) sf_warning(".");
    
    cfft2_finalize();
    return 0;
}

int ccrimg(sf_complex **img, sf_complex ***wvfld, sf_complex ***wvfld_b, float **sill, geopar geop) 
/*< cross-correlation imaging condition >*/
{
  int nx, nz, nt, snpint, wfnt;
  int ix, iz, it;
  bool illum;

  nx = geop->nx; nz = geop->nz; nt = geop->nt;
  snpint = geop->snpint;
  illum = geop->illum;

  wfnt = (int)(nt-1)/snpint+1;

  /*initialization*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
  for (ix=0; ix<nx; ix++)
    for (iz=0; iz<nz; iz++)
      img[ix][iz] = sf_cmplx(0.,0.);

  /*cross-correlation imaging condition*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,it)
#endif
  for (it=0; it<wfnt; it++) {
    for (ix=0; ix<nx; ix++) {
      for (iz=0; iz<nz; iz++) {
#ifdef SF_HAS_COMPLEX_H
	img[ix][iz] += conjf(wvfld[it][ix][iz])*wvfld_b[it][ix][iz];
#else
	img[ix][iz] = sf_cadd(img[ix][iz],sf_cmul(conjf(wvfld[it][ix][iz]),wvfld_b[it][ix][iz]));
#endif
      }
    }
  }

  /*illumination*/
  if (illum) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
    for (ix=0; ix<nx; ix++) {
      for (iz=0; iz<nz; iz++) {
#ifdef SF_HAS_COMPLEX_H
	img[ix][iz] = img[ix][iz]/(sill[ix][iz]+SF_EPS);
#else
	img[ix][iz] = sf_crmul(img[ix][iz],1./(sill[ix][iz]+SF_EPS));
#endif
      }
    }
  }

  return 0;
}
