/* 2-D FFT-based zero-offset exploding reflector modeling/migration linear operator */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include "cfft2nsps.h"
#include "timer.h"

/* wave propagation matrices */
static sf_complex **ltf, **rtf, **ltb, **rtb;

/*******************************************************/
/*wave propagation utils*/

typedef struct Geopar {
    /*geometry parameters*/
    int   nx;
    int   nz;
    float dx;
    float dz;
    float ox;
    int   gpz;
    /*time parameters*/
    int nt;
    float dt;
    int snap;
    int wfnt;
    /*miscellaneous*/
    int pad1;
    bool verb;
    /*propagation matrix*/
    int nzx2;
    int nk;
    int m2;
    /*fft parameters*/
    int nkz;
    int nkx;
    float dkz;
    float dkx;
    float kz0;
    float kx0;
    /*lowrank parameters*/
    int m; 
    int n2; // left: m*n2
    int n;  //right: n2*n
    /* wave propagation matrices*/
    /*sf_complex **ltf, **rtf, **ltb, **rtb;*/
} * geopar;
/*^*/

void lrexp_init (sf_complex **lt1,sf_complex **rt1,sf_complex **lt2,sf_complex **rt2)
/*< initialize the wave propagation matricies >*/
{
  ltf = lt1;
  rtf = rt1;
  ltb = lt2;
  rtb = rt2;
}

void lrexp_close(void)
/*< free allocated storage >*/
{
  free (*ltf); free (ltf);
  free (*rtf); free (rtf);
  free (*ltb); free (ltb);
  free (*rtb); free (rtb);
}

void lrexp(sf_complex *img, sf_complex *dat, bool adj, sf_complex **lt, sf_complex **rt, geopar geop, sf_complex ***wvfld)
/*< zero-offset exploding reflector modeling/migration >*/
{
    int it, nt, ix, nx, nx2, iz, nz, nz2, nzx2, gpz, wfnt, wfit, snap;
    int im, i, j, m2, ik, nk, pad1;
    float dt, dx, dz, ox;
    sf_complex *curr, **wave, *cwave, *cwavem, c;
    sf_complex *currm;
    bool verb;

    nx  = geop->nx;
    nz  = geop->nz;
    dx  = geop->dx;
    dz  = geop->dz;
    ox  = geop->ox;
    gpz = geop->gpz;
    nt  = geop->nt;
    dt  = geop->dt;
    snap= geop->snap;
    nzx2= geop->nzx2;
    m2  = geop->m2;
    wfnt= geop->wfnt;
    pad1= geop->pad1;
    verb= geop->verb;

    nk = cfft2_init(pad1,nz,nx,&nz2,&nx2);
    if (nk!=geop->nk) sf_error("nk discrepancy!");

    curr = sf_complexalloc(nzx2);
    cwave  = sf_complexalloc(nk);
    wave = sf_complexalloc2(nzx2,m2);
    if (adj) {
	currm  = sf_complexalloc(nzx2);
	icfft2_allocate(cwave);
    } else {
	cwavem = sf_complexalloc(nk);
	icfft2_allocate(cwavem);
    }

#ifdef _OPENMP
#pragma omp parallel for private(iz)
#endif
    for (iz=0; iz < nzx2; iz++) {
	curr[iz] = sf_cmplx(0.,0.);
    }

    if (adj) { /* migration <- read data */
        if (snap>0) wfit = (int)(nt-1)/snap; // wfnt-1
	/* time stepping */
	for (it=nt-1; it > -1; it--) {
	    if (verb) sf_warning("it=%d;",it);
	
	    /* matrix multiplication */
	    for (im = 0; im < m2; im++) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,i,j) shared(currm,lt,curr)
#endif
		for (ix = 0; ix < nx; ix++) {
		    for (iz=0; iz < nz; iz++) {
			i = iz+ix*nz;  /* original grid */
			j = iz+ix*nz2; /* padded grid */
#ifdef SF_HAS_COMPLEX_H
			currm[j] = conjf(lt[im][i])*curr[j];
#else
			currm[j] = sf_cmul(conjf(lt[im][i]), curr[j]);
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
		    c += wave[im][ik]*conjf(rt[ik][im]);
#else
		    c += sf_cmul(wave[im][ik],conjf(rt[ik][im])); //complex multiplies complex
#endif
		}
		cwave[ik] = c;
	    }

	    icfft2(curr,cwave);

#ifdef _OPENMP
#pragma omp parallel for private(ix)
#endif
	    /* data injection */
	    for (ix=0; ix < nx; ix++) {
		curr[gpz+ix*nz2] += dat[it+ix*nt];
	    }

	    if (snap > 0 && it%snap == 0) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,j)
#endif
		for ( ix = 0; ix < nx; ix++) {
		    for ( iz = 0; iz<nz; iz++ ) { 
			j = iz+ix*nz2; /* padded grid */
			wvfld[wfit][ix][iz] = curr[j];
		    }
		}
		if (snap>0) wfit--;
	    }
	} /*time iteration*/
	/*generate image*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		img[iz+ix*nz] = curr[iz+ix*nz2];
	    }
	}
    } else { /* modeling -> write data */
	/*exploding reflector*/
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz)
#endif
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		curr[iz+ix*nz2]=img[iz+ix*nz];
	    }
	}
	if (snap>0) wfit = 0;
	/* time stepping */
	for (it=0; it < nt; it++) {
	    if (verb) sf_warning("it=%d;",it);
#ifdef _OPENMP
#pragma omp parallel for private(ix)
#endif
	    /* record data */
	    for (ix=0; ix < nx; ix++) {
		dat[it+ix*nt] = curr[gpz+ix*nz2];
	    }
	    /* matrix multiplication */
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
	    for (ix = 0; ix < nx; ix++) {
		for (iz=0; iz < nz; iz++) {
		    i = iz+ix*nz;  /* original grid */
		    j = iz+ix*nz2; /* padded grid */
		    
		    c = sf_cmplx(0.,0.); /* initialize */
		    
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
	    if (snap > 0 && it%snap == 0) {
#ifdef _OPENMP
#pragma omp parallel for private(ix,iz,j)
#endif
		for ( ix = 0; ix < nx; ix++) {
		    for ( iz = 0; iz<nz; iz++ ) { 
			j = iz+ix*nz2; /* padded grid */
			wvfld[wfit][ix][iz] = curr[j];
		    }
		}
		if (snap>0) wfit++;
	    }
	}
    }
    if (verb) sf_warning(".");

    cfft2_finalize();
}

/* just for debugging!!! */
//void lrexp_op(int nx, sf_complex* x, sf_complex* y, void* mat, sf_complex ***wvfld, bool adj)
///*< lowrank onestep exploding reflector linear operator (square matrix, img to img) >*/
//{
//  geopar geop;
//  sf_complex *dat, *img;
//  geop = (geopar) mat;
//  img = x;
//  dat = y;
//  lrexp(img, dat, adj, ltf, rtf, geop, wvfld);
//}

void lrexp_op(int nx, const sf_complex* x, sf_complex* y, void* mat)
/*< lowrank onestep exploding reflector linear operator (square matrix, img to img) >*/
{
  geopar geop;
  sf_complex *dat;
  int ntx;
  /*converting the pointer to correct type*/
  geop = (geopar) mat;
  ntx = geop->nt * geop->nx;
  /* image is x, dat is y */
  dat = sf_complexalloc(ntx);
  /* to disable saving the wavefield */
  if (geop->snap > 0)
    sf_error("snap cannot be bigger than zero in this program!");

  /* exploding reflector modeling */
  lrexp(x, dat, false, ltf, rtf, geop, NULL);

  /* exploding reflector migration */
  lrexp(y, dat, true, ltb, rtb, geop, NULL);

}
