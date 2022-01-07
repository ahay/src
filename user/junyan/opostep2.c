/* Modeling of pure acoustic wave in 2-D transversely isotropic meida using optimized pseudo-Laplacian operator */
/*
  Copyright (C)2014 Institute of Geology and Geophysics, Chinese Academy of Sciences (Jun Yan) 
  2009 University of Texas at Austin
  
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
#include "opostep2.h"
#include "fft2.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static sf_complex *cwave, *cwavem;
static float **uu1, **uu2, **uu3, **uu4, **uu5, **uu6, **uu7, **curtmp, **wave;
static float **Cxxl, **Cxxr, **Czzl, **Czzr, **Cxxxxl, **Cxxxxr, **Czzzzl, **Czzzzr, **Cxzzzl,  **Cxzzzr, **Cxxxzl, **Cxxxzr, **Cxxzzl, **Cxxzzr;
static int nx, nz, nk, nx2, nz2, m2, nxzb2;


void lowrank_init2(int nzb, int nxb, int nkxz, int nkzz, int nkxx, int m, int nxzb, float **cxxl, float **cxxr, float **czzl, float **czzr, float **cxxxxl, float **cxxxxr, float **czzzzl, float **czzzzr, float **cxzzzl, float **cxzzzr, float **cxxxzl, float **cxxxzr, float **cxxzzl, float **cxxzzr)
/*< 2D lowrank initiazation >*/
{
    int ix, iz, ik;

    nz = nzb;
    nx = nxb;
    nk = nkxz;

    nz2 = nkzz;
    nx2 = nkxx;

    m2 = m;
    nxzb2 = nxzb;

    Cxxl = cxxl;
    Cxxr = cxxr;
    Czzl = czzl;
    Czzr = czzr;
    Cxxxxl = cxxxxl;
    Cxxxxr = cxxxxr;
    Czzzzl = czzzzl;
    Czzzzr = czzzzr;
    Cxzzzl = cxzzzl;
    Cxzzzr = cxzzzr;
    Cxxxzl = cxxxzl;
    Cxxxzr = cxxxzr;
    Cxxzzl = cxxzzl;
    Cxxzzr = cxxzzr;

    uu1 = sf_floatalloc2(nz,nx);
    uu2 = sf_floatalloc2(nz,nx);
    uu3 = sf_floatalloc2(nz,nx);
    uu4 = sf_floatalloc2(nz,nx);
    uu5 = sf_floatalloc2(nz,nx);
    uu6 = sf_floatalloc2(nz,nx);
    uu7 = sf_floatalloc2(nz,nx);
 
    curtmp = sf_floatalloc2(nz2,nx2);

    cwave  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);
    wave = sf_floatalloc2(nxzb2,m2);

    for (ik = 0; ik < nk; ik++) {
	cwave[ik] = sf_cmplx(0.0, 0.0);
	cwavem[ik] = sf_cmplx(0.0, 0.0);
    }
    for (ix = 0; ix < m2; ix++) {
	for (iz = 0; iz < nxzb2; iz++) {
	    wave[ix][iz] = 0.0;
	}
    }

}

void lowrank_close2()
/*< free the work space for opo >*/
{

}


void lowrank_comp2(float **uu, sf_complex *cwave, sf_complex *cwavem, float **lt, float **rt)
/*< K space mix-matrix approximated by lowrank >*/
{
    int i, j, ix, iz, ik, im;

    ifft2_allocate(cwavem);

    for (im = 0; im < m2; im++) {
	for (ik = 0; ik < nk; ik++) {
#ifdef SF_HAS_COMPLEX_H
	    cwavem[ik] = cwave[ik]*rt[ik][im];
#else
	    cwavem[ik] = sf_crmul(cwave[ik],rt[ik][im]);
#endif
	}
	ifft2(wave[im],cwavem);
    }

    for (ix = 0; ix < nx; ix++) {
	for (iz=0; iz < nz; iz++) {
	    i = iz+ix*nz;  /* original grid */
	    j = iz+ix*nz2; /* padded grid */

	    for (im = 0; im < m2; im++) {
		uu[ix][iz] += lt[im][i]*wave[im][j];
	    }
	}
    }
}  

void opostep2(float **old /*previous step*/,
	      float **cur /*current step*/,
	      int nz, int nx /*model size*/,
	      float dz, float dx /*model grid */,
	      float v0 /*reference vel*/,
	      float **v /*reference vel*/,
	      float **sigma /*reference vel*/,
	      float **delta /*reference vel*/,
	      float **seta /*reference vel*/,
	      float dt /*time step size*/)
/*< Optimized pseudo-Laplacian operator>*/
{
    int ix, iz;


    for (ix=0; ix < nx; ix++){ 
        for (iz=0; iz < nz; iz++){ 
            uu1[ix][iz] = 0.; 
            uu2[ix][iz] = 0.; 
            uu3[ix][iz] = 0.; 
            uu4[ix][iz] = 0.; 
            uu5[ix][iz] = 0.; 
            uu6[ix][iz] = 0.; 
            uu7[ix][iz] = 0.;
	}  
    }
	
    for (ix=0; ix < nx2; ix++){ 
        for (iz=0; iz < nz2; iz++){ 
	    curtmp[ix][iz] = 0.;
	}
    }

    for (ix=0; ix < nx; ix++){ 
        for (iz=0; iz < nz; iz++){ 
	    curtmp[ix][iz] = cur[ix][iz];	
	}  
    }

    //	nkxz=fft2_init(true, 1, nz, nx, &nkzz, &nkxx);

    /* Just do one FFT and severn IFFT in one time step */
    fft2(curtmp[0], cwave);

    lowrank_comp2(uu1, cwave, cwavem, Cxxl, Cxxr);
    lowrank_comp2(uu2, cwave, cwavem, Czzl, Czzr);
    lowrank_comp2(uu3, cwave, cwavem, Cxxxxl, Cxxxxr);
    lowrank_comp2(uu4, cwave, cwavem, Czzzzl, Czzzzr);
    lowrank_comp2(uu5, cwave, cwavem, Cxzzzl, Cxzzzr);
    lowrank_comp2(uu6, cwave, cwavem, Cxxxzl, Cxxxzr);
    lowrank_comp2(uu7, cwave, cwavem, Cxxzzl, Cxxzzr);
	

#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix=0; ix < nx; ix++) {  
        for (iz=0; iz < nz; iz++) {  
            old[ix][iz]  = (-1.0*dt*dt*v[ix][iz]*v[ix][iz])* ( uu1[ix][iz] + uu2[ix][iz]
							       + (2.0*sigma[ix][iz]*cosf(seta[ix][iz])*cosf(seta[ix][iz])*cosf(seta[ix][iz])*cosf(seta[ix][iz]) + 2.0*delta[ix][iz]*sinf(seta[ix][iz])*sinf(seta[ix][iz])*cosf(seta[ix][iz])*cosf(seta[ix][iz]))*uu3[ix][iz]
							       + (2.0*sigma[ix][iz]*sinf(seta[ix][iz])*sinf(seta[ix][iz])*sinf(seta[ix][iz])*sinf(seta[ix][iz]) + 2.0*delta[ix][iz]*sinf(seta[ix][iz])*sinf(seta[ix][iz])*cosf(seta[ix][iz])*cosf(seta[ix][iz]))*uu4[ix][iz]
							       + (-4.0*sigma[ix][iz]*sinf(2.0*seta[ix][iz])*cosf(seta[ix][iz])*cosf(seta[ix][iz]) + delta[ix][iz]*sinf(4.0*seta[ix][iz]))*uu6[ix][iz]
							       + (-4.0*sigma[ix][iz]*sinf(2.0*seta[ix][iz])*sinf(seta[ix][iz])*sinf(seta[ix][iz]) - delta[ix][iz]*sinf(4.0*seta[ix][iz]))*uu5[ix][iz]
							       + (3.0*sigma[ix][iz]*sinf(2.0*seta[ix][iz])*sinf(2.0*seta[ix][iz]) - delta[ix][iz]*sinf(2.0*seta[ix][iz])*sinf(2.0*seta[ix][iz]) + 2.0*delta[ix][iz]*cosf(2.0*seta[ix][iz])*cosf(2.0*seta[ix][iz]))*uu7[ix][iz] )			   + 2.0*cur[ix][iz]-old[ix][iz];
        }
    }  

}

