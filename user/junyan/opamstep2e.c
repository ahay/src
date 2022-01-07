/* 2-D opam for elastic wave modeling and vector field decompostion */
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
#include "opamstep2e.h"
#include "fft2.h"
//#include "fft2d_JYAN.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static sf_complex *cwaveu, *cwavew, *cwavem;
static float **uxx, **uzz, **uxzp, **uxzs, **wxx, **wzz, **wxzp, **wxzs, **uutmp, **wwtmp, **wave;
float **Cuxxl, **Cuxxr, **Cwxzpl, **Cwxzpr,  **Cuxzpl, **Cuxzpr, **Cwzzl, **Cwzzr, **Cuzzl, **Cuzzr, **Cwxzsl, **Cwxzsr, **Cwxxl, **Cwxxr, **Cuxzsl, **Cuxzsr;
static int nx, nz, nk, nx2, nz2, M1, M2, M3, M4, M5, M6, M7, M8,  nxzb2;

void lowrankelastic_init2(int nzb, int nxb, int nkxz, int nkzz, int nkxx, int m1, int m2, int m3, int m4, int m5, int m6, int m7, int m8, int nxzb, float **cuxxl, float **cuxxr, float **cwxzpl, float **cwxzpr, float **cuxzpl, float **cuxzpr, float **cwzzl, float **cwzzr, float **cuzzl, float **cuzzr, float **cwxzsl, float **cwxzsr, float **cwxxl, float **cwxxr, float **cuxzsl, float **cuxzsr)
/*< 2D lowrank elastic initiazation >*/
{
    int ik;

    nz = nzb;
    nx = nxb;
    nk = nkxz;

    nz2 = nkzz;
    nx2 = nkxx;

    M1 = m1;
    M2 = m2;
    M3 = m3;
    M4 = m4;
    M5 = m5;
    M6 = m6;
    M7 = m7;
    M8 = m8;

    nxzb2 = nxzb;

    Cuxxl = cuxxl;
    Cuxxr = cuxxr;
    Cwxzpl = cwxzpl;
    Cwxzpr = cwxzpr;
    Cuxzpl = cuxzpl;
    Cuxzpr = cuxzpr;
    Cwzzl = cwzzl;
    Cwzzr = cwzzr;
    Cuzzl = cuzzl;
    Cuzzr = cuzzr;
    Cwxzsl = cwxzsl;
    Cwxzsr = cwxzsr;
    Cwxxl = cwxxl;
    Cwxxr = cwxxr;
    Cuxzsl = cuxzsl;
    Cuxzsr = cuxzsr;

    uxx = sf_floatalloc2(nz,nx);
    uzz = sf_floatalloc2(nz,nx);
    uxzp = sf_floatalloc2(nz,nx);
    uxzs = sf_floatalloc2(nz,nx);
    wxx = sf_floatalloc2(nz,nx);
    wzz = sf_floatalloc2(nz,nx);
    wxzp = sf_floatalloc2(nz,nx);
    wxzs = sf_floatalloc2(nz,nx);


    uutmp = sf_floatalloc2(nz2,nx2);
    wwtmp = sf_floatalloc2(nz2,nx2);

    cwaveu  = sf_complexalloc(nk);
    cwavew  = sf_complexalloc(nk);
    cwavem = sf_complexalloc(nk);

//	wave = sf_floatalloc2(nxzb2,m2);

    for (ik = 0; ik < nk; ik++) {
	cwaveu[ik] = sf_cmplx(0.0, 0.0);
	cwavew[ik] = sf_cmplx(0.0, 0.0);
	cwavem[ik] = sf_cmplx(0.0, 0.0);
    }
//	for (ix = 0; ix < m2; ix++) {
//	    for (iz = 0; iz < nxzb2; iz++) {
//			wave[ix][iz] = 0.0;
//		}
//	}
}


void lowrankelastic_close2()
/*< free the work space >*/
{

    free(*uxx);     
    free(uxx);     
    free(*uzz);     
    free(uzz);     
    free(*uxzp);     
    free(uxzp);     
    free(*uxzs);     
    free(uxzs);     
    free(*wxx);     
    free(wxx);     
    free(*wzz);     
    free(wzz);     
    free(*wxzp);     
    free(wxzp);     
    free(*wxzs);     
    free(wxzs);     

    free(*uutmp);     
    free(uutmp);
    free(*wwtmp);     
    free(wwtmp);

    free(*wave);     
    free(wave);
    free(cwaveu);
    free(cwavew);
    free(cwavem);


}

void lowrankelastic_comp2(float **uu, sf_complex *cwave, sf_complex *cwavem, float **lt, float **rt, int mm)
/*< K space mix-matrix approximated by lowrank >*/
{
    int i, j, ix, iz, ik, im;

    ifft2_allocate(cwavem);
	
    if (wave != NULL){
	free(*wave);     
	free(wave);
    }
    wave = sf_floatalloc2(nxzb2,mm);
	
    for (ix = 0; ix < mm; ix++) {
	for (iz = 0; iz < nxzb2; iz++) {
	    wave[ix][iz] = 0.0;
	}
    }

    for (im = 0; im < mm; im++) {
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

	    for (im = 0; im < mm; im++) {
		uu[ix][iz] += lt[im][i]*wave[im][j];
	    }
	}
    }
}  


void opamstep2e(float **upold /*previous step*/,
		float **upcur /*current step*/,
		float **wpold,
		float **wpcur,
		float **usold,
		float **uscur,
		float **wsold,
		float **wscur,
		float **uu,
		float **ww,
		int nz, int nx /*model size*/,
		float dz, float dx /*model grid */,
		float **vp /*p vel*/,
		float **vs /*s vel*/,
		float dt /*time step size*/)
/*< Optimized PAM step>*/
{
    int ix, iz;

    for (ix=0; ix < nx; ix++){ 
        for (iz=0; iz < nz; iz++){ 
            uxx[ix][iz] = 0.; 
            wxzp[ix][iz] = 0.; 
            uxzp[ix][iz] = 0.; 
            wzz[ix][iz] = 0.; 
            uzz[ix][iz] = 0.; 
            wxzs[ix][iz] = 0.; 
            wxx[ix][iz] = 0.;
            uxzs[ix][iz] = 0.;
	}  
    }

    for (ix=0; ix < nx2; ix++){ 
        for (iz=0; iz < nz2; iz++){ 

	    uutmp[ix][iz] = 0.;
	    wwtmp[ix][iz] = 0.;
	}  
    }
	
    for (ix=0; ix < nx; ix++){ 
        for (iz=0; iz < nz; iz++){ 
	    uutmp[ix][iz] = uu[ix][iz];	
	    wwtmp[ix][iz] = ww[ix][iz];	
	}  
    }
   
    //	nkxz=fft2_init(true, 1, nz, nx, &nkzz, &nkxx);

    /* Just do one FFT and severn IFFT in one time step */
    fft2(uutmp[0], cwaveu);
    fft2(wwtmp[0], cwavew);
//fprintf(stderr, "I'm here...\n");
    lowrankelastic_comp2(uxx, cwaveu, cwavem, Cuxxl, Cuxxr, M1);
    lowrankelastic_comp2(wxzp, cwavew, cwavem, Cwxzpl, Cwxzpr, M2);
    lowrankelastic_comp2(uxzp, cwaveu, cwavem, Cuxzpl, Cuxzpr, M3);
    lowrankelastic_comp2(wzz, cwavew, cwavem, Cwzzl, Cwzzr, M4);
    lowrankelastic_comp2(uzz, cwaveu, cwavem, Cuzzl, Cuzzr, M5);
    lowrankelastic_comp2(wxzs, cwavew, cwavem, Cwxzsl, Cwxzsr, M6);
    lowrankelastic_comp2(wxx, cwavew, cwavem, Cwxxl, Cwxxr, M7);
    lowrankelastic_comp2(uxzs, cwaveu, cwavem, Cuxzsl, Cuxzsr, M8);

//fprintf(stderr, "I'm here...too\n");
	
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix=0; ix < nx; ix++) {  
        for (iz=0; iz < nz; iz++) {  
            upold[ix][iz] = (dt*dt*vp[ix][iz]*vp[ix][iz])*(uxx[ix][iz]+wxzp[ix][iz]) + 2.0*upcur[ix][iz]-upold[ix][iz];
	    wpold[ix][iz]  = (dt*dt*vp[ix][iz]*vp[ix][iz])*(uxzp[ix][iz]+wzz[ix][iz]) + 2.0*wpcur[ix][iz]-wpold[ix][iz];
	    usold[ix][iz]  = (dt*dt*vs[ix][iz]*vs[ix][iz])*(uzz[ix][iz]-wxzs[ix][iz]) + 2.0*uscur[ix][iz]-usold[ix][iz];
	    wsold[ix][iz]  = (dt*dt*vs[ix][iz]*vs[ix][iz])*(wxx[ix][iz]-uxzs[ix][iz]) + 2.0*wscur[ix][iz]-wsold[ix][iz];
        }
    }  

//    for (ix=0; ix < nx; ix++) {  
//        for (iz=0; iz < nz; iz++) {  
//fprintf(stderr, "%f\n", upold[ix][iz]);
//		}
//	}

}
