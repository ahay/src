/*  Pseudo-spectral method for 2-D elastic wave modeling */
/*
  Copyright (C) 2014 Institute of Geology and Geophysics, Chinese Academy of Sciences (Jun Yan) 
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
#include "psstep2e.h"
#include "fft2.h"
//#include "fft2d_JYAN.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static sf_complex **ukxx, **ukzz, **ukxz, **wkxx, **wkzz, **wkxz, **uktmp, **wktmp;
static float **uxx, **uzz, **uxz, **wxx, **wzz, **wxz,**uutmp, **wwtmp;
static float dkx, dkz;
static int nkx, nkz, opt;

void psstep2e_init(int nz, int nx /*model size*/,
                  float dz, float dx /*model grid*/,
				  int opt1)
/*< initialization >*/
{
//	nkx = fft_nk(nx, opt1);
//	nkz = fft_nk(nz, opt1);
	nkx = kiss_fft_next_fast_size(nx);
	nkz = kiss_fft_next_fast_size(nz);

    ukxx = sf_complexalloc2(nkz,nkx);
    ukzz = sf_complexalloc2(nkz,nkx);
    ukxz = sf_complexalloc2(nkz,nkx);
    wkxx = sf_complexalloc2(nkz,nkx);
    wkzz = sf_complexalloc2(nkz,nkx);
    wkxz = sf_complexalloc2(nkz,nkx);
    uktmp = sf_complexalloc2(nkz,nkx);
    wktmp = sf_complexalloc2(nkz,nkx);

    uxx = sf_floatalloc2(nkz,nkx);
    uzz = sf_floatalloc2(nkz,nkx);
    uxz = sf_floatalloc2(nkz,nkx);
    wxx = sf_floatalloc2(nkz,nkx);
    wzz = sf_floatalloc2(nkz,nkx);
    wxz = sf_floatalloc2(nkz,nkx);
    uutmp = sf_floatalloc2(nkz,nkx);
    wwtmp = sf_floatalloc2(nkz,nkx);

    dkx = 1./(nkx*dx);
    dkz = 1./(nkz*dz);
	opt = opt1;
}

void psstep2e_close(void)
/*< free memory allocation>*/
{
    free(*ukxx);     
    free(ukxx);     
    free(*ukzz);     
    free(ukzz);     
    free(*ukxz);     
    free(ukxz);     
    free(*wkxx);     
    free(wkxx);     
    free(*wkzz);     
    free(wkzz);     
    free(*wkxz);     
    free(wkxz);     

    free(*uktmp);     
    free(uktmp);  

    free(*wktmp);     
    free(wktmp); 

    free(*uxx);     
    free(uxx);     
    free(*uzz);     
    free(uzz);     
    free(*uxz);     
    free(uxz);     
    free(*wxx);     
    free(wxx);     
    free(*wzz);     
    free(wzz);     
    free(*wxz);     
    free(wxz);     

    free(*uutmp);     
    free(uutmp);
    free(*wwtmp);     
    free(wwtmp);


}


void psstep2e(float **upold /*previous step*/,
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
             float vp0 /*reference vel*/,
			 float vs0,
             float **vp /*reference vel*/,
             float **vs/*reference vel*/,
             float dt /*time step size*/)
/*< PS step>*/
{
    int ix, ikx, ikz, iz;
    float kx, kz, tmpdt, pi=SF_PI;
	float kx0,kz0;

	kx0 =-0.5/dx;
	kz0 =-0.5/dz;

    for (ix=0; ix < nkx; ix++){ 
        for (iz=0; iz < nkz; iz++){ 

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

	/* Transform the input data into wavenumber space domain */
	fft2(uutmp[0], uktmp[0]);
	fft2(wwtmp[0], wktmp[0]);
	
	for (ikx=0; ikx < nkx; ikx++) {
		for (ikz=0; ikz < nkz; ikz++) {
			ukxx[ikx][ikz] = uktmp[ikx][ikz];
			ukzz[ikx][ikz] = uktmp[ikx][ikz];
			ukxz[ikx][ikz] = uktmp[ikx][ikz];
			wkxx[ikx][ikz] = wktmp[ikx][ikz];
			wkzz[ikx][ikz] = wktmp[ikx][ikz];
			wkxz[ikx][ikz] = wktmp[ikx][ikz];
			}
	}
	

	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			tmpdt = -1.0*kx*kx;
#ifdef SF_HAS_COMPLEX_H
			ukxx[ikx][ikz] *= tmpdt; 
#else
			ukxx[ikx][ikz] = sf_crmul(ukxx[ikx][ikz],tmpdt);
#endif
		}
	}
	/* Inverse FFT*/
	ifft2(uxx[0], ukxx[0]);
	
    /* compute u2(kx,kz) */
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			tmpdt = -1.0*kz*kz;
#ifdef SF_HAS_COMPLEX_H
			ukzz[ikx][ikz] *= tmpdt; 
#else
			ukzz[ikx][ikz] = sf_crmul(ukzz[ikx][ikz],tmpdt);
#endif
		}
	}
	ifft2(uzz[0], ukzz[0]);

   
   	/* compute u3(kx,kz) */
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			tmpdt = -1.0*kx*kz;
#ifdef SF_HAS_COMPLEX_H
			ukxz[ikx][ikz] *= tmpdt; 
#else
			ukxz[ikx][ikz] = sf_crmul(ukxz[ikx][ikz],tmpdt);
#endif
		}
	}
	/* Inverse FFT*/
	ifft2(uxz[0], ukxz[0]);


    /* compute u4(kx,kz) */
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			tmpdt = -1.0*kx*kx;
#ifdef SF_HAS_COMPLEX_H
			wkxx[ikx][ikz] *= tmpdt; 
#else
			wkxx[ikx][ikz] = sf_crmul(wkxx[ikx][ikz],tmpdt);
#endif
		}
	}
	/* Inverse FFT*/
	ifft2(wxx[0], wkxx[0]);


    /* compute u5(kx,kz) */
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			tmpdt = -1.0*kz*kz;
#ifdef SF_HAS_COMPLEX_H
			wkzz[ikx][ikz] *= tmpdt;
#else
			wkzz[ikx][ikz] = sf_crmul(wkzz[ikx][ikz],tmpdt);
#endif
		}
	}
	/* Inverse FFT*/
	ifft2(wzz[0], wkzz[0]);


     /* compute u6(kx,kz) */
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			tmpdt = -1.0*kx*kz;
#ifdef SF_HAS_COMPLEX_H
			wkxz[ikx][ikz] *= tmpdt;
#else
			wkxz[ikx][ikz] = sf_crmul(wkxz[ikx][ikz],tmpdt);
#endif
		}
	}
	/* Inverse FFT*/
	ifft2(wxz[0], wkxz[0]);



#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix=0; ix < nx; ix++) {  
        for (iz=0; iz < nz; iz++) {  
            upold[ix][iz] = (dt*dt*vp[ix][iz]*vp[ix][iz])*(uxx[ix][iz]+wxz[ix][iz]) + 2.0*upcur[ix][iz]-upold[ix][iz];
			wpold[ix][iz]  = (dt*dt*vp[ix][iz]*vp[ix][iz])*(uxz[ix][iz]+wzz[ix][iz]) + 2.0*wpcur[ix][iz]-wpold[ix][iz];
  			usold[ix][iz]  = (dt*dt*vs[ix][iz]*vs[ix][iz])*(uzz[ix][iz]-wxz[ix][iz]) + 2.0*uscur[ix][iz]-usold[ix][iz];
  			wsold[ix][iz]  = (dt*dt*vs[ix][iz]*vs[ix][iz])*(wxx[ix][iz]-uxz[ix][iz]) + 2.0*wscur[ix][iz]-wsold[ix][iz];
        }
	}  

}
