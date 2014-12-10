/*  2-D elasitc wave modeling and vector field decompostion using pseudo-analytical method */
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
#include "pamstep2e.h"
#include "fft2.h"
//#include "fft2d_JYAN.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static sf_complex **ukxx, **ukzz, **ukxzp, **ukxzs, **wkxx, **wkzz, **wkxzp, **wkxzs, **uktmp, **wktmp;
static float **uxx, **uzz, **uxzp, **uxzs, **wxx, **wzz, **wxzp, **wxzs, **uutmp, **wwtmp;
static float dkx, dkz;
static int nkx, nkz, opt;

void pamstep2e_init(int nz, int nx /*model size*/,
                  float dz, float dx /*model grid*/,
				  int opt1)
/*< initialization >*/
{
//	nkx = fft_nk(nx, opt1);
//	nkz = fft_nk(nz, opt1);
	nkx = opt1? kiss_fft_next_fast_size(nx):nx;
	nkz = opt1? kiss_fft_next_fast_size(nz):nz;

    ukxx = sf_complexalloc2(nkz,nkx);
    ukzz = sf_complexalloc2(nkz,nkx);
    ukxzp = sf_complexalloc2(nkz,nkx);
    ukxzs = sf_complexalloc2(nkz,nkx);
    wkxx = sf_complexalloc2(nkz,nkx);
    wkzz = sf_complexalloc2(nkz,nkx);
    wkxzp = sf_complexalloc2(nkz,nkx);
    wkxzs = sf_complexalloc2(nkz,nkx);
    uktmp = sf_complexalloc2(nkz,nkx);
    wktmp = sf_complexalloc2(nkz,nkx);

    uxx = sf_floatalloc2(nkz,nkx);
    uzz = sf_floatalloc2(nkz,nkx);
    uxzp = sf_floatalloc2(nkz,nkx);
    uxzs = sf_floatalloc2(nkz,nkx);
    wxx = sf_floatalloc2(nkz,nkx);
    wzz = sf_floatalloc2(nkz,nkx);
    wxzp = sf_floatalloc2(nkz,nkx);
    wxzs = sf_floatalloc2(nkz,nkx);
    uutmp = sf_floatalloc2(nkz,nkx);
    wwtmp = sf_floatalloc2(nkz,nkx);

    dkx = 1./(nkx*dx);
    dkz = 1./(nkz*dz);
	opt = opt1;
}

void pamstep2e_close(void)
/*< free memory allocation>*/
{
    free(*ukxx);     
    free(ukxx);     
    free(*ukzz);     
    free(ukzz);     
    free(*ukxzp);     
    free(ukxzp);     
    free(*ukxzs);     
    free(ukxzs);     
    free(*wkxx);     
    free(wkxx);     
    free(*wkzz);     
    free(wkzz);     
    free(*wkxzp);     
    free(wkxzp);     
    free(*wkxzs);     
    free(wkxzs);     

    free(*uktmp);     
    free(uktmp);  

    free(*wktmp);     
    free(wktmp); 

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


}


void pamstep2e(float **upold /*previous step*/,
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
/*< PAM step>*/
{
    int ix, ikx, ikz, iz;
    float kx, kz, tmpdt, pi=SF_PI;
	float kx0,kz0,kxz2;
	int nkxx, nkzz, nkxz;

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
   
//	nkxz=fft2_init(true, opt, 1,  nz, nx, &nkzz, &nkxx);

	/* Just do one FFT and severn IFFT in one time step */

	/* Transform the input data into wavenumber space domain */
	fft2(uutmp[0], uktmp[0]);
	fft2(wwtmp[0], wktmp[0]);
	
	for (ikx=0; ikx < nkx; ikx++) {
		for (ikz=0; ikz < nkz; ikz++) {
			ukxx[ikx][ikz] = uktmp[ikx][ikz];
			ukzz[ikx][ikz] = uktmp[ikx][ikz];
			ukxzp[ikx][ikz] = uktmp[ikx][ikz];
			ukxzs[ikx][ikz] = uktmp[ikx][ikz];
			wkxx[ikx][ikz] = wktmp[ikx][ikz];
			wkzz[ikx][ikz] = wktmp[ikx][ikz];
			wkxzp[ikx][ikz] = wktmp[ikx][ikz];
			wkxzs[ikx][ikz] = wktmp[ikx][ikz];
			}
	}
	

	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxz2 = kx*kx + kz*kz;
			if (!kxz2)
			kxz2 +=0.000001;
			tmpdt = 1.0*kx*kx*2.0*(cosf(vp0*sqrtf(kx*kx+kz*kz)*dt)-1.0)/(vp0*vp0*dt*dt*kxz2);
			ukxx[ikx][ikz] = sf_crmul(ukxx[ikx][ikz],tmpdt);
		}
	}
	/* Inverse FFT*/
	ifft2(uxx[0], ukxx[0]);
	
    /* compute u2(kx,kz) */
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxz2 = kx*kx + kz*kz;
			if (!kxz2)
			kxz2 +=0.000001;
			tmpdt = 1.0*kz*kz*2.0*(cosf(vs0*sqrtf(kx*kx+kz*kz)*dt)-1.0)/(vs0*vs0*dt*dt*kxz2);
			ukzz[ikx][ikz] = sf_crmul(ukzz[ikx][ikz],tmpdt);
		}
	}
	ifft2(uzz[0], ukzz[0]);

   
   	/* compute u3(kx,kz) */
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxz2 = kx*kx + kz*kz;
			if (!kxz2)
			kxz2 +=0.000001;
			tmpdt = 1.0*kx*kz*2.0*(cosf(vp0*sqrtf(kx*kx+kz*kz)*dt)-1.0)/(vp0*vp0*dt*dt*kxz2);
			ukxzp[ikx][ikz] = sf_crmul(ukxzp[ikx][ikz],tmpdt);
		}
	}
	/* Inverse FFT*/
	ifft2(uxzp[0], ukxzp[0]);

   	/* compute u3(kx,kz) */
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxz2 = kx*kx + kz*kz;
			if (!kxz2)
			kxz2 +=0.000001;
			tmpdt = 1.0*kx*kz*2.0*(cosf(vs0*sqrtf(kx*kx+kz*kz)*dt)-1.0)/(vs0*vs0*dt*dt*kxz2);
			ukxzs[ikx][ikz] = sf_crmul(ukxzs[ikx][ikz],tmpdt);
		}
	}
	/* Inverse FFT*/
	ifft2(uxzs[0], ukxzs[0]);

    /* compute u4(kx,kz) */
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxz2 = kx*kx + kz*kz;
			if (!kxz2)
			kxz2 +=0.000001;
			tmpdt = 1.0*kx*kx*2.0*(cosf(vs0*sqrtf(kx*kx+kz*kz)*dt)-1.0)/(vs0*vs0*dt*dt*kxz2);
			wkxx[ikx][ikz] = sf_crmul(wkxx[ikx][ikz],tmpdt);
		}
	}
	/* Inverse FFT*/
	ifft2(wxx[0], wkxx[0]);


    /* compute u5(kx,kz) */
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxz2 = kx*kx + kz*kz;
			if (!kxz2)
			kxz2 +=0.000001;
			tmpdt = 1.0*kz*kz*2.0*(cosf(vp0*sqrtf(kx*kx+kz*kz)*dt)-1.0)/(vp0*vp0*dt*dt*kxz2);
			wkzz[ikx][ikz] = sf_crmul(wkzz[ikx][ikz],tmpdt);
		}
	}
	/* Inverse FFT*/
	ifft2(wzz[0], wkzz[0]);


     /* compute u6(kx,kz) */
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxz2 = kx*kx + kz*kz;
			if (!kxz2)
			kxz2 +=0.000001;
			tmpdt = 1.0*kx*kz*2.0*(cosf(vp0*sqrtf(kx*kx+kz*kz)*dt)-1.0)/(vp0*vp0*dt*dt*kxz2);
			wkxzp[ikx][ikz] = sf_crmul(wkxzp[ikx][ikz],tmpdt);
		}
	}
	/* Inverse FFT*/
	ifft2(wxzp[0], wkxzp[0]);


     /* compute u6(kx,kz) */
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxz2 = kx*kx + kz*kz;
			if (!kxz2)
			kxz2 +=0.000001;
			tmpdt = 1.0*kx*kz*2.0*(cosf(vs0*sqrtf(kx*kx+kz*kz)*dt)-1.0)/(vs0*vs0*dt*dt*kxz2);
			wkxzs[ikx][ikz] = sf_crmul(wkxzs[ikx][ikz],tmpdt);
		}
	}
	/* Inverse FFT*/
	ifft2(wxzs[0], wkxzs[0]);

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

}
