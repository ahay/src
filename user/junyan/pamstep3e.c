/*  3-D elasitc wave modeling and vector field decompostion using pseudo-analytical method */
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
#include "pamstep3e.h"
#include "fft3.h"
//#include "fft2d_JYAN.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static sf_complex ***ukxx, ***ukyy, ***ukzz, ***ukxyp, ***ukxys, ***ukxzp, ***ukxzs, ***vkxx,  ***vkyy,  ***vkzz,  ***vkxyp,  ***vkxys,  ***vkyzp,  ***vkyzs,  ***wkxx, ***wkyy, ***wkzz, ***wkxzp, ***wkxzs,***wkyzp, ***wkyzs, ***uktmp, ***vktmp, ***wktmp;
static float ***uxx, ***uyy, ***uzz, ***uxyp, ***uxys, ***uxzp, ***uxzs,  ***vxx,  ***vyy,  ***vzz,  ***vxyp,  ***vxys,  ***vyzp,  ***vyzs,  ***wxx, ***wyy, ***wzz, ***wxzp, ***wxzs, ***wyzp, ***wyzs,  ***uutmp, ***vvtmp, ***wwtmp;
static float dkx, dky, dkz;
static int nkx, nky, nkz, opt;

void pamstep3e_init(int nz, int nx, int ny /*model size*/,
                  float dz, float dx, float dy /*model grid*/,
				  int opt1)
/*< initialization >*/
{
//	nkx = fft_nk(nx, opt1);
//	nkz = fft_nk(nz, opt1);
	nkx = opt1? kiss_fft_next_fast_size(nx):nx;
	nky = opt1? kiss_fft_next_fast_size(ny):ny;
	nkz = opt1? kiss_fft_next_fast_size(nz):nz;

    ukxx = sf_complexalloc3(nkz,nkx,nky);
    ukyy = sf_complexalloc3(nkz,nkx,nky);
    ukzz = sf_complexalloc3(nkz,nkx,nky);
    ukxzp = sf_complexalloc3(nkz,nkx,nky);
    ukxzs = sf_complexalloc3(nkz,nkx,nky);
    ukxyp = sf_complexalloc3(nkz,nkx,nky);
    ukxys = sf_complexalloc3(nkz,nkx,nky);

    vkxx = sf_complexalloc3(nkz,nkx,nky);
    vkyy = sf_complexalloc3(nkz,nkx,nky);
    vkzz = sf_complexalloc3(nkz,nkx,nky);
    vkxyp = sf_complexalloc3(nkz,nkx,nky);
    vkxys = sf_complexalloc3(nkz,nkx,nky);
    vkyzp = sf_complexalloc3(nkz,nkx,nky);
    vkyzs = sf_complexalloc3(nkz,nkx,nky);

    wkxx = sf_complexalloc3(nkz,nkx,nky);
    wkyy = sf_complexalloc3(nkz,nkx,nky);
    wkzz = sf_complexalloc3(nkz,nkx,nky);
    wkxzp = sf_complexalloc3(nkz,nkx,nky);
    wkxzs = sf_complexalloc3(nkz,nkx,nky);
    wkyzp = sf_complexalloc3(nkz,nkx,nky);
    wkyzs = sf_complexalloc3(nkz,nkx,nky);

    uktmp = sf_complexalloc3(nkz,nkx,nky);
    vktmp = sf_complexalloc3(nkz,nkx,nky);
    wktmp = sf_complexalloc3(nkz,nkx,nky);

    uxx = sf_floatalloc3(nkz,nkx,nky);
    uyy = sf_floatalloc3(nkz,nkx,nky);
    uzz = sf_floatalloc3(nkz,nkx,nky);
    uxyp = sf_floatalloc3(nkz,nkx,nky);
    uxys = sf_floatalloc3(nkz,nkx,nky);
    uxzp = sf_floatalloc3(nkz,nkx,nky);
    uxzs = sf_floatalloc3(nkz,nkx,nky);

    vxx = sf_floatalloc3(nkz,nkx,nky);
    vyy = sf_floatalloc3(nkz,nkx,nky);
    vzz = sf_floatalloc3(nkz,nkx,nky);
    vxyp = sf_floatalloc3(nkz,nkx,nky);
    vxys = sf_floatalloc3(nkz,nkx,nky);
    vyzp = sf_floatalloc3(nkz,nkx,nky);
    vyzs = sf_floatalloc3(nkz,nkx,nky);

    wxx = sf_floatalloc3(nkz,nkx,nky);
    wyy = sf_floatalloc3(nkz,nkx,nky);
    wzz = sf_floatalloc3(nkz,nkx,nky);
    wxzp = sf_floatalloc3(nkz,nkx,nky);
    wxzs = sf_floatalloc3(nkz,nkx,nky);
    wyzp = sf_floatalloc3(nkz,nkx,nky);
    wyzs = sf_floatalloc3(nkz,nkx,nky);

    uutmp = sf_floatalloc3(nkz,nkx,nky);
    vvtmp = sf_floatalloc3(nkz,nkx,nky);
    wwtmp = sf_floatalloc3(nkz,nkx,nky);

    dkx = 1./(nkx*dx);
    dky = 1./(nky*dy);
    dkz = 1./(nkz*dz);
	opt = opt1;
}

void pamstep3e_close(void)
/*< free memory allocation>*/
{
    free(**ukxx);     
    free(**ukyy);     
    free(**ukzz);     
    free(**ukxzp);     
    free(**ukxzs);     
    free(**ukxyp);     
    free(**ukxys);
    free(**uxx);     
    free(**uyy);     
    free(**uzz);     
    free(**uxzp);     
    free(**uxzs);     
    free(**uxyp);     
    free(**uxys);

    free(**vkxx);     
    free(**vkyy);     
    free(**vkzz);     
    free(**vkyzp);     
    free(**vkyzs);     
    free(**vkxyp);     
    free(**vkxys);
    free(**vxx);     
    free(**vyy);     
    free(**vzz);     
    free(**vyzp);     
    free(**vyzs);     
    free(**vxyp);     
    free(**vxys);

    free(**wkxx);     
    free(**wkyy);     
    free(**wkzz);     
    free(**wkxzp);     
    free(**wkxzs);     
    free(**wkyzp);     
    free(**wkyzs);
    free(**wxx);     
    free(**wyy);     
    free(**wzz);     
    free(**wxzp);     
    free(**wxzs);     
    free(**wyzp);     
    free(**wyzs);

    free(*ukxx);     
    free(*ukyy);     
    free(*ukzz);     
    free(*ukxzp);     
    free(*ukxzs);     
    free(*ukxyp);     
    free(*ukxys);
    free(*uxx);     
    free(*uyy);     
    free(*uzz);     
    free(*uxzp);     
    free(*uxzs);     
    free(*uxyp);     
    free(*uxys);

    free(*vkxx);     
    free(*vkyy);     
    free(*vkzz);     
    free(*vkyzp);     
    free(*vkyzs);     
    free(*vkxyp);     
    free(*vkxys);
    free(*vxx);     
    free(*vyy);     
    free(*vzz);     
    free(*vyzp);     
    free(*vyzs);     
    free(*vxyp);     
    free(*vxys);

    free(*wkxx);     
    free(*wkyy);     
    free(*wkzz);     
    free(*wkxzp);     
    free(*wkxzs);     
    free(*wkyzp);     
    free(*wkyzs);
    free(*wxx);     
    free(*wyy);     
    free(*wzz);     
    free(*wxzp);     
    free(*wxzs);     
    free(*wyzp);     
    free(*wyzs);

    free(vkxx);     
    free(vkyy);     
    free(vkzz);     
    free(vkyzp);     
    free(vkyzs);     
    free(vkxyp);     
    free(vkxys);
    free(vxx);     
    free(vyy);     
    free(vzz);     
    free(vyzp);     
    free(vyzs);     
    free(vxyp);     
    free(vxys);

    free(wkxx);     
    free(wkyy);     
    free(wkzz);     
    free(wkxzp);     
    free(wkxzs);     
    free(wkyzp);     
    free(wkyzs);
    free(wxx);     
    free(wyy);     
    free(wzz);     
    free(wxzp);     
    free(wxzs);     
    free(wyzp);     
    free(wyzs);


    free(**uutmp);     
    free(**vvtmp);
    free(**wwtmp);     
    free(**uktmp);
    free(**vktmp);
    free(**wktmp);
    free(*uutmp);     
    free(*vvtmp);
    free(*wwtmp);     
    free(*uktmp);
    free(*vktmp);
    free(*wktmp);
    free(uutmp);     
    free(vvtmp);
    free(wwtmp);     
    free(uktmp);
    free(vktmp);
    free(wktmp);
}


void pamstep3e(float ***upold /*previous step*/,
             float ***upcur /*current step*/,
			 float ***vpold, 
			 float ***vpcur, 
			 float ***wpold,
			 float ***wpcur,
			 float ***usold,
			 float ***uscur,
			 float ***vsold, 
			 float ***vscur, 
			 float ***wsold,
			 float ***wscur,
			 float ***uu,
			 float ***vv, 
			 float ***ww,
             int nz, int nx, int ny /*model size*/,
			 float dz, float dx , float dy/*model grid */,
             float vp0 /*reference vel*/,
			 float vs0,
             float ***vp /*reference vel*/,
             float ***vs/*reference vel*/,
             float dt /*time step size*/)
/*< PAM step>*/
{
    int ix, ikx, iy, iky, ikz, iz;
    float kx, ky, kz, tmpdt, pi=SF_PI;
	float kx0,ky0,kz0,kxyz2;

	kx0 =-0.5/dx;
	ky0 =-0.5/dy;
	kz0 =-0.5/dz;

    for (iy=0; iy < nky; iy++){ 
    for (ix=0; ix < nkx; ix++){ 
        for (iz=0; iz < nkz; iz++){ 
				uutmp[iy][ix][iz] = 0.;
				vvtmp[iy][ix][iz] = 0.;
				wwtmp[iy][ix][iz] = 0.;
         }  
	}
	}
	
//fprintf(stderr, "I'm here1...\n");
    for (iy=0; iy < ny; iy++){ 
    for (ix=0; ix < nx; ix++){ 
        for (iz=0; iz < nz; iz++){ 
		        uutmp[iy][ix][iz] = uu[iy][ix][iz];	
		        vvtmp[iy][ix][iz] = vv[iy][ix][iz];	
		        wwtmp[iy][ix][iz] = ww[iy][ix][iz];	
         }  
	}
	}
   
//	nkxz=fft2_init(true, opt, 1,  nz, nx, &nkzz, &nkxx);

	/* Just do one FFT and severn IFFT in one time step */

	/* Transform the input data into wavenumber space domain */
	fft3(uutmp[0][0], uktmp[0][0]);
	fft3(vvtmp[0][0], vktmp[0][0]);
	fft3(wwtmp[0][0], wktmp[0][0]);
	
	for (iky=0; iky < nky; iky++) {
	for (ikx=0; ikx < nkx; ikx++) {
		for (ikz=0; ikz < nkz; ikz++) {
			ukxx[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukyy[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukzz[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukxzp[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukxzs[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukxyp[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukxys[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			vkxx[iky][ikx][ikz] = vktmp[iky][ikx][ikz];
			vkyy[iky][ikx][ikz] = vktmp[iky][ikx][ikz];
			vkzz[iky][ikx][ikz] = vktmp[iky][ikx][ikz];
			vkyzp[iky][ikx][ikz] = vktmp[iky][ikx][ikz];
			vkyzs[iky][ikx][ikz] = vktmp[iky][ikx][ikz];
			vkxyp[iky][ikx][ikz] = vktmp[iky][ikx][ikz];
			vkxys[iky][ikx][ikz] = vktmp[iky][ikx][ikz];
			wkxx[iky][ikx][ikz] = wktmp[iky][ikx][ikz];
			wkyy[iky][ikx][ikz] = wktmp[iky][ikx][ikz];
			wkzz[iky][ikx][ikz] = wktmp[iky][ikx][ikz];
			wkxzp[iky][ikx][ikz] = wktmp[iky][ikx][ikz];
			wkxzs[iky][ikx][ikz] = wktmp[iky][ikx][ikz];
			wkyzp[iky][ikx][ikz] = wktmp[iky][ikx][ikz];
			wkyzs[iky][ikx][ikz] = wktmp[iky][ikx][ikz];
		}
	}
	}
	

    /* compute uxx */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*kx*kx*2.0*(cosf(vp0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vp0*vp0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			ukxx[iky][ikx][ikz] *= tmpdt; 
#else
			ukxx[iky][ikx][ikz] = sf_crmul(ukxx[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uxx[0][0], ukxx[0][0]);
	

    /* compute uyy */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*ky*ky*2.0*(cosf(vs0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vs0*vs0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			ukyy[iky][ikx][ikz] *= tmpdt; 
#else
			ukyy[iky][ikx][ikz] = sf_crmul(ukyy[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uyy[0][0], ukyy[0][0]);
	

    /* compute uzz */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*kz*kz*2.0*(cosf(vs0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vs0*vs0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			ukzz[iky][ikx][ikz] *= tmpdt;
#else
			ukzz[iky][ikx][ikz] = sf_crmul(ukzz[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uzz[0][0], ukzz[0][0]);
	

    /* compute uxyp */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*kx*ky*2.0*(cosf(vp0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vp0*vp0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			ukxyp[iky][ikx][ikz] *= tmpdt; 
#else
			ukxyp[iky][ikx][ikz] = sf_crmul(ukxyp[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uxyp[0][0], ukxyp[0][0]);
	

    /* compute uxys */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*kx*ky*2.0*(cosf(vs0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vs0*vs0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			ukxys[iky][ikx][ikz] *= tmpdt;
#else
			ukxys[iky][ikx][ikz] = sf_crmul(ukxys[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uxys[0][0], ukxys[0][0]);
	

    /* compute uxzp */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*kx*kz*2.0*(cosf(vp0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vp0*vp0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			ukxzp[iky][ikx][ikz] *= tmpdt;
#else
			ukxzp[iky][ikx][ikz] = sf_crmul(ukxzp[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uxzp[0][0], ukxzp[0][0]);
	
    /* compute uxzs */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*kx*kz*2.0*(cosf(vs0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vs0*vs0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			ukxzs[iky][ikx][ikz] *= tmpdt; 
#else
			ukxzs[iky][ikx][ikz] = sf_crmul(ukxzs[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uxzs[0][0], ukxzs[0][0]);
	

    /* compute vxx */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*kx*kx*2.0*(cosf(vs0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vs0*vs0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			vkxx[iky][ikx][ikz] *= tmpdt; 
#else
			vkxx[iky][ikx][ikz] = sf_crmul(vkxx[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(vxx[0][0], vkxx[0][0]);
	
    /* compute vyy */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*ky*ky*2.0*(cosf(vp0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vp0*vp0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			vkyy[iky][ikx][ikz] *= tmpdt; 
#else
			vkyy[iky][ikx][ikz] = sf_crmul(vkyy[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(vyy[0][0], vkyy[0][0]);
	

    /* compute vzz */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*kz*kz*2.0*(cosf(vs0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vs0*vs0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			vkzz[iky][ikx][ikz] *= tmpdt; 
#else
			vkzz[iky][ikx][ikz] = sf_crmul(vkzz[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(vzz[0][0], vkzz[0][0]);
	

    /* compute vxyp */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*kx*ky*2.0*(cosf(vp0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vp0*vp0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			vkxyp[iky][ikx][ikz] *= tmpdt;
#else
			vkxyp[iky][ikx][ikz] = sf_crmul(vkxyp[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(vxyp[0][0], vkxyp[0][0]);
	

    /* compute vxys */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*kx*ky*2.0*(cosf(vs0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vs0*vs0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			vkxys[iky][ikx][ikz] *= tmpdt; 
#else
			vkxys[iky][ikx][ikz] = sf_crmul(vkxys[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(vxys[0][0], vkxys[0][0]);
	
    /* compute vyzp */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*ky*kz*2.0*(cosf(vp0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vp0*vp0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			vkyzp[iky][ikx][ikz] *= tmpdt;
#else
			vkyzp[iky][ikx][ikz] = sf_crmul(vkyzp[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(vyzp[0][0], vkyzp[0][0]);
	
    /* compute vyzs */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*ky*kz*2.0*(cosf(vs0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vs0*vs0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			vkyzs[iky][ikx][ikz] *= tmpdt;
#else
			vkyzs[iky][ikx][ikz] = sf_crmul(vkyzs[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(vyzs[0][0], vkyzs[0][0]);
	
    /* compute wxx */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*kx*kx*2.0*(cosf(vs0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vs0*vs0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			wkxx[iky][ikx][ikz] *= tmpdt;
#else
			wkxx[iky][ikx][ikz] = sf_crmul(wkxx[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(wxx[0][0], wkxx[0][0]);
	
    /* compute wyy */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*ky*ky*2.0*(cosf(vs0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vs0*vs0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			wkyy[iky][ikx][ikz] *= tmpdt;
#else
			wkyy[iky][ikx][ikz] = sf_crmul(wkyy[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(wyy[0][0], wkyy[0][0]);
	
    /* compute wzz*/
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*kz*kz*2.0*(cosf(vs0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vs0*vs0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			wkzz[iky][ikx][ikz] *= tmpdt; 
#else
			wkzz[iky][ikx][ikz] = sf_crmul(wkzz[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(wzz[0][0], wkzz[0][0]);

    /* compute wxzp */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*kx*kz*2.0*(cosf(vp0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vp0*vp0*dt*dt*kxyz2);
			//dehf(ikx,nkx,ax,factor)*dehf(ikz,nkz,az,factor);
#ifdef SF_HAS_COMPLEX_H
			wkxzp[iky][ikx][ikz] *= tmpdt; 
#else
			wkxzp[iky][ikx][ikz] = sf_crmul(wkxzp[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(wxzp[0][0], wkxzp[0][0]);
	

    /* compute wxzs */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*kx*kz*2.0*(cosf(vs0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vs0*vs0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			wkxzs[iky][ikx][ikz] *= tmpdt; 
#else
			wkxzs[iky][ikx][ikz] = sf_crmul(wkxzs[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(wxzs[0][0], wkxzs[0][0]);
	
    /* compute wyzp */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*ky*kz*2.0*(cosf(vp0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vp0*vp0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			wkyzp[iky][ikx][ikz] *= tmpdt; 
#else
			wkyzp[iky][ikx][ikz] = sf_crmul(wkyzp[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(wyzp[0][0], wkyzp[0][0]);
	
    /* compute wyzs */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (!kxyz2)
			kxyz2 +=0.000001;
			tmpdt = 1.0*ky*kz*2.0*(cosf(vs0*sqrtf(kx*kx+ky*ky+kz*kz)*dt)-1.0)/(vs0*vs0*dt*dt*kxyz2);
#ifdef SF_HAS_COMPLEX_H
			wkyzs[iky][ikx][ikz] *= tmpdt; 
#else
			wkyzs[iky][ikx][ikz] = sf_crmul(wkyzs[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(wyzs[0][0], wkyzs[0][0]);


#ifdef _OPENMP
#pragma omp parallel for private(ix, iy, iz)
#endif
    for (iy=0; iy < ny; iy++) {  
    for (ix=0; ix < nx; ix++) {  
        for (iz=0; iz < nz; iz++) {  
            upold[iy][ix][iz] = (dt*dt*vp[iy][ix][iz]*vp[iy][ix][iz]) *
			(uxx[iy][ix][iz]+vxyp[iy][ix][iz]+wxzp[iy][ix][iz]) + 2.0*upcur[iy][ix][iz]-upold[iy][ix][iz];
            vpold[iy][ix][iz] = (dt*dt*vp[iy][ix][iz]*vp[iy][ix][iz]) *
			(uxyp[iy][ix][iz]+vyy[iy][ix][iz]+wyzp[iy][ix][iz]) + 2.0*vpcur[iy][ix][iz]-vpold[iy][ix][iz];
            wpold[iy][ix][iz] = (dt*dt*vp[iy][ix][iz]*vp[iy][ix][iz]) *
			(uxzp[iy][ix][iz]+vyzp[iy][ix][iz]+wzz[iy][ix][iz]) + 2.0*wpcur[iy][ix][iz]-wpold[iy][ix][iz];

            usold[iy][ix][iz] = (dt*dt*vs[iy][ix][iz]*vs[iy][ix][iz]) *
			(uyy[iy][ix][iz]+uzz[iy][ix][iz]-vxys[iy][ix][iz]-wxzs[iy][ix][iz]) + 2.0*uscur[iy][ix][iz]-usold[iy][ix][iz];
            vsold[iy][ix][iz] = (dt*dt*vs[iy][ix][iz]*vs[iy][ix][iz]) *
			(vxx[iy][ix][iz]+vzz[iy][ix][iz]-uxys[iy][ix][iz]-wyzs[iy][ix][iz]) + 2.0*vscur[iy][ix][iz]-vsold[iy][ix][iz];
            wsold[iy][ix][iz] = (dt*dt*vs[iy][ix][iz]*vs[iy][ix][iz]) *
			(wxx[iy][ix][iz]+wyy[iy][ix][iz]-uxzs[iy][ix][iz]-vyzs[iy][ix][iz]) + 2.0*wscur[iy][ix][iz]-wsold[iy][ix][iz];
		}
	}
	}

}
