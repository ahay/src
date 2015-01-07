/* Modeling of pure acoustic wave in 3-D transversely isotropic meida using psuedospectral method */
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
#include "psstep3.h"
#include "fft3.h"
//#include "fft2d_JYAN.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static sf_complex ***ukxx, ***ukyy, ***ukzz, ***ukxy, ***ukxz, ***ukyz, ***uktmp;
static sf_complex ***ukxxxx, ***ukyyyy, ***ukzzzz, ***ukxxxy, ***ukxxxz, ***ukxyyy, ***ukyyyz, ***ukxzzz, ***ukyzzz, ***ukxxyy, ***ukxxzz, ***ukyyzz, ***ukxxyz, ***ukxyyz, ***ukxyzz;
static float ***uuxx, ***uuyy, ***uuzz, ***uuxy, ***uuxz, ***uuyz, ***curtmp;
static float ***uxxxx, ***uyyyy, ***uzzzz, ***uxxxy, ***uxxxz, ***uxyyy, ***uyyyz, ***uxzzz, ***uyzzz, ***uxxyy, ***uxxzz, ***uyyzz, ***uxxyz, ***uxyyz, ***uxyzz;
static float dkx, dky, dkz;
static int nkx, nky, nkz, opt;

void psstep3_init(int nz, int nx , int ny /*model size*/,
                  float dz, float dx, float dy /*model grid*/,
				  int opt1)
/*< initialization >*/
{
//	nkx = fft_nk(nx, opt1);
//	nky = fft_nk(ny, opt1);
//	nkz = fft_nk(nz, opt1);
	nkx = kiss_fft_next_fast_size(nx);
	nky = kiss_fft_next_fast_size(ny);
	nkz = kiss_fft_next_fast_size(nz);

    ukxx = sf_complexalloc3(nkz,nkx,nky);
    ukyy = sf_complexalloc3(nkz,nkx,nky);
    ukzz = sf_complexalloc3(nkz,nkx,nky);
    ukxy = sf_complexalloc3(nkz,nkx,nky);
    ukxz = sf_complexalloc3(nkz,nkx,nky);
    ukyz = sf_complexalloc3(nkz,nkx,nky);
    uktmp = sf_complexalloc3(nkz,nkx,nky);
    
	uuxx = sf_floatalloc3(nkz,nkx,nky);
    uuyy = sf_floatalloc3(nkz,nkx,nky);
    uuzz = sf_floatalloc3(nkz,nkx,nky);
    uuxy = sf_floatalloc3(nkz,nkx,nky);
    uuxz = sf_floatalloc3(nkz,nkx,nky);
    uuyz = sf_floatalloc3(nkz,nkx,nky);
    curtmp = sf_floatalloc3(nkz,nkx,nky);

    ukxxxx = sf_complexalloc3(nkz,nkx,nky);
    ukyyyy = sf_complexalloc3(nkz,nkx,nky);
    ukzzzz = sf_complexalloc3(nkz,nkx,nky);
    ukxxxy = sf_complexalloc3(nkz,nkx,nky);
    ukxxxz = sf_complexalloc3(nkz,nkx,nky);
    ukxyyy = sf_complexalloc3(nkz,nkx,nky);
    ukyyyz = sf_complexalloc3(nkz,nkx,nky);
    ukxzzz = sf_complexalloc3(nkz,nkx,nky);
    ukyzzz = sf_complexalloc3(nkz,nkx,nky);
    ukxxyy = sf_complexalloc3(nkz,nkx,nky);
    ukxxzz = sf_complexalloc3(nkz,nkx,nky);
    ukyyzz = sf_complexalloc3(nkz,nkx,nky);
    ukxxyz = sf_complexalloc3(nkz,nkx,nky);
    ukxyyz = sf_complexalloc3(nkz,nkx,nky);
    ukxyzz = sf_complexalloc3(nkz,nkx,nky);

    uxxxx = sf_floatalloc3(nkz,nkx,nky);
    uyyyy = sf_floatalloc3(nkz,nkx,nky);
    uzzzz = sf_floatalloc3(nkz,nkx,nky);
    uxxxy = sf_floatalloc3(nkz,nkx,nky);
    uxxxz = sf_floatalloc3(nkz,nkx,nky);
    uxyyy = sf_floatalloc3(nkz,nkx,nky);
    uyyyz = sf_floatalloc3(nkz,nkx,nky);
    uxzzz = sf_floatalloc3(nkz,nkx,nky);
    uyzzz = sf_floatalloc3(nkz,nkx,nky);
    uxxyy = sf_floatalloc3(nkz,nkx,nky);
    uxxzz = sf_floatalloc3(nkz,nkx,nky);
    uyyzz = sf_floatalloc3(nkz,nkx,nky);
    uxxyz = sf_floatalloc3(nkz,nkx,nky);
    uxyyz = sf_floatalloc3(nkz,nkx,nky);
    uxyzz = sf_floatalloc3(nkz,nkx,nky);

    dkx = 1./(nkx*dx);
    dky = 1./(nky*dy);
    dkz = 1./(nkz*dz);
	opt = opt1;
}

void psstep3_close(void)
/*< free memory allocation>*/
{

	free(**uktmp);
	free(*uktmp);
	free(uktmp);
	free(**curtmp);
	free(*curtmp);
	free(curtmp);

    free(**ukxx);     
    free(**ukyy);     
    free(**ukzz);     
    free(**ukxy);     
    free(**ukxz);     
    free(**ukyz);     
    free(**ukxxxx);    
    free(**ukyyyy);    
    free(**ukzzzz);    
    free(**ukxxxy);    
    free(**ukxxxz);    
    free(**ukxyyy);    
    free(**ukyyyz);    
    free(**ukxzzz);    
    free(**ukyzzz);    
    free(**ukxxyy);    
    free(**ukxxzz);    
    free(**ukyyzz);    
    free(**ukxxyz);    
    free(**ukxyyz);    
    free(**ukxyzz);    

    free(*ukxx);     
    free(*ukyy);     
    free(*ukzz);     
    free(*ukxy);     
    free(*ukxz);     
    free(*ukyz);     
    free(*ukxxxx);    
    free(*ukyyyy);    
    free(*ukzzzz);    
    free(*ukxxxy);    
    free(*ukxxxz);    
    free(*ukxyyy);    
    free(*ukyyyz);    
    free(*ukxzzz);    
    free(*ukyzzz);    
    free(*ukxxyy);    
    free(*ukxxzz);    
    free(*ukyyzz);    
    free(*ukxxyz);    
    free(*ukxyyz);    
    free(*ukxyzz);   

    free(ukxx);     
    free(ukyy);     
    free(ukzz);     
    free(ukxy);     
    free(ukxz);     
    free(ukyz);     
    free(ukxxxx);    
    free(ukyyyy);    
    free(ukzzzz);    
    free(ukxxxy);    
    free(ukxxxz);    
    free(ukxyyy);    
    free(ukyyyz);    
    free(ukxzzz);    
    free(ukyzzz);    
    free(ukxxyy);    
    free(ukxxzz);    
    free(ukyyzz);    
    free(ukxxyz);    
    free(ukxyyz);    
    free(ukxyzz);


    free(**uuxx);     
    free(**uuyy);     
    free(**uuzz);     
    free(**uuxy);     
    free(**uuxz);     
    free(**uuyz);     
    free(**uxxxx);    
    free(**uyyyy);    
    free(**uzzzz);    
    free(**uxxxy);    
    free(**uxxxz);    
    free(**uxyyy);    
    free(**uyyyz);    
    free(**uxzzz);    
    free(**uyzzz);    
    free(**uxxyy);    
    free(**uxxzz);    
    free(**uyyzz);    
    free(**uxxyz);    
    free(**uxyyz);    
    free(**uxyzz);    


    free(*uuxx);     
    free(*uuyy);     
    free(*uuzz);     
    free(*uuxy);     
    free(*uuxz);     
    free(*uuyz);     
    free(*uxxxx);    
    free(*uyyyy);    
    free(*uzzzz);    
    free(*uxxxy);    
    free(*uxxxz);    
    free(*uxyyy);    
    free(*uyyyz);    
    free(*uxzzz);    
    free(*uyzzz);    
    free(*uxxyy);    
    free(*uxxzz);    
    free(*uyyzz);    
    free(*uxxyz);    
    free(*uxyyz);    
    free(*uxyzz);   


    free(uuxx);     
    free(uuyy);     
    free(uuzz);     
    free(uuxy);     
    free(uuxz);     
    free(uuyz);     
    free(uxxxx);    
    free(uyyyy);    
    free(uzzzz);    
    free(uxxxy);    
    free(uxxxz);    
    free(uxyyy);    
    free(uyyyz);    
    free(uxzzz);    
    free(uyzzz);    
    free(uxxyy);    
    free(uxxzz);    
    free(uyyzz);    
    free(uxxyz);    
    free(uxyyz);    
    free(uxyzz);


}


void psstep3(float ***old /*previous step*/,
             float ***cur /*current step*/,
             int nz, int nx, int ny /*model size*/,
			 float dz, float dx, float dy /*model grid */,
             float v0 /*reference vel*/,
             float ***v /*reference vel*/,
             float ***sigma /*reference vel*/,
             float ***delta /*reference vel*/,
             float ***seta /*reference vel*/,
             float ***phi /*reference vel*/,
             float dt /*time step size*/)
/*< PS step>*/
{
    int ix, ikx, iy, iky, ikz, iz;
    float kx, ky, kz, tmpdt, pi=SF_PI;
	float kx0,ky0,kz0,kxyz2;
	int nkxx,nkyy, nkzz, nkxyz;

	kx0 =-0.5/dx;
	ky0 =-0.5/dy;
	kz0 =-0.5/dz;

    for (iky=0; iky < nky; iky++){ 
    for (ikx=0; ikx < nkx; ikx++){ 
        for (ikz=0; ikz < nkz; ikz++){ 
			ukxx[iky][ikx][ikz] = sf_cmplx(0., 0.);
			ukyy[iky][ikx][ikz] = sf_cmplx(0., 0.);
			ukzz[iky][ikx][ikz] = sf_cmplx(0., 0.);
			ukxy[iky][ikx][ikz] = sf_cmplx(0., 0.);
			ukxz[iky][ikx][ikz] = sf_cmplx(0., 0.);
			ukyz[iky][ikx][ikz] = sf_cmplx(0., 0.);
			
			ukxxxx[iky][ikx][ikz] = sf_cmplx(0., 0.);
			ukyyyy[iky][ikx][ikz] = sf_cmplx(0., 0.);
			ukzzzz[iky][ikx][ikz] = sf_cmplx(0., 0.);
			ukxxxy[iky][ikx][ikz] = sf_cmplx(0., 0.);
			ukxxxz[iky][ikx][ikz] = sf_cmplx(0., 0.);
			ukxyyy[iky][ikx][ikz] = sf_cmplx(0., 0.);
			ukyyyz[iky][ikx][ikz] = sf_cmplx(0., 0.);
			ukxzzz[iky][ikx][ikz] = sf_cmplx(0., 0.);
			ukyzzz[iky][ikx][ikz] = sf_cmplx(0., 0.);
			ukxxyy[iky][ikx][ikz] = sf_cmplx(0., 0.);
			ukxxzz[iky][ikx][ikz] = sf_cmplx(0., 0.);
			ukyyzz[iky][ikx][ikz] = sf_cmplx(0., 0.);
			ukxxyz[iky][ikx][ikz] = sf_cmplx(0., 0.);
			ukxyyz[iky][ikx][ikz] = sf_cmplx(0., 0.);
			ukxyzz[iky][ikx][ikz] = sf_cmplx(0., 0.);
			
			uktmp[iky][ikx][ikz] = sf_cmplx(0., 0.);
			
			uuxx[iky][ikx][ikz] = 0.;
			uuyy[iky][ikx][ikz] = 0.;
			uuzz[iky][ikx][ikz] = 0.;
			uuxy[iky][ikx][ikz] = 0.;
			uuxz[iky][ikx][ikz] = 0.;
			uuyz[iky][ikx][ikz] = 0.;
			
			uxxxx[iky][ikx][ikz] = 0.;
			uyyyy[iky][ikx][ikz] = 0.;
			uzzzz[iky][ikx][ikz] = 0.;
			uxxxy[iky][ikx][ikz] = 0.;
			uxxxz[iky][ikx][ikz] = 0.;
			uxyyy[iky][ikx][ikz] = 0.;
			uyyyz[iky][ikx][ikz] = 0.;
			uxzzz[iky][ikx][ikz] = 0.;
			uyzzz[iky][ikx][ikz] = 0.;
			uxxyy[iky][ikx][ikz] = 0.;
			uxxzz[iky][ikx][ikz] = 0.;
			uyyzz[iky][ikx][ikz] = 0.;
			uxxyz[iky][ikx][ikz] = 0.;
			uxyyz[iky][ikx][ikz] = 0.;
			uxyzz[iky][ikx][ikz] = 0.;

			curtmp[iky][ikx][ikz] = 0.;
         }  
	}
	}
	
    for (iy=0; iy < ny; iy++){ 
    for (ix=0; ix < nx; ix++){ 
        for (iz=0; iz < nz; iz++){ 
		    curtmp[iy][ix][iz] = cur[iy][ix][iz];	
         }  
	}
	}

	/* Just do one FFT and severn IFFT in one time step */

	/* Transform the input data into wavenumber space domain */
	fft3(curtmp[0][0], uktmp[0][0]);
/*	
fprintf(stderr, "Begin FFT \n");
	for (iky=0; iky<nky;iky++) {
	for (ikx=0; ikx<nkx;ikx++) {
		for (ikz=0;ikz < nkz;ikz++){
			fprintf(stderr, "%f ", uktmp[iky][ikx][ikz].r);
			if (ikz==nkz-1)
				fprintf(stderr, "\n");
		}
	}
	}
fprintf(stderr, "End FFT\n");
*/

	for (iky=0; iky < nky; iky++) {
	for (ikx=0; ikx < nkx; ikx++) {
		for (ikz=0; ikz < nkz; ikz++) {
			ukxx[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukyy[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukzz[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukxy[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukxz[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukyz[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			
			ukxxxx[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukyyyy[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukzzzz[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukxxxy[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukxxxz[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukxyyy[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukyyyz[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukxzzz[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukyzzz[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukxxyy[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukxxzz[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukyyzz[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukxxyz[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukxyyz[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
			ukxyzz[iky][ikx][ikz] = uktmp[iky][ikx][ikz];
		}
	}
	}
	

    /* compute ukxx(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = kx*kx;
#ifdef SF_HAS_COMPLEX_H
			ukxx[iky][ikx][ikz] *= tmpdt; 
#else
			ukxx[iky][ikx][ikz] = sf_crmul(ukxx[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uuxx[0][0], ukxx[0][0]);
	
    /* compute ukyy(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = ky*ky;
#ifdef SF_HAS_COMPLEX_H
			ukyy[iky][ikx][ikz] *= tmpdt; 
#else
			ukyy[iky][ikx][ikz] = sf_crmul(ukyy[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uuyy[0][0], ukyy[0][0]);


    /* compute ukzz(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = kz*kz;
#ifdef SF_HAS_COMPLEX_H
			ukzz[iky][ikx][ikz] *= tmpdt; 
#else
			ukzz[iky][ikx][ikz] = sf_crmul(ukzz[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uuzz[0][0], ukzz[0][0]);


    /* compute ukxy(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = kx*ky;
#ifdef SF_HAS_COMPLEX_H
			ukxy[iky][ikx][ikz] *= tmpdt; 
#else
			ukxy[iky][ikx][ikz] = sf_crmul(ukxy[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uuxy[0][0], ukxy[0][0]);


    /* compute ukxz(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = kx*kz;
#ifdef SF_HAS_COMPLEX_H
			ukxz[iky][ikx][ikz] *= tmpdt; 
#else
			ukxz[iky][ikx][ikz] = sf_crmul(ukxz[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uuxz[0][0], ukxz[0][0]);


    /* compute ukyz(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = ky*kz;
#ifdef SF_HAS_COMPLEX_H
			ukyz[iky][ikx][ikz] *= tmpdt;  
#else
			ukyz[iky][ikx][ikz] = sf_crmul(ukyz[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uuyz[0][0], ukyz[0][0]);


    /* compute ukxxxx(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = kx*kx*kx*kx/kxyz2;
#ifdef SF_HAS_COMPLEX_H
			ukxxxx[iky][ikx][ikz] *= tmpdt; 
#else
			ukxxxx[iky][ikx][ikz] = sf_crmul(ukxxxx[iky][ikx][ikz],tmpdt);
#endif
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uxxxx[0][0], ukxxxx[0][0]);


    /* compute ukyyyy(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = ky*ky*ky*ky/kxyz2;

			ukyyyy[iky][ikx][ikz] = sf_crmul(ukyyyy[iky][ikx][ikz],tmpdt);
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uyyyy[0][0], ukyyyy[0][0]);


    /* compute ukzzzz(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = kz*kz*kz*kz/kxyz2;
			ukzzzz[iky][ikx][ikz] = sf_crmul(ukzzzz[iky][ikx][ikz],tmpdt);
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uzzzz[0][0], ukzzzz[0][0]);


    /* compute ukxxxy(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = kx*kx*kx*ky/kxyz2;
			ukxxxy[iky][ikx][ikz] = sf_crmul(ukxxxy[iky][ikx][ikz],tmpdt);
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uxxxy[0][0], ukxxxy[0][0]);


    /* compute ukxxxz(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = kx*kx*kx*kz/kxyz2;
			ukxxxz[iky][ikx][ikz] = sf_crmul(ukxxxz[iky][ikx][ikz],tmpdt);
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uxxxz[0][0], ukxxxz[0][0]);


    /* compute ukxyyy(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = kx*ky*ky*ky/kxyz2;
			ukxyyy[iky][ikx][ikz] = sf_crmul(ukxyyy[iky][ikx][ikz],tmpdt);
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uxyyy[0][0], ukxyyy[0][0]);


    /* compute ukyyyz(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = ky*ky*ky*kz/kxyz2;
			ukyyyz[iky][ikx][ikz] = sf_crmul(ukyyyz[iky][ikx][ikz],tmpdt);
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uyyyz[0][0], ukyyyz[0][0]);


    /* compute ukxzzz(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = kx*kz*kz*kz/kxyz2;
			ukxzzz[iky][ikx][ikz] = sf_crmul(ukxzzz[iky][ikx][ikz],tmpdt);
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uxzzz[0][0], ukxzzz[0][0]);


    /* compute ukyzzz(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = ky*kz*kz*kz/kxyz2;
			ukyzzz[iky][ikx][ikz] = sf_crmul(ukyzzz[iky][ikx][ikz],tmpdt);
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uyzzz[0][0], ukyzzz[0][0]);


    /* compute ukxxyy(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = kx*kx*ky*ky/kxyz2;
			ukxxyy[iky][ikx][ikz] = sf_crmul(ukxxyy[iky][ikx][ikz],tmpdt);
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uxxyy[0][0], ukxxyy[0][0]);


    /* compute ukxxzz(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = kx*kx*kz*kz/kxyz2;
			ukxxzz[iky][ikx][ikz] = sf_crmul(ukxxzz[iky][ikx][ikz],tmpdt);
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uxxzz[0][0], ukxxzz[0][0]);


    /* compute ukyyzz(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = ky*ky*kz*kz/kxyz2;
			ukyyzz[iky][ikx][ikz] = sf_crmul(ukyyzz[iky][ikx][ikz],tmpdt);
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uyyzz[0][0], ukyyzz[0][0]);


    /* compute ukxxyz(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = kx*kx*ky*kz/kxyz2;
			ukxxyz[iky][ikx][ikz] = sf_crmul(ukxxyz[iky][ikx][ikz],tmpdt);
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uxxyz[0][0], ukxxyz[0][0]);


    /* compute ukxyyz(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = kx*ky*ky*kz/kxyz2;
			ukxyyz[iky][ikx][ikz] = sf_crmul(ukxyyz[iky][ikx][ikz],tmpdt);
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uxyyz[0][0], ukxyyz[0][0]);


    /* compute ukxyzz(ky,kx,kz) */
	for (iky=0; iky < nky; iky++) {
		ky = (ky0+iky*dky)*2.0*pi;
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxyz2 = kx*kx + ky*ky + kz*kz;
			if (kxyz2==0) {
			kxyz2 +=0.000001;
			}
			tmpdt = kx*ky*kz*kz/kxyz2;
			ukxyzz[iky][ikx][ikz] = sf_crmul(ukxyzz[iky][ikx][ikz],tmpdt);
		}
	}
	}
	/* Inverse FFT*/
	ifft3(uxyzz[0][0], ukxyzz[0][0]);


#ifdef _OPENMP
#pragma omp parallel for private(ix, iy, iz)
#endif
    for (iy=0; iy < ny; iy++) {  
    for (ix=0; ix < nx; ix++) {  
        for (iz=0; iz < nz; iz++) {  
            old[iy][ix][iz]  = (-1.0*dt*dt*v[iy][ix][iz]*v[iy][ix][iz])* (
			   	( 1.0+2.0*sigma[iy][ix][iz] + (2.0*delta[iy][ix][iz]-4.0*sigma[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*cosf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz]))*uuxx[iy][ix][iz]
			+	( 1.0+2.0*sigma[iy][ix][iz] + (2.0*delta[iy][ix][iz]-4.0*sigma[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz]))*uuyy[iy][ix][iz]
			+	( 1.0+2.0*sigma[iy][ix][iz] + (2.0*delta[iy][ix][iz]-4.0*sigma[iy][ix][iz])*cosf(seta[iy][ix][iz])*cosf(seta[iy][ix][iz]))*uuzz[iy][ix][iz]
			+	( (2.0*delta[iy][ix][iz]-4.0*sigma[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(2.0*phi[iy][ix][iz]))*uuxy[iy][ix][iz]
			+	( (2.0*delta[iy][ix][iz]-4.0*sigma[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*cosf(phi[iy][ix][iz]))*uuxz[iy][ix][iz]
			+	( (2.0*delta[iy][ix][iz]-4.0*sigma[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(phi[iy][ix][iz]))*uuyz[iy][ix][iz]
			+	( 2.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*cosf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz]))*uxxxx[iy][ix][iz]
			+	( 2.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz]))*uyyyy[iy][ix][iz]
			+	( 2.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*cosf(seta[iy][ix][iz])*cosf(seta[iy][ix][iz])*cosf(seta[iy][ix][iz])*cosf(seta[iy][ix][iz])*cosf(phi[iy][ix][iz]))*uzzzz[iy][ix][iz]
			+	( 4.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(2.0*phi[iy][ix][iz])*cosf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz]))*uxxxy[iy][ix][iz]
			+	( 4.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*cosf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz]))*uxxxz[iy][ix][iz]
			+	( 4.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(2.0*phi[iy][ix][iz])*sinf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz]))*uxyyy[iy][ix][iz]
			+	( 4.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz]))*uyyyz[iy][ix][iz]
			+	( 4.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*cosf(seta[iy][ix][iz])*cosf(seta[iy][ix][iz])*cosf(phi[iy][ix][iz]))*uxzzz[iy][ix][iz]
			+	( 4.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*cosf(seta[iy][ix][iz])*cosf(seta[iy][ix][iz])*sinf(phi[iy][ix][iz]))*uyzzz[iy][ix][iz]
			+	( 3.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(2.0*phi[iy][ix][iz])*sinf(2.0*phi[iy][ix][iz]))*uxxyy[iy][ix][iz]
			+	( 3.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*cosf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz]))*uxxzz[iy][ix][iz]
			+	( 3.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz]))*uyyzz[iy][ix][iz]
			+	( 12.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz])*cosf(phi[iy][ix][iz]))*uxxyz[iy][ix][iz]
			+	( 12.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*sinf(seta[iy][ix][iz])*cosf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz])*sinf(phi[iy][ix][iz]))*uxyyz[iy][ix][iz]
			+	( 3.0*(sigma[iy][ix][iz]-delta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(2.0*seta[iy][ix][iz])*sinf(2.0*phi[iy][ix][iz]))*uxyzz[iy][ix][iz]
			) + 2.0*cur[iy][ix][iz]-old[iy][ix][iz];

		}
	}  
	}


}

