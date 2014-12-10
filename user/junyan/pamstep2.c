/* Modeling of pure acoustic wave in 2-D transversely isotropic media with psuedo-analytic method */
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
#include "pamstep2.h"
#include "fft2d_JYAN.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static kiss_fft_cpx **uk1, **uk2, **uk3, **uk4, **uk5, **uk6, **uk7;
static float **uu1, **uu2, **uu3, **uu4, **uu5, **uu6, **uu7;
static float dkx, dkz;
static int nkx, nkz, opt;

void pamstep2_init(int nz, int nx /*model size*/,
                  float dz, float dx /*model grid*/,
				  int opt1)
/*< initialization >*/
{
	nkx = fft_nk(nx, opt1);
	nkz = fft_nk(nz, opt1);
    uk1 = (kiss_fft_cpx **) sf_complexalloc2(nkz,nkx);
    uk2 = (kiss_fft_cpx **) sf_complexalloc2(nkz,nkx);
    uk3 = (kiss_fft_cpx **) sf_complexalloc2(nkz,nkx);
    uk4 = (kiss_fft_cpx **) sf_complexalloc2(nkz,nkx);
    uk5 = (kiss_fft_cpx **) sf_complexalloc2(nkz,nkx);
    uk6 = (kiss_fft_cpx **) sf_complexalloc2(nkz,nkx);
    uk7 = (kiss_fft_cpx **) sf_complexalloc2(nkz,nkx);

    uu1 = sf_floatalloc2(nz,nx);
    uu2 = sf_floatalloc2(nz,nx);
    uu3 = sf_floatalloc2(nz,nx);
    uu4 = sf_floatalloc2(nz,nx);
    uu5 = sf_floatalloc2(nz,nx);
    uu6 = sf_floatalloc2(nz,nx);
    uu7 = sf_floatalloc2(nz,nx);


    dkx = 1./(nkx*dx);
    dkz = 1./(nkz*dz);
	opt = opt1;
}

void pamstep2_close(void)
/*< free memory allocation>*/
{
    free(*uk1);     
    free(uk1);     
    free(*uk2);     
    free(uk2);     
    free(*uk3);     
    free(uk3);     
    free(*uk4);     
    free(uk4);     
    free(*uk5);     
    free(uk5);     
    free(*uk6);     
    free(uk6);     
    free(*uk7);     
    free(uk7);  

    free(*uu1);     
    free(uu1);     
    free(*uu2);     
    free(uu2);     
    free(*uu3);     
    free(uu3);     
    free(*uu4);     
    free(uu4);     
    free(*uu5);     
    free(uu5);     
    free(*uu6);     
    free(uu6);     
    free(*uu7);     
    free(uu7);    
}

void pamstep2(float **old /*previous step*/,
             float **cur /*current step*/,
             int nz, int nx /*model size*/,
			 float dz, float dx /*model grid */,
             float v0 /*reference vel*/,
             float **v /*reference vel*/,
             float **sigma /*reference vel*/,
             float **delta /*reference vel*/,
             float **seta /*reference vel*/,
             float dt /*time step size*/)
/*< PAM step>*/
{
    int ix, ikx, ikz, iz;
    float kx, kz, tmpdt, pi=SF_PI;
	float kx0,kz0,kxz2;

	kx0 =-0.5/dx;
	kz0 =-0.5/dz;

    for (ix=0; ix < nkx; ix++){ 
        for (iz=0; iz < nkz; iz++){ 
            uk1[ix][iz].r = 0.; 
            uk2[ix][iz].r = 0.; 
            uk3[ix][iz].r = 0.; 
            uk4[ix][iz].r = 0.; 
            uk5[ix][iz].r = 0.; 
            uk6[ix][iz].r = 0.; 
            uk7[ix][iz].r = 0.;
 
            uk1[ix][iz].i = 0.; 
            uk2[ix][iz].i = 0.; 
            uk3[ix][iz].i = 0.; 
            uk4[ix][iz].i = 0.; 
            uk5[ix][iz].i = 0.; 
            uk6[ix][iz].i = 0.; 
            uk7[ix][iz].i = 0.; 
         }  
	}
    for (ix=0; ix < nx; ix++){ 
        for (iz=0; iz < nz; iz++){ 
            uk1[ix][iz].r = cur[ix][iz]; 
            uk2[ix][iz].r = cur[ix][iz]; 
            uk3[ix][iz].r = cur[ix][iz]; 
            uk4[ix][iz].r = cur[ix][iz]; 
            uk5[ix][iz].r = cur[ix][iz]; 
            uk6[ix][iz].r = cur[ix][iz]; 
            uk7[ix][iz].r = cur[ix][iz]; 
         }  
	}

    /* compute u1(kx,kz) */
	fft2d_init(nz, nx, opt);
	fft2d_JYAN(false, nz, nx, uk1, NULL);
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxz2 = kx*kx + kz*kz;
			if (kxz2==0) {
			kxz2 +=0.000001;
			fprintf(stderr, "\n zero here \n");
			}
			tmpdt = -1.0*kx*kx*2.0*(cosf(v0*sqrtf(kx*kx+kz*kz)*dt)-1.0)/(v0*v0*dt*dt*kxz2);
			uk1[ikx][ikz] = sf_crmul(uk1[ikx][ikz],tmpdt);
		}
	}
	/* Inverse FFT*/
	fft2d_JYAN(true, nz, nx, uk1, uu1);
	fft2d_close();
	
    /* compute u2(kx,kz) */
	fft2d_init(nz, nx, opt);
	fft2d_JYAN(false, nz, nx, uk2, NULL);
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxz2 = kx*kx + kz*kz;
			if (!kxz2)
			kxz2 +=0.000001;
			tmpdt = -1.0*kz*kz*2.0*(cosf(v0*sqrtf(kx*kx+kz*kz)*dt)-1.0)/(v0*v0*dt*dt*kxz2);
			uk2[ikx][ikz] = sf_crmul(uk2[ikx][ikz],tmpdt);
		}
	}
	/* Inverse FFT*/
	fft2d_JYAN(true, nz, nx, uk2, uu2);
	fft2d_close();

    /* compute u3(kx,kz) */
	fft2d_init(nz, nx, opt);
    fft2d_JYAN(false, nz, nx, uk3, NULL);
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxz2 = kx*kx + kz*kz;
			if (!kxz2)
			kxz2 +=0.000001;
			tmpdt = -1.0*kx*kx*kx*kx/kxz2*2.0*(cosf(v0*sqrtf(kx*kx+kz*kz)*dt)-1.0)/(v0*v0*dt*dt*kxz2);
			uk3[ikx][ikz] = sf_crmul(uk3[ikx][ikz],tmpdt);
		}
	}
	/* Inverse FFT*/
    fft2d_JYAN(true, nz, nx, uk3, uu3);
	fft2d_close();

    /* compute u4(kx,kz) */
	fft2d_init(nz, nx, opt);
    fft2d_JYAN(false, nz, nx, uk4, NULL);
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxz2 = kx*kx + kz*kz;
			if (!kxz2)
			kxz2 +=0.000001;
			tmpdt = -1.0*kz*kz*kz*kz/kxz2*2.0*(cosf(v0*sqrtf(kx*kx+kz*kz)*dt)-1.0)/(v0*v0*dt*dt*kxz2);
			uk4[ikx][ikz] = sf_crmul(uk4[ikx][ikz],tmpdt);
		}
	}
	/* Inverse FFT*/
    fft2d_JYAN(true, nz, nx, uk4, uu4);
	fft2d_close();

    /* compute u5(kx,kz) */
	fft2d_init(nz, nx, opt);
    fft2d_JYAN(false, nz, nx, uk5, NULL);
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxz2 = kx*kx + kz*kz;
			if (!kxz2)
			kxz2 +=0.000001;
			tmpdt = -1.0*kx*kx*kx*kz/kxz2*2.0*(cosf(v0*sqrtf(kx*kx+kz*kz)*dt)-1.0)/(v0*v0*dt*dt*kxz2);
			uk5[ikx][ikz] = sf_crmul(uk5[ikx][ikz],tmpdt);
		}
	}
	/* Inverse FFT*/
    fft2d_JYAN(true, nz, nx, uk5, uu5);
	fft2d_close();


     /* compute u6(kx,kz) */
	fft2d_init(nz, nx, opt);
    fft2d_JYAN(false, nz, nx, uk6, NULL);
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxz2 = kx*kx + kz*kz;
			if (!kxz2)
			kxz2 +=0.000001;
			tmpdt = -1.0*kx*kz*kz*kz/kxz2*2.0*(cosf(v0*sqrtf(kx*kx+kz*kz)*dt)-1.0)/(v0*v0*dt*dt*kxz2);
			uk6[ikx][ikz] = sf_crmul(uk6[ikx][ikz],tmpdt);
		}
	}
	/* Inverse FFT*/
    fft2d_JYAN(true, nz, nx, uk6, uu6);
	fft2d_close();


    /* compute u7(kx,kz) */
	fft2d_init(nz, nx, opt);
    fft2d_JYAN(false, nz, nx, uk7, NULL);
	for (ikx=0; ikx < nkx; ikx++) {
		kx = (kx0+ikx*dkx)*2.0*pi;
		for (ikz=0; ikz < nkz; ikz++) {
			kz = (kz0+ikz*dkz)*2.0*pi;
			kxz2 = kx*kx + kz*kz;
			if (!kxz2)
			kxz2 +=0.000001;
			tmpdt = -1.0*kx*kx*kz*kz/kxz2*2.0*(cosf(v0*sqrtf(kx*kx+kz*kz)*dt)-1.0)/(v0*v0*dt*dt*kxz2);
			uk7[ikx][ikz] = sf_crmul(uk7[ikx][ikz],tmpdt);
		}
	}
	/* Inverse FFT*/
    fft2d_JYAN(true, nz, nx, uk7, uu7);
	fft2d_close();


#ifdef _OPENMP
#pragma omp parallel for private(ix, iz)
#endif
    for (ix=0; ix < nx; ix++) {  
        for (iz=0; iz < nz; iz++) {  
            old[ix][iz]  = (-1.0*dt*dt*v[ix][iz]*v[ix][iz])* ( uu1[ix][iz] + uu2[ix][iz]
					     + (2.0*sigma[ix][iz]*cosf(seta[ix][iz])*cosf(seta[ix][iz])*cosf(seta[ix][iz])*cosf(seta[ix][iz]) + 2.0*delta[ix][iz]*sinf(seta[ix][iz])*sinf(seta[ix][iz])*cosf(seta[ix][iz])*cosf(seta[ix][iz]))*uu3[ix][iz]
					     + (2.0*sigma[ix][iz]*sinf(seta[ix][iz])*sinf(seta[ix][iz])*sinf(seta[ix][iz])*sinf(seta[ix][iz]) + 2.0*delta[ix][iz]*sinf(seta[ix][iz])*sinf(seta[ix][iz])*cosf(seta[ix][iz])*cosf(seta[ix][iz]))*uu4[ix][iz]
					     + (-4.0*sigma[ix][iz]*sinf(2.0*seta[ix][iz])*cosf(seta[ix][iz])*cosf(seta[ix][iz]) + delta[ix][iz]*sinf(4.0*seta[ix][iz]))*uu5[ix][iz]
					     + (-4.0*sigma[ix][iz]*sinf(2.0*seta[ix][iz])*sinf(seta[ix][iz])*sinf(seta[ix][iz]) - delta[ix][iz]*sinf(4.0*seta[ix][iz]))*uu6[ix][iz]
					     + (3.0*sigma[ix][iz]*sinf(2.0*seta[ix][iz])*sinf(2.0*seta[ix][iz]) - delta[ix][iz]*sinf(2.0*seta[ix][iz])*sinf(2.0*seta[ix][iz]) + 2.0*delta[ix][iz]*cosf(2.0*seta[ix][iz])*cosf(2.0*seta[ix][iz]))*uu7[ix][iz] )			   + 2.0*cur[ix][iz]-old[ix][iz];
        }
	}  

}

