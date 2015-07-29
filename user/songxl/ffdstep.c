/* 2-D Fourier finite-difference wave extrapolation */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
#include "ffdstep.h"
static float **uk, dkx, dkz;
static int nkx, nkz;

void ffdstep_init(int nx, int nz /*model size*/,
                  float dx, float dz /*model grid*/)
/*< initialization >*/
{
    uk  =  sf_floatalloc2(nx,nz);
    nkx = nx;
    nkz = nz;
    dkx = 1./(2.0*kiss_fft_next_fast_size(nx-1)*dx);
    dkz = 1./(2.0*kiss_fft_next_fast_size(nz-1)*dz);
}

void ffdstep_close(void)
/*< free memory allocation>*/
{
    free(*uk);     
    free(uk);     
}

void ffdstep(float **old /*previous step*/,
             float **cur /*current step*/,
             float ***aa /*precomputed coefficients*/,
             int nx, int nz /*model size*/,
             float v0 /*reference vel*/,
             float dt /*time step size*/)
/*< FFD step>*/
{
    int ix, ikx, ikz, iz;
    float kx, kz, tmpdt, pi=SF_PI;
    for (iz=0; iz < nz; iz++){
        for (ix=0; ix < nx; ix++){ 
            uk[iz][ix] = cur[iz][ix]; 
            }
         }  

    /* compute u(kx,kz) */
    sf_cosft_init(nx);
    for (iz=0; iz < nz; iz++){
        /* Fourier transform x to kx */
        sf_cosft_frw(uk[iz],0,1);
    }
    sf_cosft_close();

    sf_cosft_init(nz);
    for (ikx=0; ikx < nkx; ikx++){
        /* Fourier transform z to kz */
        sf_cosft_frw(uk[0],ikx,nx);
    }
    sf_cosft_close();

    for (ikz=0; ikz < nkz; ikz++) {
        kz = ikz* dkz*2.0*pi;

        for (ikx=0; ikx < nkx; ikx++) {
            kx = ikx * dkx*2.0*pi;
            tmpdt = 2.0*(cosf(v0*sqrtf(kx*kx+kz*kz)*dt)-1.0)/(v0*v0);
            uk[ikz][ikx] = uk[ikz][ikx]*tmpdt;
        }

    }   
/* Inverse FFT*/
    sf_cosft_init(nz);
    for (ikx=0; ikx < nkx; ikx++){
        /* Inverse Fourier transform kz to z */
        sf_cosft_inv(uk[0],ikx,nx);
    }
    sf_cosft_close();
    sf_cosft_init(nx);
    for (iz=0; iz < nz; iz++){
    /* Inverse Fourier transform kx to x */
        sf_cosft_inv(uk[iz],0,1);
    }
    sf_cosft_close();

    for (iz=1; iz < nz-1; iz++) {  
        for (ix=1; ix < nx-1; ix++) {  
            old[iz][ix]  = uk[iz][ix]*aa[iz][ix][0]
                         + (uk[iz][ix-1]+uk[iz][ix+1])*aa[iz][ix][1]
                         + (uk[iz-1][ix]+uk[iz+1][ix])*aa[iz][ix][2]
                         + 2.0*cur[iz][ix]-old[iz][ix];
         }
    }  

}           



           
float dehf(int k /*current frequency*/,
          int kn /*highest frequency*/,
          float a /*suppress factor*/,
          float factor /*propotion*/)
/*< high frequency depressing>*/
{
    int kmax;
    float depress;
    kmax = (int) (kn*factor);
    if (k < kmax) {
       depress = 1.0;
       }
    else
	/* depress =cosf(((float)(k-kmax))/((float)(kn-kmax))*pi/2.0); */
	/* depress = exp(-a*(float)((k-kmax)*(k-kmax))/((float)((kn-kmax)*(kn-kmax)))); */
       depress = exp(-a*(float)((k-kmax)*(k-kmax)*(k-kmax))/((float)((kn-kmax)*(kn-kmax)*(kn-kmax))));
    return(depress);
}

void ffdstep_dehf(float **old /*previous step*/,
             float **cur /*current step*/,
             float ***aa /*precomputed coefficients*/,
             int nx, int nz /*model size*/,
             float v0 /*reference vel*/,
             float dt /*time step size*/,
             float ax, float az /*suppress factor*/,
             float factor /*propotion*/)
/*< FFD step>*/
{
    int ix, ikx, ikz, iz;
    float kx, kz, tmpdt, pi=SF_PI;
    for (iz=0; iz < nz; iz++){
        for (ix=0; ix < nx; ix++){ 
            uk[iz][ix] = cur[iz][ix]; 
            }
         }  

    /* compute u(kx,kz) */
    sf_cosft_init(nx);
    for (iz=0; iz < nz; iz++){
        /* Fourier transform x to kx */
        sf_cosft_frw(uk[iz],0,1);
    }
    sf_cosft_close();

    sf_cosft_init(nz);
    for (ikx=0; ikx < nkx; ikx++){
        /* Fourier transform z to kz */
        sf_cosft_frw(uk[0],ikx,nx);
    }
    sf_cosft_close();


    for (ikz=0; ikz < nkz; ikz++) {
        kz = ikz* dkz*2.0*pi;

        for (ikx=0; ikx < nkx; ikx++) {
            kx = ikx * dkx*2.0*pi;
            tmpdt = 2.0*(cosf(v0*sqrtf(kx*kx+kz*kz)*dt)-1.0)/(v0*v0)*dehf(ikx,nkx,ax,factor)*dehf(ikz,nkz,az,factor);
            uk[ikz][ikx] = uk[ikz][ikx]*tmpdt;
        }

    }   
/* Inverse FFT*/
    sf_cosft_init(nz);
    for (ikx=0; ikx < nkx; ikx++){
        /* Inverse Fourier transform kz to z */
        sf_cosft_inv(uk[0],ikx,nx);
    }
    sf_cosft_close();
    sf_cosft_init(nx);
    for (iz=0; iz < nz; iz++){
    /* Inverse Fourier transform kx to x */
        sf_cosft_inv(uk[iz],0,1);
    }
    sf_cosft_close();

    for (iz=1; iz < nz-1; iz++) {  
        for (ix=1; ix < nx-1; ix++) {  
            old[iz][ix]  = uk[iz][ix]*aa[iz][ix][0]
                         + (uk[iz][ix-1]+uk[iz][ix+1])*aa[iz][ix][1]
                         + (uk[iz-1][ix]+uk[iz+1][ix])*aa[iz][ix][2]
                         + 2.0*cur[iz][ix]-old[iz][ix];
         }
    }  

}           
