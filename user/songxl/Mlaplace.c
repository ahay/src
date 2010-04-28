/* 2-D Laplace operator */
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
#include <math.h>
#include <limits.h>
#include "abcpass.h"
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    int nx,  nz,  ix, iz;
    int nkx, nkz, ikx, ikz;
    float dx, dkx, kx, dz, dkz, kz, kx0, kz0, pi=SF_PI, o1, o2;
    float **old;
    kiss_fft_cpx **uk, *ctracex, *ctracez;
    kiss_fft_cfg cfgx, cfgxi, cfgz, cfgzi;

    sf_file out, img;
     

    sf_init(argc,argv);
    out = sf_output("out");
    img= sf_input("in");   /* input file*/

    if (SF_FLOAT != sf_gettype(img)) sf_error("Need float input");
    if (!sf_histint(img,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(img,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histint(img,"n2",&nz)) sf_error("No n2= in input");
    if (!sf_histfloat(img,"d2",&dz)) sf_error("No d2= in input");

    if (!sf_histfloat(img,"o1",&o1)) o1=0.0;
    if (!sf_histfloat(img,"o2",&o2)) o2=0.0;

    sf_putint(out,"n1",nx);
    sf_putfloat(out,"d1",dx);
    sf_putfloat(out,"o1",o1);

    sf_putint(out,"n2",nz);
    sf_putfloat(out,"d2",dz);
    sf_putfloat(out,"o2",o2);

    nkx = kiss_fft_next_fast_size(nx);
    nkz = kiss_fft_next_fast_size(nz);
    if (nkx != nx) sf_warning("nkx padded to %d",nkx);
    if (nkz != nz) sf_warning("nkz padded to %d",nkz);
//    dkx =  1./(2.0*kiss_fft_next_fast_size(nx-1)*dx);
//    dkz =  1./(2.0*kiss_fft_next_fast_size(nz-1)*dz);
    dkx = 1./(nkx*dx);
    kx0 = -0.5/dx;
    dkz = 1./(nkz*dz);
    kz0 = -0.5/dz;

    cfgx = kiss_fft_alloc(nkx,0,NULL,NULL);
    cfgxi = kiss_fft_alloc(nkx,1,NULL,NULL);
    cfgz = kiss_fft_alloc(nkz,0,NULL,NULL);
    cfgzi = kiss_fft_alloc(nkz,1,NULL,NULL);


    uk = (kiss_fft_cpx **) sf_complexalloc2(nkx,nkz);
    ctracex = (kiss_fft_cpx *) sf_complexalloc(nkx);
    ctracez = (kiss_fft_cpx *) sf_complexalloc(nkz);

           
           
           
    old     =  sf_floatalloc2(nx,nz);
    sf_floatread(old[0],nx*nz,img);

    for (iz=0; iz < nz; iz++){
        for (ix=0; ix < nx; ix++){
             uk[iz][ix].r = old[iz][ix];
             uk[iz][ix].i = 0.0;
           }
    }

/* compute  u(kx,kz) */
    for (iz=0; iz < nz; iz++){
         /* Fourier transform x to kx */
            for (ix=1; ix < nx; ix+=2){
                uk[iz][ix] = sf_cneg(uk[iz][ix]);
                }
            kiss_fft_stride(cfgx,uk[iz],ctracex,1);
            for (ikx=0; ikx<nkx; ikx++) uk[iz][ikx] = ctracex[ikx];
         }
     for (ikx=0; ikx < nkx; ikx++){
         /* Fourier transform z to kz */
            for (ikz=1; ikz<nkz; ikz+=2){
                uk[ikz][ikx] = sf_cneg(uk[ikz][ikx]);
                }
            kiss_fft_stride(cfgz,uk[0]+ikx,ctracez,nkx);
            for (ikz=0; ikz<nkz; ikz++) uk[ikz][ikx] = ctracez[ikz];
           }

/*    #ifdef _OPENMP
    #pragma omp parallel for private(ik,ix,x,k,tmp,tmpex,tmpdt) 
    #endif
*/        
          
         for (ikz=0; ikz < nkz; ikz++) {
              kz = ikz* dkz*2.0*pi;
           
              for (ikx=0; ikx < nkx; ikx++) {
                  kx = ikx * dkx*2.0*pi;
                  uk[ikz][ikx] = sf_crmul(uk[ikz][ikx],(kx*kx+kz*kz));
              }
           
         }   
/* Inverse FFT*/
         for (ikx=0; ikx < nkx; ikx++){
         /* Inverse Fourier transform kz to z */
             kiss_fft_stride(cfgzi,(kiss_fft_cpx *)uk[0]+ikx,(kiss_fft_cpx *)ctracez,nkx);
             for (ikz=0; ikz < nkz; ikz++) uk[ikz][ikx] = sf_crmul(ctracez[ikz],ikz%2?-1.0:1.0);
              }
             for (ikz=0; ikz < nkz; ikz++){
             /* Inverse Fourier transform kx to x */
                 kiss_fft_stride(cfgxi,(kiss_fft_cpx *)uk[ikz],(kiss_fft_cpx *)ctracex,1);
                 for (ikx=0; ikx < nkx; ikx++) uk[ikz][ikx] = sf_crmul(ctracex[ikx],ikx%2?-1.0:1.0);
             }

         for (iz=0; iz < nz; iz++){
             for (ix=0; ix < nx; ix++){
                  old[iz][ix] = sf_crealf(uk[iz][ix]);
                  old[iz][ix] /= (nkx*nkz);
                }
         }


    sf_floatwrite(old[0],nx*nz,out);
    free(*old);     
    free(old);     
    free(*uk);
    free(uk);
    free(ctracex);
    free(ctracez);
    
    exit(0); 
}           
           
