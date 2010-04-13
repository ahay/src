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
#include <math.h>
#include <limits.h>
#include "abcpass.h"
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    int nx,  nz,  iz;
    int nkx, nkz, ikx, ikz;
    float dx, dkx, kx, dz, dkz, kz, pi=SF_PI, o1, o2;
    float **old;
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

    nkx = nx;
    nkz = nz;
    dkx = 1./(2.0*kiss_fft_next_fast_size(nx-1)*dx);
    dkz = 1./(2.0*kiss_fft_next_fast_size(nz-1)*dz);





    old    =  sf_floatalloc2(nx,nz);

    
    sf_floatread(old[0],nx*nz,img);
/* compute u(kx,kz) */
         sf_cosft_init(nx);
         for (iz=0; iz < nz; iz++){
             /* Fourier transform x to kx */
             sf_cosft_frw(old[iz],0,1);
         }
         sf_cosft_close();

         sf_cosft_init(nz);
         for (ikx=0; ikx < nkx; ikx++){
             /* Fourier transform z to kz */
             sf_cosft_frw(old[0],ikx,nx);
         }
         sf_cosft_close();

/*    #ifdef _OPENMP
    #pragma omp parallel for private(ik,ix,x,k,tmp,tmpex,tmpdt) 
    #endif
*/

         for (ikz=0; ikz < nkz; ikz++) {
             kz = ikz* dkz*2.0*pi;

             for (ikx=0; ikx < nkx; ikx++) {
                 kx = ikx * dkx*2.0*pi;
                 old[ikz][ikx] = old[ikz][ikx]*(kx*kx+kz*kz);
             }

         }   
/* Inverse FFT*/
         sf_cosft_init(nz);
         for (ikx=0; ikx < nkx; ikx++){
             /* Inverse Fourier transform kz to z */
             sf_cosft_inv(old[0],ikx,nx);
         }
         sf_cosft_close();
         sf_cosft_init(nx);
         for (iz=0; iz < nz; iz++){
             /* Inverse Fourier transform kx to x */
              sf_cosft_inv(old[iz],0,1);
         }
         sf_cosft_close();

    sf_floatwrite(old[0],nx*nz,out);
    free(*old);     
    free(old);     
    exit(0); 
}           
           
