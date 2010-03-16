/* 1-D finite-difference wave extrapolation */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    int nx, nkx, nkz,  ix, ikx, ikz, nz, iz;
    float dx, dkx, kx, dz, dkz, kz;
    float **data;
    sf_file out, in;
     

    sf_init(argc,argv);
    out = sf_output("out");
    in = sf_input("in");   /* velocity */

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histint(in,"n2",&nz)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&dz)) sf_error("No d2= in input");

    nkx = nx;
    nkz = nz;
    sf_warning("%d",kiss_fft_next_fast_size(nx-1));
    dkx = 1./(2.0*kiss_fft_next_fast_size(nx-1)*dx);
    dkz = 1./(2.0*kiss_fft_next_fast_size(nz-1)*dz);

    sf_putint(out,"n1",nx);
    sf_putfloat(out,"d1",dkx);
    sf_putint(out,"n2",nz);
    sf_putfloat(out,"d2",dkz);
    sf_putfloat(out,"o1",0.);
    sf_putfloat(out,"o2",0.);




    data = sf_floatalloc2(nx,nz);


    sf_floatread(data[0],nx*nz,in);

/* compute u(kx,kz) */
    sf_cosft_init(nx);
    for (iz=0; iz < nz; iz++){
        /* Fourier transform x to kx */
        sf_cosft_frw(data[iz],0,1);
    }
    sf_cosft_close();

    sf_cosft_init(nz);
    for (ikx=0; ikx < nkx; ikx++){
    /* Fourier transform z to kz */
    sf_cosft_frw(data[0],ikx,nx);
    }
    sf_cosft_close();
//
/* Inverse FFT*/
 //   sf_cosft_init(nz);
  //  for (ikx=0; ikx < nkx; ikx++){
    /* Inverse Fourier transform kz to z */
  //  sf_cosft_inv(data[0],ikx,nx);
//    }
//    sf_cosft_close();
//    sf_cosft_init(nx);
//    for (iz=0; iz < nz; iz++){
    /* Inverse Fourier transform kx to x */
//        sf_cosft_inv(data[iz],0,1);
 //   }
//    sf_cosft_close();

    sf_floatwrite(data[0],nz*nx,out);

    exit(0); 
}           
           
