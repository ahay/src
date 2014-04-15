/* multiply in 3D k-domain (kx, ky, kz)  */
/*
  Copyright (C) 2012 Tongji University (Jiubing Cheng) 
 
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
#include "_cjb.h"

#include <fftw3.h>
/* head files aumatically produced from *.c */

void spec3dmultiply(float ***d, float ***f, int nx, int ny, int nz, 
                    int *ijkx, int *ijky, int *ijkz, int iflag) 
/*< spec3dmultiply: multiply in 3D k-domain (kx, ky, kz) >*/
{
     int ix, iy, iz, ixf, iyf, izf, i;

#ifdef SF_HAS_FFTW  // using FFTW in Madagascar

       sf_warning("============= using SF_HAS_FFTW ====");

       sf_complex *xin, *xout;

       fftwf_plan xf;
       fftwf_plan yr;

       int m = nx*ny*nz;

       xin=sf_complexalloc(m);
       xout=sf_complexalloc(m);

       xf=fftwf_plan_dft_3d(ny,nx,nz, (fftwf_complex *) xin, (fftwf_complex *) xout,
                            FFTW_FORWARD,FFTW_ESTIMATE);

       yr=fftwf_plan_dft_3d(ny,nx,nz, (fftwf_complex *) xin, (fftwf_complex *) xout,
                            FFTW_BACKWARD,FFTW_ESTIMATE);

       /* FFT: from (x,y,z) to (kx, ky, kz) domain */
       i=0;
       for(iy=0;iy<ny;iy++)
       for(ix=0;ix<nx;ix++)
       for(iz=0;iz<nz;iz++){
            xin[i]=sf_cmplx(d[iy][ix][iz], 0.0);
            i++;
       }

       fftwf_execute(xf);

       int ii, jj, kk;

       i=0;
       for(iy=0;iy<ny;iy++){
          iyf=ijky[iy];
          ii = iy*nx*nz;
          for(ix=0;ix<nx;ix++){
            ixf=ijkx[ix];
            jj = ix*nz;
            for(iz=0;iz<nz;iz++){
              izf=ijkz[iz];
              kk = ii + jj + iz;
              //if(iy!=ny/2&&ix!=nx/2&&iz!=nz/2){
                 if(iflag==1)
                  xin[kk] = xout[kk]*sf_cmplx(f[iyf][ixf][izf], 0.0);    
                 else
                  xin[kk] = xout[kk]*sf_cmplx(0.0, f[iyf][ixf][izf]);    
              //}else
              //   xin[kk] = xout[kk];    
/*
                  xin[i] = xout[i];    
                  i++;
*/
            }// iz loop
          }// ix loop
       }//iy loop

       fftwf_execute(yr);

       i=0;
       for(iy=0;iy<ny;iy++)
       for(ix=0;ix<nx;ix++)
       for(iz=0;iz<nz;iz++){
            d[iy][ix][iz] = creal(xout[i])/m;
            i++;
       }

       fftwf_destroy_plan(xf);
       fftwf_destroy_plan(yr);

       free(xin);
       free(xout);

#endif
}
