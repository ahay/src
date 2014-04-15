/* separating 3D wave-modes based on low-rank decomposition using OpenMP */
/*
  Copyright (C) 2012 Tongji University (Jiubing Cheng) 
  and The University of Texas at Austin (Sergey Fomel)
 
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
#ifdef _OPENMP
#include <omp.h>
#endif
#include "_cjb.h"

#ifdef SF_HAS_FFTW
#include <fftw3.h>
#endif

/*****************************************************************************************/
void seplowrank3domp(float *ldata,float *rdata,float *fmid, float *p, int *ijkx, int *ijky, int *ijkz,
                      int nx, int ny, int nz, int m, int n, int m2, int n2, int iflag)
/*< seplowrank3domp: 3D wave-mode separation based on low-rank decomposition using OpenMP >*/
{

       //sf_warning("m2= %d n2=%d",m2,n2);

#ifdef SF_HAS_FFTW  // using FFTW in Madagascar
 
       //sf_warning("============= using SF_HAS_FFTW ====");

       sf_complex *xx, *xin, *xout;

       fftwf_plan xp;
       fftwf_plan xpi;

#endif
       float *wp;

        wp = sf_floatalloc(m*n2);

#ifdef SF_HAS_FFT
       xin=sf_complexalloc(m);
       xout=sf_complexalloc(n);
       xx=sf_complexalloc(n);

       xp=fftwf_plan_dft_3d(ny,nx,nz, (fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_FORWARD,FFTW_ESTIMATE);

       xpi=fftwf_plan_dft_3d(ny,nx,nz,(fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_BACKWARD,FFTW_ESTIMATE);

       // FFT: from (x,z) to (kx, kz) domain 

       if(iflag==1)
           for(i=0;i<m;i++) xin[i]=sf_cmplx(p[i], 0.);
       else 
           for(i=0;i<m;i++) xin[i]=sf_cmplx(0.0, p[i]);

       fftwf_execute(xp);
           
       for(i=0;i<n;i++) xx[i] = xout[i];

       // n2 IFFT from (kx, kz) to (x, z) domain
	   int jn2, ikx, iky, ikz, nxz, jn2n, iynxz, ii, io1, ixnz, iii, io2;
       nxz=nx*nz;

       for(jn2=0;jn2<n2;jn2++)
       {
           jn2n=jn2*n;
#ifdef _OPENMP
#pragma omp parallel for private(iky,ikx,ikz,iynxz,i,ii,iii,io1,io2) \
	   schedule(dynamic) \
	   shared(rdata, xx, xin, ijky, ijkx, ijkz, jn2n, ny, nx, nz, nxz)
#endif
           for(iky=0;iky<ny;iky++)
           {
              iynxz=ijky[iky]*nxz;
              ii=jn2n+iynxz;
			  io1=iky*nxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 ixnz=ijkx[ikx]*nz;
                 iii=ii+ixnz;
				 io2=io1+ikx*nz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
				   
                   i=io2+ikz;
                   xin[i]=rdata[iii+ijkz[ikz]]*xx[i];
                 } // ikz loop
             } // ikx loop
           }// iky loop

           // (kx,kz) to (x, z) domain
           fftwf_execute(xpi);

           for(im=0;im<m;im++)
             wp[jn2*m+im] = creal(xout[im])/n;
      }

     fftwf_destroy_plan(xp);
     fftwf_destroy_plan(xpi);

     free(xx);
     free(xin);
     free(xout);

   // Matrix multiplication in space-domain 
   for(im=0;im<m;im++)
   {
         sum1=0.0;
         for(im2=0;im2<m2;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2;jn2++)
              sum2 += fmid[im2*n2+jn2]*wp[jn2*m+im];

           sum1 += ldata[im*m2+im2]*sum2;
        }//im2 loop
        p[im] = sum1;
  } 

#endif

  free(wp);
}
