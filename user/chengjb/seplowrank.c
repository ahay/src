/* separating wave-modes based on low-rank decomposition */
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
#include "_cjb.h"

#include <fftw3.h>

/*****************************************************************************************/
void seplowrank(float *ldata,float *rdata,float *fmid, float *x, int *ijkx, int *ijkz,
                int nx,int nz,int m,int n,int m2,int n2, int iflag)
/*< seplowrank: separating wave-modes based on low-rank decomposition >*/
{
       int i, im, im2, jn2, ikx, ikz;
       float sum1, sum2, *wp;

       wp = sf_floatalloc(m*n2);

/*
 * Note:  m=nx*nz; n=nkx*nkz;
 *        
 *        x[nx*nz]:      1D array for input and output 2D wavefield
 *        ldata[m*m2]:   1D array for left matrix from low-rank decomposition
 *        rdata[n2*n]:   1D array for right matrix from low-rank decomposition
 *        fmid[m2*n2]:   1D array for mid matrix from low-rank decomposition
 */

#ifdef SF_HAS_FFTW  // using FFTW in Madagascar
 
       sf_warning("============= using SF_HAS_FFTW ====");

       sf_complex *xin, *xout;

       fftwf_plan xp;
       fftwf_plan xpi;

       xin=sf_complexalloc(m);
       xout=sf_complexalloc(n);

       xp=fftwf_plan_dft_2d(nx,nz, (fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_FORWARD,FFTW_ESTIMATE);

       xpi=fftwf_plan_dft_2d(nx,nz,(fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_BACKWARD,FFTW_ESTIMATE);

       /* FFT: from (x,z) to (kx, kz) domain */

       if(iflag==1)
          for(i=0;i<m;i++) xin[i]=sf_cmplx(x[i], 0.);
       else
          for(i=0;i<m;i++) xin[i]=sf_cmplx(0.0, x[i]);

       fftwf_execute(xp);
           
       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       /* n2 IFFT from (kx, kz) to (x, z) domain*/

/* for test
       i=0;
       for(ikx=0;ikx<nx;ikx++)
          for(ikz=0;ikz<nz;ikz++)
          {
                 xin[i]=xout[i];          
                 i++;
          }
       // (kx,kz) to (x, z) domain
       fftwf_execute(xpi);
       for(im=0;im<m;im++)
                x[im] = creal(xout[im])/n;
*/

       for(jn2=0;jn2<n2;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(ikx=0;ikx<nx;ikx++)
           {
              // Note: Spectrum of the operator is differently orderred as the spectrum after FFT
              int ixnz=ijkx[ikx]*nz;
              int ii=jn2n+ixnz;
              for(ikz=0;ikz<nz;ikz++)
              {
                 xin[i]=rdata[ii+ijkz[ikz]]*xout[i];          
                 i++;
              }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }
       fftwf_destroy_plan(xp);
       fftwf_destroy_plan(xpi);
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
         x[im] = sum1;
       } 


#else  // using FFTW in user's own computer
       sf_warning("============= using user installed FFTW ====");

       fftw_complex *xin, *xout;

       fftw_plan xp;
       fftw_plan xpi;

       xin=fftw_complexalloc(m);
       xout=fftw_complexalloc(n);

       xp=fftw_plan_dft_2d(nx,nz, (fftw_complex *) xin, (fftw_complex *) xout,
			    FFTW_FORWARD,FFTW_ESTIMATE);
       xpi=fftw_plan_dft_2d(nx,nz,(fftw_complex *) xin, (fftw_complex *) xout,
			    FFTW_BACKWARD,FFTW_ESTIMATE);

       /* FFT: from (x,z) to (kx, kz) domain */
       for(i=0;i<m;i++)
       {
          xin[i][0]=x[i];
          xin[i][1]=0.;
       }

       fftw_execute(xp);
           
       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       for(jn2=0;jn2<n2;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(ikx=0;ikx<nx;ikx++)
           {
              /* Note: Spectrum of the operator is differently orderred as the spectrum after FFT */ 
              int ixnz=ijkx[ikx]*nz;
              int ii=jn2n+ixnz;
              for(ikz=0;ikz<nz;ikz++)
              {
                 xin[i]=rdata[ii+ijkz[ikz]]*xout[i];          
                 i++;
              }
            }

            /* (kx,kz) to (x, z) domain */
            fftw_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = xout[im][0]/n;
       }
       fftw_destroy_plan(xp);
       fftw_destroy_plan(xpi);
       free(xin);
       free(xout);

       /* Matrix multiplication in space-domain */
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
         x[im] = sum1;
       } 
#endif
       free(wp);
}

/*****************************************************************************************/
void sep(float *w, float *x, int *ijkx, int *ijkz, int nx,int nz,int m,int n, int iflag)
/*< sep: separating wave-modes by filtering in wavenumber domain >*/
{
       int i, ikx, ikz, im;

#ifdef SF_HAS_FFTW  // using FFTW in Madagascar
 
       sf_warning("============= using SF_HAS_FFTW ====");

// test FFT
/*
       int ntest=30;
       sf_complex *xi, *xo;
       xi=sf_complexalloc(ntest);
       xo=sf_complexalloc(ntest);
       fftwf_plan xf;
       xf=fftwf_plan_dft_1d(ntest,(fftwf_complex *) xi, (fftwf_complex *) xo,
			    FFTW_FORWARD,FFTW_ESTIMATE);
       for(i=0;i<ntest;i++)
          xi[i] = sf_cmplx(0.0, 0.0);
       for(i=ntest/2-3;i<ntest/2+3;i++)
          xi[i] = sf_cmplx(sin(i*1.0), 0.0);
       fftwf_execute(xf);
       for(i=0;i<ntest;i++)
          sf_warning("xo[%d]=(%f, %f)",i,creal(xo[i]),cimag(xo[i]));
       free(xi);
       free(xo);
       fftwf_destroy_plan(xf);
       exit(0);
*/

       sf_complex *xin, *xout;

       fftwf_plan xp;
       fftwf_plan xpi;

       xin=sf_complexalloc(m);
       xout=sf_complexalloc(n);

       xp=fftwf_plan_dft_2d(nx,nz, (fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_FORWARD,FFTW_ESTIMATE);

       xpi=fftwf_plan_dft_2d(nx,nz,(fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_BACKWARD,FFTW_ESTIMATE);

       /* (x,z) to (kx, kz) domain */
       if(iflag==1)
          for(i=0;i<m;i++) xin[i]=sf_cmplx(x[i], 0.);
       else
          for(i=0;i<m;i++) xin[i]=sf_cmplx(0.0, x[i]);


       fftwf_execute(xp);
           
       /* Filtering in (kx, kz) domain */
       i=0;
       for(ikx=0;ikx<nx;ikx++)
       {
          /* Note: Spectrum of the operator is differently orderred as the spectrum after FFT */ 
          int ixnz=ijkx[ikx]*nz;
          for(ikz=0;ikz<nz;ikz++)
          {
             //sf_warning("xout[%d,%d]= (%f,%f)",ikx,ikz,creal(xout[i]),cimag(xout[i]));             
             xin[i]=w[ixnz+ijkz[ikz]]*xout[i];          
             i++;
          }
       }

       /* IFFT from (kx,kz) to (x, z) domain */
       fftwf_execute(xpi);

       for(im=0;im<m;im++)
             x[im] = creal(xout[im])/n;
       
       fftwf_destroy_plan(xp);
       fftwf_destroy_plan(xpi);
       free(xin);
       free(xout);

#else  // using FFTW in user's own computer
       sf_warning("============= using user installed FFTW ====");

       fftw_complex *xin, *xout;

       fftw_plan xp;
       fftw_plan xpi;

       xin=fftw_complexalloc(m);
       xout=fftw_complexalloc(n);

       xp=fftw_plan_dft_2d(nx,nz, (fftw_complex *) xin, (fftw_complex *) xout,
			    FFTW_FORWARD,FFTW_ESTIMATE);
       xpi=fftw_plan_dft_2d(nx,nz,(fftw_complex *) xin, (fftw_complex *) xout,
			    FFTW_BACKWARD,FFTW_ESTIMATE);

       /* (x,z) to (kx, kz) domain */
       for(i=0;i<m;i++)
       {
          xin[i][0]=x[i];
          xin[i][1]=0.;
       }

       fftw_execute(xp);
           
       /* Filtering in (kx, kz) domain */
       i=0;
       for(ikx=0;ikx<nx;ikx++)
       {
          /* Note: Spectrum of the operator is differently orderred as the spectrum after FFT */ 
          int ixnz=ijkx[ikx]*nz;
          for(ikz=0;ikz<nz;ikz++)
          {
             xin[i]=w[ixnz+ijkz[ikz]]*xout[i];          
             i++;
          }
       }

       /* IFFT from (kx,kz) to (x, z) domain */
       fftw_execute(xpi);

       for(im=0;im<m;im++)
           x[im] = xout[im][0]/n;
      
       fftw_destroy_plan(xp);
       fftw_destroy_plan(xpi);
       free(xin);
       free(xout);

#endif
}
