/* wave-modes vector decomposition based on low-rank decomposition */
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

#include <fftw3.h>

/*****************************************************************************************/
void decomplowrank2d(float *ldataxx,float *rdataxx,float *fmidxx,
                      float *ldataxz,float *rdataxz,float *fmidxz,
                      float *ldatazz,float *rdatazz,float *fmidzz,
                      float *px, float *pz, int *ijkx, int *ijkz,
                      int nx, int nz, int m, int n, int MM,
                      int m2xx, int n2xx, int m2xz, int n2xz, int m2zz, int n2zz)
/*< decomplowrank2d: vector decomposition based on low-rank decomposition >*/
{
       int   i, im, im2, jn2, ikx, ikz;
       float sum1, sum2, *wp;

#ifdef SF_HAS_FFTW  // using FFTW in Madagascar
 
       sf_warning("============= using SF_HAS_FFTW ====");

       sf_complex *pxx, *xin, *xout;
       sf_complex *pzz;

       fftwf_plan xp;
       fftwf_plan xpi;

       xin=sf_complexalloc(m);
       xout=sf_complexalloc(n);
       pxx=sf_complexalloc(n);
       pzz=sf_complexalloc(n);

       xp=fftwf_plan_dft_2d(nx,nz, (fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_FORWARD,FFTW_ESTIMATE);

       xpi=fftwf_plan_dft_2d(nx,nz,(fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_BACKWARD,FFTW_ESTIMATE);

       /* FFT: from (x,z) to (kx, kz) domain */
       for(i=0;i<m;i++){
           xin[i]=sf_cmplx(px[i], 0.);
           px[i] = 0.0;
       }
       fftwf_execute(xp);
       for(i=0;i<n;i++) pxx[i] = xout[i];

       for(i=0;i<m;i++){
           xin[i]=sf_cmplx(pz[i], 0.);
           pz[i] = 0.0;
       }
       fftwf_execute(xp);
       for(i=0;i<n;i++) pzz[i] = xout[i];

       ///////////////////////////////////////////////////////// P-wave's x-component
       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2xx);
       for(jn2=0;jn2<n2xx;jn2++)
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
                xin[i]=rdataxx[ii+ijkz[ikz]]*pxx[i];          
                i++;
              }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2xx;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2xx;jn2++)
              sum2 += fmidxx[im2*n2xx+jn2]*wp[jn2*m+im];

           sum1 += ldataxx[im*m2xx+im2]*sum2;
         }//im2 loop
         px[im] = sum1;
       } 
       free(wp);

       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2xz);
       for(jn2=0;jn2<n2xz;jn2++)
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
                xin[i]=rdataxz[ii+ijkz[ikz]]*pzz[i];          
                i++;
              }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2xz;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2xz;jn2++)
              sum2 += fmidxz[im2*n2xz+jn2]*wp[jn2*m+im];

           sum1 += ldataxz[im*m2xz+im2]*sum2;
         }//im2 loop
         px[im] += sum1;
       } 
       free(wp);
       ///////////////////////////////////////////////////////// P-wave's z-component
       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2zz);
       for(jn2=0;jn2<n2zz;jn2++)
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
                xin[i]=rdatazz[ii+ijkz[ikz]]*pzz[i];          
                i++;
              }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2zz;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2zz;jn2++)
              sum2 += fmidzz[im2*n2zz+jn2]*wp[jn2*m+im];

           sum1 += ldatazz[im*m2zz+im2]*sum2;
         }//im2 loop
         pz[im] = sum1;
       } 
       free(wp);

       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2xz);
       for(jn2=0;jn2<n2xz;jn2++)
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
                xin[i]=rdataxz[ii+ijkz[ikz]]*pxx[i];          
                i++;
              }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2xz;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2xz;jn2++)
              sum2 += fmidxz[im2*n2xz+jn2]*wp[jn2*m+im];

           sum1 += ldataxz[im*m2xz+im2]*sum2;
         }//im2 loop
         pz[im] += sum1;
       } 
       free(wp);

       fftwf_destroy_plan(xp);
       fftwf_destroy_plan(xpi);

       free(pxx);
       free(pzz);
       free(xin);
       free(xout);
#else  // using FFTW in user's own computer
       sf_warning("============= using user installed FFTW ====");
#endif
}

/*****************************************************************************************/
void decomplowrank2ds(float *ldataxx,float *rdataxx,float *fmidxx,
                      float *ldataxz,float *rdataxz,float *fmidxz,
                      float *ldatazz,float *rdatazz,float *fmidzz,
                      float *px, float *pz, int *ijkx, int *ijkz,
                      int nx, int nz, int m, int n, int MM,
                      int m2xx, int n2xx, int m2xz, int n2xz, int m2zz, int n2zz)
/*< decomplowrank2ds: SV-wave vector decomposition based on low-rank decomposition >*/
{
       int   i, im, im2, jn2, ikx, ikz;
       float sum1, sum2, *wp;

#ifdef SF_HAS_FFTW  // using FFTW in Madagascar
 
       sf_warning("============= using SF_HAS_FFTW ====");

       sf_complex *pxx, *xin, *xout;
       sf_complex *pzz;

       fftwf_plan xp;
       fftwf_plan xpi;

       xin=sf_complexalloc(m);
       xout=sf_complexalloc(n);
       pxx=sf_complexalloc(n);
       pzz=sf_complexalloc(n);

       xp=fftwf_plan_dft_2d(nx,nz, (fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_FORWARD,FFTW_ESTIMATE);

       xpi=fftwf_plan_dft_2d(nx,nz,(fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_BACKWARD,FFTW_ESTIMATE);

       /* FFT: from (x,z) to (kx, kz) domain */
       for(i=0;i<m;i++){
           xin[i]=sf_cmplx(px[i], 0.);
           px[i] = 0.0;
       }
       fftwf_execute(xp);
       for(i=0;i<n;i++) pxx[i] = xout[i];

       for(i=0;i<m;i++){
           xin[i]=sf_cmplx(pz[i], 0.);
           pz[i] = 0.0;
       }
       fftwf_execute(xp);
       for(i=0;i<n;i++) pzz[i] = xout[i];

       ///////////////////////////////////////////////////////// SV-wave's x-component
       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2zz);
       for(jn2=0;jn2<n2zz;jn2++)
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
                xin[i]=rdatazz[ii+ijkz[ikz]]*pxx[i];          
                i++;
              }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2zz;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2zz;jn2++)
              sum2 += fmidzz[im2*n2zz+jn2]*wp[jn2*m+im];

           sum1 += ldatazz[im*m2zz+im2]*sum2;
         }//im2 loop
         px[im] = sum1;
       } 
       free(wp);

       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2xz);
       for(jn2=0;jn2<n2xz;jn2++)
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
                xin[i]=rdataxz[ii+ijkz[ikz]]*pzz[i];          
                i++;
              }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2xz;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2xz;jn2++)
              sum2 += fmidxz[im2*n2xz+jn2]*wp[jn2*m+im];

           sum1 += ldataxz[im*m2xz+im2]*sum2;
         }//im2 loop
         px[im] -= sum1;
       } 
       free(wp);
       ///////////////////////////////////////////////////////// SV-wave's z-component
       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2xx);
       for(jn2=0;jn2<n2xx;jn2++)
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
                xin[i]=rdataxx[ii+ijkz[ikz]]*pzz[i];          
                i++;
              }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2xx;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2xx;jn2++)
              sum2 += fmidxx[im2*n2xx+jn2]*wp[jn2*m+im];

           sum1 += ldataxx[im*m2xx+im2]*sum2;
         }//im2 loop
         pz[im] = sum1;
       } 
       free(wp);

       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2xz);
       for(jn2=0;jn2<n2xz;jn2++)
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
                xin[i]=rdataxz[ii+ijkz[ikz]]*pxx[i];          
                i++;
              }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2xz;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2xz;jn2++)
              sum2 += fmidxz[im2*n2xz+jn2]*wp[jn2*m+im];

           sum1 += ldataxz[im*m2xz+im2]*sum2;
         }//im2 loop
         pz[im] -= sum1;
       } 
       free(wp);

       fftwf_destroy_plan(xp);
       fftwf_destroy_plan(xpi);

       free(pxx);
       free(pzz);
       free(xin);
       free(xout);
#else  // using FFTW in user's own computer
       sf_warning("============= using user installed FFTW ====");
#endif
}

/*****************************************************************************************/
void decomplowrank3dp(float *ldataxx,float *rdataxx,float *fmidxx,
                      float *ldatayy,float *rdatayy,float *fmidyy,
                      float *ldatazz,float *rdatazz,float *fmidzz,
                      float *ldataxy,float *rdataxy,float *fmidxy,
                      float *ldataxz,float *rdataxz,float *fmidxz,
                      float *ldatayz,float *rdatayz,float *fmidyz,
                      float *px, float *py, float *pz, int *ijkx, int *ijky, int *ijkz,
                      int nx, int ny, int nz, int m, int n,
                      int m2xx, int n2xx, int m2yy, int n2yy, int m2zz, int n2zz,
                      int m2xy, int n2xy, int m2xz, int n2xz, int m2yz, int n2yz)
/*< decomplowrank3dp: P-wave vector decomposition based on low-rank decomposition >*/
{
       int   i, im, im2, jn2, ikx, iky, ikz, nxz;
       float sum1, sum2, *wp;

       nxz=nx*nz;

#ifdef SF_HAS_FFTW  // using FFTW in Madagascar
 
       sf_warning("============= using SF_HAS_FFTW ====");

       sf_complex *xin, *xout;
       sf_complex *pxx, *pyy, *pzz; 

       fftwf_plan xp;
       fftwf_plan xpi;

       xin=sf_complexalloc(m);
       xout=sf_complexalloc(n);
       pxx=sf_complexalloc(n);
       pyy=sf_complexalloc(n);
       pzz=sf_complexalloc(n);

       xp=fftwf_plan_dft_3d(ny,nx,nz, (fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_FORWARD,FFTW_ESTIMATE);

       xpi=fftwf_plan_dft_3d(ny,nx,nz,(fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_BACKWARD,FFTW_ESTIMATE);

       /* FFT: from (y,x,z) to (ky, kx, kz) domain */
       for(i=0;i<m;i++){
           xin[i]=sf_cmplx(px[i], 0.);
           px[i] = 0.0;
       }
       fftwf_execute(xp);
       for(i=0;i<n;i++) pxx[i] = xout[i];

       for(i=0;i<m;i++){
           xin[i]=sf_cmplx(py[i], 0.);
           py[i] = 0.0;
       }
       fftwf_execute(xp);
       for(i=0;i<n;i++) pyy[i] = xout[i];

       for(i=0;i<m;i++){
           xin[i]=sf_cmplx(pz[i], 0.);
           pz[i] = 0.0;
       }
       fftwf_execute(xp);
       for(i=0;i<n;i++) pzz[i] = xout[i];

       ///////////////////////////////////////////////////////// P-wave's x-component
       /* n2 IFFT from (ky, kx, kz) to (y, x, z) domain*/
       wp = sf_floatalloc(m*n2xx);
       for(jn2=0;jn2<n2xx;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdataxx[iii+ijkz[ikz]]*pxx[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2xx;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2xx;jn2++)
              sum2 += fmidxx[im2*n2xx+jn2]*wp[jn2*m+im];

           sum1 += ldataxx[im*m2xx+im2]*sum2;
         }//im2 loop
         px[im] = sum1;
       } 
       free(wp);

       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2xy);
       for(jn2=0;jn2<n2xy;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdataxy[iii+ijkz[ikz]]*pyy[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2xy;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2xy;jn2++)
              sum2 += fmidxy[im2*n2xy+jn2]*wp[jn2*m+im];

           sum1 += ldataxy[im*m2xy+im2]*sum2;
         }//im2 loop
         px[im] += sum1;
       } 
       free(wp);

       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2xz);
       for(jn2=0;jn2<n2xz;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdataxz[iii+ijkz[ikz]]*pzz[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2xz;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2xz;jn2++)
              sum2 += fmidxz[im2*n2xz+jn2]*wp[jn2*m+im];

           sum1 += ldataxz[im*m2xz+im2]*sum2;
         }//im2 loop
         px[im] += sum1;
       } 
       free(wp);

       ///////////////////////////////////////////////////////// P-wave's y-component
       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2yy);
       for(jn2=0;jn2<n2yy;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdatayy[iii+ijkz[ikz]]*pyy[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2yy;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2yy;jn2++)
              sum2 += fmidyy[im2*n2yy+jn2]*wp[jn2*m+im];

           sum1 += ldatayy[im*m2yy+im2]*sum2;
         }//im2 loop
         py[im] = sum1;
       } 
       free(wp);

       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2yz);
       for(jn2=0;jn2<n2yz;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdatayz[iii+ijkz[ikz]]*pzz[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2yz;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2yz;jn2++)
              sum2 += fmidyz[im2*n2yz+jn2]*wp[jn2*m+im];

           sum1 += ldatayz[im*m2yz+im2]*sum2;
         }//im2 loop
         py[im] += sum1;
       } 
       free(wp);

       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2xy);
       for(jn2=0;jn2<n2xy;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdataxy[iii+ijkz[ikz]]*pxx[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2xy;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2xy;jn2++)
              sum2 += fmidxy[im2*n2xy+jn2]*wp[jn2*m+im];

           sum1 += ldataxy[im*m2xy+im2]*sum2;
         }//im2 loop
         py[im] += sum1;
       } 
       free(wp);
       ///////////////////////////////////////////////////////// P-wave's z-component
       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2zz);
       for(jn2=0;jn2<n2zz;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdatazz[iii+ijkz[ikz]]*pzz[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2zz;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2zz;jn2++)
              sum2 += fmidzz[im2*n2zz+jn2]*wp[jn2*m+im];

           sum1 += ldatazz[im*m2zz+im2]*sum2;
         }//im2 loop
         pz[im] = sum1;
       } 
       free(wp);

       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2yz);
       for(jn2=0;jn2<n2yz;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdatayz[iii+ijkz[ikz]]*pyy[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2yz;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2yz;jn2++)
              sum2 += fmidyz[im2*n2yz+jn2]*wp[jn2*m+im];

           sum1 += ldatayz[im*m2yz+im2]*sum2;
         }//im2 loop
         pz[im] += sum1;
       } 
       free(wp);

       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2xz);
       for(jn2=0;jn2<n2xz;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdataxz[iii+ijkz[ikz]]*pxx[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2xz;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2xz;jn2++)
              sum2 += fmidxz[im2*n2xz+jn2]*wp[jn2*m+im];

           sum1 += ldataxz[im*m2xz+im2]*sum2;
         }//im2 loop
         pz[im] += sum1;
       } 
       free(wp);

       fftwf_destroy_plan(xp);
       fftwf_destroy_plan(xpi);

       free(pxx);
       free(pyy);
       free(pzz);
       free(xin);
       free(xout);
#else  // using FFTW in user's own computer
       sf_warning("============= using user installed FFTW ====");
#endif
}

/*****************************************************************************************/
void decomplowrank3ds(float *ldataxx,float *rdataxx,float *fmidxx,
                      float *ldatayy,float *rdatayy,float *fmidyy,
                      float *ldatazz,float *rdatazz,float *fmidzz,
                      float *ldataxy,float *rdataxy,float *fmidxy,
                      float *ldataxz,float *rdataxz,float *fmidxz,
                      float *ldatayz,float *rdatayz,float *fmidyz,
                      float *px, float *py, float *pz, int *ijkx, int *ijky, int *ijkz,
                      int nx, int ny, int nz, int m, int n,
                      int m2xx, int n2xx, int m2yy, int n2yy, int m2zz, int n2zz,
                      int m2xy, int n2xy, int m2xz, int n2xz, int m2yz, int n2yz)
/*< decomplowrank3ds: S-wave vector decomposition based on low-rank decomposition >*/
{
       int   i, im, im2, jn2, ikx, iky, ikz, nxz;
       float sum1, sum2, *wp;

       nxz=nx*nz;

#ifdef SF_HAS_FFTW  // using FFTW in Madagascar
 
       sf_warning("============= using SF_HAS_FFTW ====");

       sf_complex *xin, *xout;
       sf_complex *pxx, *pyy, *pzz; 

       fftwf_plan xp;
       fftwf_plan xpi;

       xin=sf_complexalloc(m);
       xout=sf_complexalloc(n);
       pxx=sf_complexalloc(n);
       pyy=sf_complexalloc(n);
       pzz=sf_complexalloc(n);

       xp=fftwf_plan_dft_3d(ny,nx,nz, (fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_FORWARD,FFTW_ESTIMATE);

       xpi=fftwf_plan_dft_3d(ny,nx,nz,(fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_BACKWARD,FFTW_ESTIMATE);

       /* FFT: from (y,x,z) to (ky, kx, kz) domain */
       for(i=0;i<m;i++){
           xin[i]=sf_cmplx(px[i], 0.);
           px[i] = 0.0;
       }
       fftwf_execute(xp);
       for(i=0;i<n;i++) pxx[i] = xout[i];

       for(i=0;i<m;i++){
           xin[i]=sf_cmplx(py[i], 0.);
           py[i] = 0.0;
       }
       fftwf_execute(xp);
       for(i=0;i<n;i++) pyy[i] = xout[i];

       for(i=0;i<m;i++){
           xin[i]=sf_cmplx(pz[i], 0.);
           pz[i] = 0.0;
       }
       fftwf_execute(xp);
       for(i=0;i<n;i++) pzz[i] = xout[i];

       ///////////////////////////////////////////////////////// S-wave's x-component
       /* n2 IFFT from (ky, kx, kz) to (y, x, z) domain*/
       wp = sf_floatalloc(m*n2xx);
       for(jn2=0;jn2<n2xx;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdataxx[iii+ijkz[ikz]]*pxx[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2xx;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2xx;jn2++)
              sum2 += fmidxx[im2*n2xx+jn2]*wp[jn2*m+im];

           sum1 += ldataxx[im*m2xx+im2]*sum2;
         }//im2 loop
         px[im] = sum1;
       } 
       free(wp);

       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2xy);
       for(jn2=0;jn2<n2xy;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdataxy[iii+ijkz[ikz]]*pyy[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2xy;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2xy;jn2++)
              sum2 += fmidxy[im2*n2xy+jn2]*wp[jn2*m+im];

           sum1 += ldataxy[im*m2xy+im2]*sum2;
         }//im2 loop
         px[im] -= sum1;  /* Bxy of qS-wave is negative of Bxy of qP-wave*/
       } 
       free(wp);

       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2xz);
       for(jn2=0;jn2<n2xz;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdataxz[iii+ijkz[ikz]]*pzz[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2xz;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2xz;jn2++)
              sum2 += fmidxz[im2*n2xz+jn2]*wp[jn2*m+im];

           sum1 += ldataxz[im*m2xz+im2]*sum2;
         }//im2 loop
         px[im] -= sum1;
       } 
       free(wp);

       ///////////////////////////////////////////////////////// S-wave's y-component
       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2yy);
       for(jn2=0;jn2<n2yy;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdatayy[iii+ijkz[ikz]]*pyy[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2yy;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2yy;jn2++)
              sum2 += fmidyy[im2*n2yy+jn2]*wp[jn2*m+im];

           sum1 += ldatayy[im*m2yy+im2]*sum2;
         }//im2 loop
         py[im] = sum1;
       } 
       free(wp);

       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2yz);
       for(jn2=0;jn2<n2yz;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdatayz[iii+ijkz[ikz]]*pzz[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2yz;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2yz;jn2++)
              sum2 += fmidyz[im2*n2yz+jn2]*wp[jn2*m+im];

           sum1 += ldatayz[im*m2yz+im2]*sum2;
         }//im2 loop
         py[im] -= sum1;
       } 
       free(wp);

       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2xy);
       for(jn2=0;jn2<n2xy;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdataxy[iii+ijkz[ikz]]*pxx[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2xy;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2xy;jn2++)
              sum2 += fmidxy[im2*n2xy+jn2]*wp[jn2*m+im];

           sum1 += ldataxy[im*m2xy+im2]*sum2;
         }//im2 loop
         py[im] -= sum1;
       } 
       free(wp);
       ///////////////////////////////////////////////////////// S-wave's z-component
       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2zz);
       for(jn2=0;jn2<n2zz;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdatazz[iii+ijkz[ikz]]*pzz[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2zz;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2zz;jn2++)
              sum2 += fmidzz[im2*n2zz+jn2]*wp[jn2*m+im];

           sum1 += ldatazz[im*m2zz+im2]*sum2;
         }//im2 loop
         pz[im] = sum1;
       } 
       free(wp);

       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2yz);
       for(jn2=0;jn2<n2yz;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdatayz[iii+ijkz[ikz]]*pyy[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2yz;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2yz;jn2++)
              sum2 += fmidyz[im2*n2yz+jn2]*wp[jn2*m+im];

           sum1 += ldatayz[im*m2yz+im2]*sum2;
         }//im2 loop
         pz[im] -= sum1;
       } 
       free(wp);

       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2xz);
       for(jn2=0;jn2<n2xz;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdataxz[iii+ijkz[ikz]]*pxx[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2xz;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2xz;jn2++)
              sum2 += fmidxz[im2*n2xz+jn2]*wp[jn2*m+im];

           sum1 += ldataxz[im*m2xz+im2]*sum2;
         }//im2 loop
         pz[im] -= sum1;
       } 
       free(wp);

       fftwf_destroy_plan(xp);
       fftwf_destroy_plan(xpi);

       free(pxx);
       free(pyy);
       free(pzz);
       free(xin);
       free(xout);
#else  // using FFTW in user's own computer
       sf_warning("============= using user installed FFTW ====");
#endif
}

/*****************************************************************************************/
void decomplowrank3dshvti(float *ldataxx,float *rdataxx,float *fmidxx,
                      float *ldatayy,float *rdatayy,float *fmidyy,
                      float *ldataxy,float *rdataxy,float *fmidxy,
                      float *px, float *py, int *ijkx, int *ijky, int *ijkz,
                      int nx, int ny, int nz, int m, int n,
                      int m2xx, int n2xx, int m2yy, int n2yy, int m2zz, int n2zz,
                      int m2xy, int n2xy)
/*< decomplowrank3dshvti: VTI SH-wave vector decomposition based on low-rank decomposition >*/
{
       int   i, im, im2, jn2, ikx, iky, ikz, nxz;
       float sum1, sum2, *wp;

       nxz=nx*nz;

#ifdef SF_HAS_FFTW  // using FFTW in Madagascar
 
       sf_warning("============= using SF_HAS_FFTW ====");

       sf_complex *xin, *xout;
       sf_complex *pxx, *pyy; 

       fftwf_plan xp;
       fftwf_plan xpi;

       xin=sf_complexalloc(m);
       xout=sf_complexalloc(n);
       pxx=sf_complexalloc(n);
       pyy=sf_complexalloc(n);

       xp=fftwf_plan_dft_3d(ny,nx,nz, (fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_FORWARD,FFTW_ESTIMATE);

       xpi=fftwf_plan_dft_3d(ny,nx,nz,(fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_BACKWARD,FFTW_ESTIMATE);

       /* FFT: from (y,x,z) to (ky, kx, kz) domain */
       for(i=0;i<m;i++){
           xin[i]=sf_cmplx(px[i], 0.);
           px[i] = 0.0;
       }
       fftwf_execute(xp);
       for(i=0;i<n;i++) pxx[i] = xout[i];

       for(i=0;i<m;i++){
           xin[i]=sf_cmplx(py[i], 0.);
           py[i] = 0.0;
       }
       fftwf_execute(xp);
       for(i=0;i<n;i++) pyy[i] = xout[i];

       ///////////////////////////////////////////////////////// P-wave's x-component
       /* n2 IFFT from (ky, kx, kz) to (y, x, z) domain*/
       wp = sf_floatalloc(m*n2xx);
       for(jn2=0;jn2<n2xx;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdataxx[iii+ijkz[ikz]]*pxx[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2xx;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2xx;jn2++)
              sum2 += fmidxx[im2*n2xx+jn2]*wp[jn2*m+im];

           sum1 += ldataxx[im*m2xx+im2]*sum2;
         }//im2 loop
         px[im] = sum1;
       } 
       free(wp);

       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2xy);
       for(jn2=0;jn2<n2xy;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdataxy[iii+ijkz[ikz]]*pyy[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2xy;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2xy;jn2++)
              sum2 += fmidxy[im2*n2xy+jn2]*wp[jn2*m+im];

           sum1 += ldataxy[im*m2xy+im2]*sum2;
         }//im2 loop
         px[im] += sum1;
       } 
       free(wp);

       ///////////////////////////////////////////////////////// P-wave's y-component
       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2yy);
       for(jn2=0;jn2<n2yy;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdatayy[iii+ijkz[ikz]]*pyy[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2yy;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2yy;jn2++)
              sum2 += fmidyy[im2*n2yy+jn2]*wp[jn2*m+im];

           sum1 += ldatayy[im*m2yy+im2]*sum2;
         }//im2 loop
         py[im] = sum1;
       } 
       free(wp);

       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       wp = sf_floatalloc(m*n2xy);
       for(jn2=0;jn2<n2xy;jn2++)
       {
           i=0;
           int jn2n=jn2*n;
           for(iky=0;iky<ny;iky++)
           {
              int iynxz=ijky[iky]*nxz;
              int ii=jn2n+iynxz;
              for(ikx=0;ikx<nx;ikx++)
              {
                 int ixnz=ijkx[ikx]*nz;
                 int iii=ii+ixnz;
                 for(ikz=0;ikz<nz;ikz++)
                 {
                   xin[i]=rdataxy[iii+ijkz[ikz]]*pxx[i];          
                   i++;
                 }
               }
            }
            // (kx,kz) to (x, z) domain
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }

       // Matrix multiplication in space-domain 
       for(im=0;im<m;im++)
       {
         sum1=0.0;
         for(im2=0;im2<m2xy;im2++)
         {
           sum2=0.0;
           for(jn2=0;jn2<n2xy;jn2++)
              sum2 += fmidxy[im2*n2xy+jn2]*wp[jn2*m+im];

           sum1 += ldataxy[im*m2xy+im2]*sum2;
         }//im2 loop
         py[im] += sum1;
       } 
       free(wp);

       fftwf_destroy_plan(xp);
       fftwf_destroy_plan(xpi);

       free(pxx);
       free(pyy);
       free(xin);
       free(xout);
#else  // using FFTW in user's own computer
       sf_warning("============= using user installed FFTW ====");
#endif
}

