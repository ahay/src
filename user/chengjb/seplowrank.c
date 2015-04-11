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
void seplowrank2d(float *ldata,float *rdata,float *fmid, float *x, int *ijkx, int *ijkz,
                int nx,int nz,int m,int n,int m2,int n2, int iflag)
/*< seplowrank2d: separating wave-modes based on low-rank decomposition >*/
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

#ifdef SF_HAS_FFTW  /* using FFTW in Madagascar */

       sf_complex *xx, *xin, *xout;

       fftwf_plan xp;
       fftwf_plan xpi;

       xin=sf_complexalloc(m);
       xout=sf_complexalloc(n);
       xx=sf_complexalloc(n);

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
           
       for(i=0;i<n;i++) xx[i] = xout[i];

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
                xin[i]=rdata[ii+ijkz[ikz]]*xx[i];          
                i++;
              }
            }
            /* (kx,kz) to (x, z) domain */
            fftwf_execute(xpi);

            for(im=0;im<m;im++)
                wp[jn2*m+im] = creal(xout[im])/n;
       }
       fftwf_destroy_plan(xp);
       fftwf_destroy_plan(xpi);
       free(xx);
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
         }/*im2 loop*/
         x[im] = sum1;
       } 

#else  /* using FFTW in user's own computer */

       fftw_complex *xx, *xin, *xout;

       fftw_plan xp;
       fftw_plan xpi;

       xin=fftw_complexalloc(m);
       xout=fftw_complexalloc(n);
       xx=fftw_complexalloc(n);

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
           
       if(iflag!=1) for(i=0;i<n;i++) xout[i] *= sf_cmplx(0.0, 1.0);

       for(i=0;i<n;i++) xx[i] = xout[i];

       /* n2 IFFT from (kx, kz) to (x, z) domain*/
       for(jn2=0;jn2<n2;jn2++)
       {
           int jn2n=jn2*n;
           i=0;
           for(ikx=0;ikx<nx;ikx++)
           {
              /* Note: Spectrum of the operator is differently orderred as the spectrum after FFT */ 
              int ixnz=ijkx[ikx]*nz;
              int ii=jn2n+ixnz;
              for(ikz=0;ikz<nz;ikz++)
              {
                xin[i]=rdata[ii+ijkz[ikz]]*xx[i];          
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
       free(xx);
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
         }/*im2 loop*/
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

#ifdef SF_HAS_FFTW  /* using FFTW in Madagascar */

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

#else  /* using FFTW in user's own computer */

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

void  reconstruct(float **w, float *ldata, float *fmid, float *rdata, int m, int n, int m2, int n2)
/*< re-construct matrix using the lowrank decomposed matrixs >*/
{
     int  im, in, im2, in2;
     float sum1, sum2;

     for(im=0;im<m;im++)
     for(in=0;in<n;in++)
     {
        sum1=0.0;
        for(im2=0;im2<m2;im2++)
        {

           sum2=0.0;
           for(in2=0;in2<n2;in2++)
              sum2 += fmid[im2*n2+in2]*rdata[in2*n+in];
           sum1 += ldata[im*m2+im2]*sum2;
        }
        w[im][in]=sum1;
     }
}

void  reconstruct1(float *w, float *ldata, float *fmid, float *rdata, int m, int n, int m2, int n2, int im)
/*< re-construct matrix using the lowrank decomposed matrixs >*/
{
     int   in, im2, in2;
     float sum1, sum2;

     for(in=0;in<n;in++)
     {
        sum1=0.0;
        for(im2=0;im2<m2;im2++)
        {

           sum2=0.0;
           for(in2=0;in2<n2;in2++)
              sum2 += fmid[im2*n2+in2]*rdata[in2*n+in];
           sum1 += ldata[im*m2+im2]*sum2;
        }
        w[in]=sum1;
     }
}

/*****************************************************************************************/
void seplowrank3d(float *ldata,float *rdata,float *fmid, float *p, int *ijkx, int *ijky, int *ijkz,
                      int nx, int ny, int nz, int m, int n, int m2, int n2, int iflag)
/*< seplowrank3d: wave-mode separation based on low-rank decomposition >*/
{
       int i, im, im2, jn2, ikx, iky, ikz, nxz;
       float sum1, sum2, *wp;

       wp = sf_floatalloc(m*n2);

       nxz=nx*nz;

#ifdef SF_HAS_FFTW  /* using FFTW in Madagascar */
 
       sf_complex *xx, *xin, *xout;

       fftwf_plan xp;
       fftwf_plan xpi;

       xin=sf_complexalloc(m);
       xout=sf_complexalloc(n);
       xx=sf_complexalloc(n);

       xp=fftwf_plan_dft_3d(ny,nx,nz, (fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_FORWARD,FFTW_ESTIMATE);

       xpi=fftwf_plan_dft_3d(ny,nx,nz,(fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_BACKWARD,FFTW_ESTIMATE);

       /* FFT: from (x,z) to (kx, kz) domain */

       if(iflag==1)
           for(i=0;i<m;i++) xin[i]=sf_cmplx(p[i], 0.);
       else 
           for(i=0;i<m;i++) xin[i]=sf_cmplx(0.0, p[i]);

       fftwf_execute(xp);
           
       for(i=0;i<n;i++) xx[i] = xout[i];

       /* n2 IFFT from (kx, kz) to (x, z) domain */
       for(jn2=0;jn2<n2;jn2++)
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
                   xin[i]=rdata[iii+ijkz[ikz]]*xx[i];
                   i++;
                 }
             }
        }

       /* (kx,kz) to (x, z) domain */
       fftwf_execute(xpi);

       for(im=0;im<m;im++)
           wp[jn2*m+im] = creal(xout[im])/n;
    }

   fftwf_destroy_plan(xp);
   fftwf_destroy_plan(xpi);

   free(xx);
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
        }/*im2 loop*/
        p[im] = sum1;
  } 

#endif

  free(wp);
}
