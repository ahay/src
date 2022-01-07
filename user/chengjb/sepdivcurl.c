/* separating wave-modes based on divergence and curl operations*/
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

void sepdiv2d(float *rk, float *x, int *ijkx, int *ijkz, int nx,int nz,int m,int n, int iflag)
/*< sepdiv2d: separating wave-modes based on divergence >*/
{
       int i, im, ikx, ikz;

#ifdef SF_HAS_FFTW  // using FFTW in Madagascar
 
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
           
       /* IFFT from (kx, kz) to (x, z) domain*/
	   i=0;
       for(ikx=0;ikx<nx;ikx++)
       {
          // Note: Spectrum of the operator is differently orderred as the spectrum after FFT
          int ixnz=ijkx[ikx]*nz;
          for(ikz=0;ikz<nz;ikz++)
          {
                xin[i]=rk[ixnz+ijkz[ikz]]*xout[i];          
                i++;
          }
	   }
       // (kx,kz) to (x, z) domain
       fftwf_execute(xpi);
       for(im=0;im<m;im++)
          x[im] = creal(xout[im])/n;

       fftwf_destroy_plan(xp);
       fftwf_destroy_plan(xpi);
       free(xin);
       free(xout);

#else  // using FFTW in user's own computer
       //sf_warning("============= using user installed FFTW ====");

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
           
       if(iflag!=1) for(i=0;i<n;i++) xout[i] *= sf_cmplx(0.0, 1.0);

       /* IFFT from (kx, kz) to (x, z) domain*/
       i=0;
       for(ikx=0;ikx<nx;ikx++)
       {
              /* Note: Spectrum of the operator is differently orderred as the spectrum after FFT */ 
              int ixnz=ijkx[ikx]*nz;
              for(ikz=0;ikz<nz;ikz++)
              {
                xin[i]=rk[ixnz+ijkz[ikz]]*xout[i];          
                i++;
              }
	   }
       /* (kx,kz) to (x, z) domain */
       fftw_execute(xpi);

       for(im=0;im<m;im++)
          x[im] = xout[im][0]/n;
      
       fftw_destroy_plan(xp);
       fftw_destroy_plan(xpi);
       free(xin);
       free(xout);

#endif
}

void sepdiv3d(float *rk, float *x, int *ijkx, int *ijky, int *ijkz, int nx,int ny,int nz,int m,int n, int iflag)
/*< sepdiv3d: separating wave-modes based on divergence >*/
{
       int i, im, ikx, iky, ikz, nxz;

	   nxz=nx*nz;
#ifdef SF_HAS_FFTW  // using FFTW in Madagascar
 
       sf_complex *xin, *xout;

       fftwf_plan xp;
       fftwf_plan xpi;

       xin=sf_complexalloc(m);
       xout=sf_complexalloc(n);

       xp=fftwf_plan_dft_3d(ny,nx,nz, (fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_FORWARD,FFTW_ESTIMATE);

       xpi=fftwf_plan_dft_3d(ny,nx,nz,(fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_BACKWARD,FFTW_ESTIMATE);

       /* FFT: from (y,x,z) to (ky, kx, kz) domain */

       if(iflag==1)
           for(i=0;i<m;i++) xin[i]=sf_cmplx(x[i], 0.);
       else 
           for(i=0;i<m;i++) xin[i]=sf_cmplx(0.0, x[i]);

       fftwf_execute(xp);
           
       /* IFFT from (ky, kx, kz) to (y, x, z) domain*/
       // Note: Spectrum of the operator is differently orderred as the spectrum after FFT
	   i=0;
       for(iky=0;iky<ny;iky++)
       {
          int iyxnz=ijky[iky]*nxz;
          for(ikx=0;ikx<nx;ikx++)
          {
             int ixnz=iyxnz+ijkx[ikx]*nz;
             for(ikz=0;ikz<nz;ikz++)
             {
                xin[i]=rk[ixnz+ijkz[ikz]]*xout[i];          
                i++;
             }
		  }
	   }
       // (kx,kz) to (x, z) domain
       fftwf_execute(xpi);
       for(im=0;im<m;im++)
          x[im] = creal(xout[im])/n;

       fftwf_destroy_plan(xp);
       fftwf_destroy_plan(xpi);
       free(xin);
       free(xout);

#else  // using FFTW in user's own computer
       //sf_warning("============= using user installed FFTW ====");

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
           
       if(iflag!=1) for(i=0;i<n;i++) xout[i] *= sf_cmplx(0.0, 1.0);

       /* Note: Spectrum of the operator is differently orderred as the spectrum after FFT */ 
       /* IFFT from (ky, kx, kz) to (y, x, z) domain*/
	   i=0;
       for(iky=0;iky<ny;iky++)
       {
          int iyxnz=ijky[iky]*nxz;
          for(ikx=0;ikx<nx;ikx++)
          {
              int ixnz=iyxnz+ijkx[ikx]*nz;
              for(ikz=0;ikz<nz;ikz++)
              {
                xin[i]=rk[ixnz+ijkz[ikz]]*xout[i];          
                i++;
              }
		  }
	   }
       /* (kx,kz) to (x, z) domain */
       fftw_execute(xpi);

       for(im=0;im<m;im++)
          x[im] = xout[im][0]/n;
      
       fftw_destroy_plan(xp);
       fftw_destroy_plan(xpi);
       free(xin);
       free(xout);

#endif
}

void sepdiv3dD(double *rk, float *x, int *ijkx, int *ijky, int *ijkz, int nx,int ny,int nz,int m,int n, int iflag)
/*< sepdiv3d: separating wave-modes based on divergence >*/
{
       int i, im, ikx, iky, ikz, nxz;

	   nxz=nx*nz;
#ifdef SF_HAS_FFTW  // using FFTW in Madagascar
 
       sf_complex *xin, *xout;

       fftwf_plan xp;
       fftwf_plan xpi;

       xin=sf_complexalloc(m);
       xout=sf_complexalloc(n);

       xp=fftwf_plan_dft_3d(ny,nx,nz, (fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_FORWARD,FFTW_ESTIMATE);

       xpi=fftwf_plan_dft_3d(ny,nx,nz,(fftwf_complex *) xin, (fftwf_complex *) xout,
			    FFTW_BACKWARD,FFTW_ESTIMATE);

       /* FFT: from (y,x,z) to (ky, kx, kz) domain */

       if(iflag==1)
           for(i=0;i<m;i++) xin[i]=sf_cmplx(x[i], 0.);
       else 
           for(i=0;i<m;i++) xin[i]=sf_cmplx(0.0, x[i]);

       fftwf_execute(xp);
           
       /* IFFT from (ky, kx, kz) to (y, x, z) domain*/
       // Note: Spectrum of the operator is differently orderred as the spectrum after FFT
	   i=0;
       for(iky=0;iky<ny;iky++)
       {
          int iyxnz=ijky[iky]*nxz;
          for(ikx=0;ikx<nx;ikx++)
          {
             int ixnz=iyxnz+ijkx[ikx]*nz;
             for(ikz=0;ikz<nz;ikz++)
             {
                xin[i]=(float)rk[ixnz+ijkz[ikz]]*xout[i];          
                i++;
             }
		  }
	   }
       // (kx,kz) to (x, z) domain
       fftwf_execute(xpi);
       for(im=0;im<m;im++)
          x[im] = creal(xout[im])/n;

       fftwf_destroy_plan(xp);
       fftwf_destroy_plan(xpi);
       free(xin);
       free(xout);

#else  // using FFTW in user's own computer
       //sf_warning("============= using user installed FFTW ====");

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
           
       if(iflag!=1) for(i=0;i<n;i++) xout[i] *= sf_cmplx(0.0, 1.0);

       /* Note: Spectrum of the operator is differently orderred as the spectrum after FFT */ 
       /* IFFT from (ky, kx, kz) to (y, x, z) domain*/
	   i=0;
       for(iky=0;iky<ny;iky++)
       {
          int iyxnz=ijky[iky]*nxz;
          for(ikx=0;ikx<nx;ikx++)
          {
              int ixnz=iyxnz+ijkx[ikx]*nz;
              for(ikz=0;ikz<nz;ikz++)
              {
                xin[i]=(float)rk[ixnz+ijkz[ikz]]*xout[i];          
                i++;
              }
		  }
	   }
       /* (kx,kz) to (x, z) domain */
       fftw_execute(xpi);

       for(im=0;im<m;im++)
          x[im] = xout[im][0]/n;
      
       fftw_destroy_plan(xp);
       fftw_destroy_plan(xpi);
       free(xin);
       free(xout);

#endif
}

