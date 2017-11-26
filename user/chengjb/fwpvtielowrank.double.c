/* apply low-rank propagators to wavefield component */
/*
  Copyright (C) 2014 Tongji University (Jiubing Cheng) 
  and King Abdulah University of Science and Technology (Zedong Wu and Tariq Alkhalifah)
 
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

#include <complex.h>

#ifdef SF_HAS_FFTW
#include <fftw3.h>
#endif

void fwpvti2delowranksvd_double(double *ldata,double *rdata,double *fmid, double *y, double *x, int *ijkx, int *ijkz,
				int nx,int nz,int m,int n,int m2,int n2, float kxm, float kzm, float kxzm, double *akx, double *akz)
/*< fwpvti2delowranksvd_double: apply low-rank decomposed propagator considering stability to the wavefield component >*/
{
    int i, im, im2, jn2, ikx, ikz;
    long double sum1, sum2,temp1,temp2;
    long double *wp;

    wp = (long double*)malloc(sizeof(long double)*m*n2);

#ifdef SF_HAS_FFTW
    fftw_complex *xx, *xin, *xout;
    double kx, kz, kxz;


    fftw_plan xp;
    fftw_plan xpi;

    xin=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m);
    xout=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);
    xx=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);

    xp=fftw_plan_dft_2d(nx,nz, (fftw_complex *) xin, (fftw_complex *) xout,
			FFTW_FORWARD,FFTW_ESTIMATE);

    xpi=fftw_plan_dft_2d(nx,nz,(fftw_complex *) xin, (fftw_complex *) xout,
			 FFTW_BACKWARD,FFTW_ESTIMATE);

    /* FFT: from (x,z) to (kx, kz) domain */
    for(i=0;i<m;i++) xin[i]=x[i];

    fftw_execute(xp);
           
    for(i=0;i<n;i++) xx[i] = xout[i];

    /* n2 IFFT from (kx, kz) to (x, z) domain*/
    int index,jn2n,ixnz,ii;

/*
  #pragma omp parallel for private(jn2,ikx,im,ikz,i,ixnz,ii,index,kx,kz,kxz) \
  schedule(dynamic) \
  shared(akx, akz, rdata, wp, xin, xx, xout, ijkx, ijkz,  nx, nz )
*/
    for(jn2=0;jn2<n2;jn2++)
    {
	jn2n=jn2*n;
	for(ikx=0;ikx<nx;ikx++)
	{
	    // Note: Spectrum of the operator is differently orderred as the spectrum after FFT
	    ixnz=ijkx[ikx]*nz;
	    ii=jn2n+ixnz;

	    for(ikz=0;ikz<nz;ikz++)
	    {
		i = ikz + ikx*nz;
		index = ixnz + ijkz[ikz];
		kx = akx[index]; 
		kz = akz[index]; 
		kxz = sqrt(kx*kx+kz*kz);

		if (kx>kxm || kz>kzm || kxz>kxzm)
		    xin[i] = 0.0;
		else if (kx<-kxm || kz<-kzm)
		    xin[i] = 0.0;
		else
		    xin[i]=rdata[ii+ijkz[ikz]]*xx[i];          
               
	    }
	}
	// (kx,kz) to (x, z) domain
	fftw_execute(xpi);

	for(im=0;im<m;im++)
	    wp[jn2*m+im] = creal(xout[im])/n;
    }

    fftw_destroy_plan(xp);
    fftw_destroy_plan(xpi);
    free(xx);
    free(xin);
    free(xout);

    /* Matrix multiplication in space-domain */
/* #pragma omp parallel for private(im,im2,jn2,temp1,temp2,sum1,sum2)	\
   schedule(dynamic) \
   shared(wp, fmid, ldata, y, m, m2, n2) */
    for(im=0;im<m;im++)
    {
	sum1=0.0;
	for(im2=0;im2<m2;im2++)
	{
	    sum2=0.0;
	    for(jn2=0;jn2<n2;jn2++)
	    {
		temp1 = fmid[im2*n2+jn2];
		sum2 += temp1*wp[jn2*m+im];
	    }

	    temp2 = ldata[im*m2+im2];
	    sum1 += temp2*sum2;
	}//im2 loop
	y[im] = (double)sum1;
    } 

#endif

    free(wp);
}

/*****************************************************************************************/

void fwpvti2delowranksvdkspace_double(double *ldata,double *rdata,double *fmid, double *y, double *x, int *ijkx, int *ijkz,
				      int nx,int nz,int m,int n,int m2,int n2, float kxm, float kzm, float kxzm, double *akx, double *akz,
				      sf_complex *source, float amps)
/*< fwpvti2delowranksvdkspace_double: apply k-space low-rank decomposed propagator considering stability to the wavefield component >*/
{
    int i, im, im2, jn2, ikx, ikz;
    long double sum1, sum2, *wp;
    double kx, kz, kxz;

    wp = (long double*)malloc(sizeof(long double)*m*n2);

#ifdef SF_HAS_FFTW       

    fftw_complex *xx, *xin, *xout;
    xin=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m);
    xout=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);
    xx=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);

    fftw_plan xp;
    fftw_plan xpi;

    xp=fftw_plan_dft_2d(nx,nz, (fftw_complex *) xin, (fftw_complex *) xout,
			FFTW_FORWARD,FFTW_ESTIMATE);

    xpi=fftw_plan_dft_2d(nx,nz,(fftw_complex *) xin, (fftw_complex *) xout,
			 FFTW_BACKWARD,FFTW_ESTIMATE);

    /* FFT: from (x,z) to (kx, kz) domain */
    for(i=0;i<m;i++) xin[i]=sf_cmplx(x[i], 0.);

    fftw_execute(xp);

    //for(i=0;i<n;i++) xx[i] = xout[i];

    int ixnz, index, jn2n, ii;
 
/* #pragma omp parallel for private(ikx,ixnz,ikz,index)	\
   schedule(dynamic) \
   shared(ijkx,ijkz, source, xx, xout, nx, nz, amps) */
    for(ikx=0;ikx<nx;ikx++)
    {
	ixnz=ijkx[ikx]*nz;
	for(ikz=0;ikz<nz;ikz++)
	{
	    index = ixnz + ijkz[ikz];
	    i = ikz+ikx*nz;
	    xx[i] = xout[i] + amps*source[index];
	}
    }

    /* n2 IFFT from (kx, kz) to (x, z) domain*/
/* #pragma omp parallel for private(jn2,jn2n,ikx,ikz,i,ixnz,ii,index,kx,kz,kxz) \
   schedule(dynamic) \
   shared(akx, akz, rdata, wp, xin, xx, xout, ijkx, ijkz,  nx, nz, kxm, kzm, kxzm) */
    for(jn2=0;jn2<n2;jn2++)
    {
	i=0;
	jn2n=jn2*n;
	for(ikx=0;ikx<nx;ikx++)
	{
	    // Note: Spectrum of the operator is differently orderred as the spectrum after FFT
	    ixnz=ijkx[ikx]*nz;
	    ii=jn2n+ixnz;

	    for(ikz=0;ikz<nz;ikz++)
	    {
                index = ixnz + ijkz[ikz];
                kx = akx[index];
                kz = akz[index];
                kxz = sqrt(kx*kx+kz*kz);

                if (kx>kxm || kz>kzm || kxz>kxzm)
                    xin[i] = 0.0;
                else if (kx<-kxm || kz<-kzm)
                    xin[i] = 0.0;
                else
                    xin[i]=rdata[ii+ijkz[ikz]]*xx[i];
                i++;
	    }
	}
	// (kx,kz) to (x, z) domain
	fftw_execute(xpi);

	for(im=0;im<m;im++)
	    wp[jn2*m+im] = creal(xout[im])/n;
    }
    printf("End omp!\n");

    fftw_destroy_plan(xp);
    fftw_destroy_plan(xpi);
    free(xx);
    free(xin);
    free(xout);

    /* Matrix multiplication in space-domain */
/* #pragma omp parallel for private(im,im2,jn2,sum1,sum2)	\
   schedule(dynamic) \
   shared(wp, fmid, ldata, y, m, m2, n2) */
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
	y[im] = sum1;
    }

#endif

}

