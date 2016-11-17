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

#ifdef SF_HAS_FFTW
#include <fftw3.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

/*****************************************************************************************/
void fwpvti2delowrank(float *ldata,float *rdata,float *fmid, float *y, float *x, int *ijkx, int *ijkz,
                      int nx,int nz,int m,int n,int m2,int n2)
/*< fwpvti2delowrank: apply low-rank decomposed propagator to the wavefield component >*/
{
    int i, im, im2, jn2, ikx, ikz;
    float sum1, sum2, *wp;

    wp = sf_floatalloc(m*n2);

#ifdef SF_HAS_FFTW  // using FFTW in Madagascar
 
    //sf_warning("============= using SF_HAS_FFTW ====");

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
    for(i=0;i<m;i++) xin[i]=sf_cmplx(x[i], 0.);

    fftwf_execute(xp);
           
    for(i=0;i<n;i++) xx[i] = xout[i];

    /* n2 IFFT from (kx, kz) to (x, z) domain*/
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
                xin[i]=rdata[ii+ijkz[ikz]]*xx[i];          
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
	y[im] = sum1;
    } 
#endif

    free(wp);
}

/*****************************************************************************************/
void fwpvti2delowranksvd(float *ldata,float *rdata,float *fmid, float *y, float *x, int *ijkx, int *ijkz,
			 int nx,int nz,int m,int n,int m2,int n2, float kxm, float kzm, float kxzm, float *akx, float *akz)
/*< fwpvti2delowranksvd: apply low-rank decomposed propagator considering stability to the wavefield component >*/
{
    int i, im, im2, jn2, ikx, ikz;
    float sum1, sum2, *wp;
    sf_complex *xx, *xin, *xout;
    float kx, kz, kxz;

    wp = sf_floatalloc(m*n2);

#ifdef SF_HAS_FFTW  // using FFTW in Madagascar

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
    for(i=0;i<m;i++) xin[i]=sf_cmplx(x[i], 0.);

    fftwf_execute(xp);
           
    for(i=0;i<n;i++) xx[i] = xout[i];

    /* n2 IFFT from (kx, kz) to (x, z) domain*/
    int index;
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
	y[im] = sum1;
    } 

#endif

    free(wp);
}

/*****************************************************************************************/
void fwpvti2delowranksvdkspace(float *ldata,float *rdata,float *fmid, float *y, float *x, int *ijkx, int *ijkz,
			       int nx,int nz,int m,int n,int m2,int n2, float kxm, float kzm, float kxzm, float *akx, float *akz,
			       sf_complex *source, float amps)
/*< fwpvti2delowranksvdkspace: apply k-space low-rank decomposed propagator considering stability to the wavefield component >*/
{
    int i, im, im2, jn2, ikx, ikz;
    float sum1, sum2, *wp;
    float kx, kz, kxz;

    wp = sf_floatalloc(m*n2);

    sf_complex *xx, *xin, *xout;

#ifdef SF_HAS_FFTW  // using FFTW in Madagascar

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
    for(i=0;i<m;i++) xin[i]=sf_cmplx(x[i], 0.);

    fftwf_execute(xp);
           
    //for(i=0;i<n;i++) xx[i] = xout[i];
	   
    int ixnz, index;

    i=0;
    for(ikx=0;ikx<nx;ikx++)
    {
	ixnz=ijkx[ikx]*nz;
	for(ikz=0;ikz<nz;ikz++)
	{
	    index = ixnz + ijkz[ikz];
	    xx[i] = xout[i] + amps*source[index];
	    i++;
	}
    }

    /* n2 IFFT from (kx, kz) to (x, z) domain*/
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
	y[im] = sum1;
    } 
#endif
}

/*****************************************************************************************/
void fwpvti2despectal(float *y, float *x, float *C13C44, float *kxkz, int *ijkx, int *ijkz,
                      int nx,int nz, int m,int n, float dt2)
/*< fwpvti2despectral: apply spectral propagator to the wavefield component >*/
{
    int i, im, ikx, ikz;

#ifdef SF_HAS_FFTW  // using FFTW in Madagascar
 
    //sf_warning("============= using SF_HAS_FFTW ====");

    sf_complex *xx, *xin, *xout;

    fftwf_plan xp;
    fftwf_plan xpi;

    xin=sf_complexalloc(m);
    xout=sf_complexalloc(n);

    xp=fftwf_plan_dft_2d(nx,nz, (fftwf_complex *) xin, (fftwf_complex *) xout,
			 FFTW_FORWARD,FFTW_ESTIMATE);

    xpi=fftwf_plan_dft_2d(nx,nz,(fftwf_complex *) xin, (fftwf_complex *) xout,
			  FFTW_BACKWARD,FFTW_ESTIMATE);

    /* FFT: from (x,z) to (kx, kz) domain */
    for(i=0;i<m;i++) xin[i]=sf_cmplx(x[i], 0.);

    fftwf_execute(xp);
           
    i=0;
    for(ikx=0;ikx<nx;ikx++)
    {
	int ixnz=ijkx[ikx]*nz;
	for(ikz=0;ikz<nz;ikz++)
	{
	    xin[i]=kxkz[ixnz+ijkz[ikz]]*xout[i];          
	    i++;
	}
    }

    // (kx,kz) to (x, z) domain
    fftwf_execute(xpi);

    for(im=0;im<m;im++)
	y[im] = dt2*C13C44[im]*creal(xout[im])/n;
       
    fftwf_destroy_plan(xp);
    fftwf_destroy_plan(xpi);
    free(xin);
    free(xout);

#endif
}

/*****************************************************************************************/
void fwpvti3delowrank(float *ldata,float *rdata,float *fmid, float *y, float *x, int *ijkx, int *ijky, int *ijkz,
                      int nx,int ny, int nz,int m,int n,int m2,int n2)
/*< fwpvti3delowrank: apply low-rank decomposed propagator to the wavefield component >*/
{
    int i, im, im2, jn2, ikx, iky, ikz, nxz;
    float sum1, sum2, *wp;

    nxz = nx*nz;

    wp = sf_floatalloc(m*n2);

#ifdef SF_HAS_FFTW  // using FFTW in Madagascar
 
    //  sf_warning("============= using SF_HAS_FFTW ====");

    sf_complex *xx, *xin, *xout;

    fftwf_plan xp;
    fftwf_plan xpi;

    xin=sf_complexalloc(m);
    xout=sf_complexalloc(n);
    xx=sf_complexalloc(n);

    fftwf_plan_with_nthreads(omp_get_max_threads());
    xp=fftwf_plan_dft_3d(ny,nx,nz, (fftwf_complex *) xin, (fftwf_complex *) xout,
			 FFTW_FORWARD,FFTW_ESTIMATE);

    fftwf_plan_with_nthreads(omp_get_max_threads());
    xpi=fftwf_plan_dft_3d(ny,nx,nz,(fftwf_complex *) xin, (fftwf_complex *) xout,
			  FFTW_BACKWARD,FFTW_ESTIMATE);

    /* FFT: from (y,x,z) to (ky, kx, kz) domain */
    for(i=0;i<m;i++) xin[i]=sf_cmplx(x[i], 0.);

    fftwf_execute(xp);
           
    for(i=0;i<n;i++) xx[i] = xout[i];

    /* n2 IFFT from (ky, kx, kz) to (y, x, z) domain*/
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
	y[im] = sum1;
    } 

#endif
    free(wp);
}

/*****************************************************************************************/
void fwpvti3delowrank_double(double *ldata,double *rdata,double *fmid, double *y, double *x, int *ijkx, int *ijky, int *ijkz,
			     int nx,int ny, int nz,int m,int n,int m2,int n2)
/*< fwpvti3delowrank: apply low-rank decomposed propagator to the wavefield component >*/
{
    int i, im, im2, jn2, ikx, iky, ikz, nxz;
    long double sum1, sum2, *wp;

    nxz = nx*nz;

    wp = (long double*)malloc(sizeof(long double)*m*n2);

#ifdef SF_HAS_FFTW  // using FFTW in Madagascar

    fftw_complex *xx, *xin, *xout;

    fftw_plan xp;
    fftw_plan xpi;

    xin=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m);
    xout=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);
    xx=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);

    fftw_plan_with_nthreads(omp_get_max_threads());
    xp=fftw_plan_dft_3d(ny,nx,nz, (fftw_complex *) xin, (fftw_complex *) xout,
			FFTW_FORWARD,FFTW_ESTIMATE);

    fftw_plan_with_nthreads(omp_get_max_threads());
    xpi=fftw_plan_dft_3d(ny,nx,nz,(fftw_complex *) xin, (fftw_complex *) xout,
			 FFTW_BACKWARD,FFTW_ESTIMATE);

    /* FFT: from (y,x,z) to (ky, kx, kz) domain */
    for(i=0;i<m;i++) xin[i]=x[i];

    fftw_execute(xp);
           
    for(i=0;i<n;i++) xx[i] = xout[i];

    /* n2 IFFT from (ky, kx, kz) to (y, x, z) domain*/
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
	y[im] = sum1;
    } 

#endif

    free(wp);
}

