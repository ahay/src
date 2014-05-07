/* Linear/parabolic radon transform frequency domain implementation 
Also referred to as high-resolution radon transform
Note: I borrowed a lot from /system/seismic/radon+Mradon.c. The distinction:
	I am using FFTW because I am inexperienced in invoking kiss_fft. 
*/
/*
  Copyright (C) 2014  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

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

  References: 
	1) Kostov C., 1990. "Toeplitz structure in Slant-Stack inversion":
	 SEG Extended Abstracts, 1647-1650.
	2) Sacchi, Mauricio D., and Milton Porsani. "Fast high resolution 
	parabolic Radon transform." Society of Exploration Geophysicists 
	69th Annual International Meeting, SPRO P. Vol. 1. No. 1. 1999.
	3) Vogel, Curtis R. Computational methods for inverse problems. 
	Vol. 23. Siam, 2002.	[Chapter 5.2.]
*/

#include <rsf.h>
#include <time.h>
#include <complex.h>
#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "myradon2.h"
#include "radonlop.h"


static bool inv;
static int nt, np, nx, nw, nfft;
static float x0, dt, dp, eps;
static float *p, *xx, *xxtmp;
static float **dd, **mm, *tmpr;
static sf_complex *cdd, *cmm;
fftwf_complex *tmpc;
fftwf_plan fft1, ifft1;

void sf_radon2_init(bool inv_, float *p_, float *xx_, int nt_, int np_, int nx_, float x0_, float dt_, float dp_, float eps_)
/*< initialize radon2 operator >*/
{
	inv=inv_;
	
	nt=nt_;
	np=np_;
	nx=nx_;
	x0=x0_;
	dt=dt_;
	dp=dp_;
	eps=eps_;

	nfft=2*kiss_fft_next_fast_size(nt);
	nw=nfft/2+1;

	p=p_;
	xx=xx_;
	xxtmp=sf_floatalloc(nx);
	dd=sf_floatalloc2(nt, nx);
	mm=sf_floatalloc2(nt, np);
	cdd=(sf_complex*)malloc(nw*nx*sizeof(sf_complex));
	cmm=(sf_complex*)malloc(nw*np*sizeof(sf_complex));
	tmpr=(float*)fftwf_malloc(nfft*sizeof(float));
    	tmpc=(fftwf_complex*)fftwf_malloc(nw*sizeof(fftwf_complex));
    	fft1=fftwf_plan_dft_r2c_1d(nfft,tmpr,tmpc,FFTW_MEASURE);	
   	ifft1=fftwf_plan_dft_c2r_1d(nfft,tmpc,tmpr,FFTW_MEASURE);
}

void sf_radon2_set(bool par)
/*< set it to parabolic/linear radon transform >*/
{
	int ix;
	for (ix=0; ix < nx; ix++)/* normalize offsets */
	{
		if (par) xxtmp[ix] = (xx[ix]*xx[ix])/(x0*x0);
		else if (x0!=1.) xxtmp[ix]=xx[ix]/x0;
	}
}

void matrix_transpose(sf_complex *matrix, int nx, int nz)
/*< matrix transpose >*/
{
	int ix, iz;
	sf_complex *tmp=(sf_complex*)malloc(nx*nz*sizeof(sf_complex));
	if (tmp==NULL) {printf("out of memory!\n"); exit(1);}
	for(iz=0; iz<nz; iz++)
	for(ix=0; ix<nx; ix++)
		tmp[iz+nz*ix]=matrix[ix+nx*iz];

	memcpy(matrix, tmp, nx*nz*sizeof(sf_complex));
	free(tmp);
}

void sf_radon2_lop(bool adj, bool add, int nm, int nd, float *xxmm, float *yydd)
/*< radon operator: forward or inverse transform >*/
{
	int ix, ip, iw;
	float w;

    	if (nm != np*nt || nd != nx*nt) sf_error("%s: mismatched data sizes",__FILE__);
    	if(adj){
		memcpy(dd[0], yydd, nd*sizeof(float));
		if(!add) memset(mm[0],0,nm*sizeof(float));
	}else{
		memcpy(mm[0], xxmm, nm*sizeof(float));
		if(!add) memset(dd[0],0,nd*sizeof(float));
	}

	if(adj){// m(tau,p)=sum_{i=0}^{nx} d(t=tau+p*x_i,x_i)
		for(ix=0; ix<nx; ix++) // loop over offsets
		{
			memset(tmpr, 0, nfft*sizeof(float));
			memcpy(tmpr, dd[ix], nt*sizeof(float));
		 	fftwf_execute(fft1);// FFT: dd-->cdd
			memcpy(&cdd[ix*nw], tmpc, nw*sizeof(sf_complex));
		}
		matrix_transpose(cdd, nw, nx);
	}else{	// d(t,h)=sum_{i=0}^{np} m(tau=t-p_i*h,p_i)
		for(ip=0; ip<np; ip++) // loop over slopes
		{
			memset(tmpr, 0, nfft*sizeof(float));
			memcpy(tmpr, mm[ip], nt*sizeof(float));
		 	fftwf_execute(fft1);// FFT: mm-->cmm
			memcpy(&cmm[ip*nw], tmpc, nw*sizeof(float));			
		}
		matrix_transpose(cmm, nw, np);
	}

	myradon2_init(np, nx, dp, p, xxtmp);
	for(iw=0; iw<nw; iw++) 
	{
		w=2.*SF_PI*iw/(nfft*dt);
		myradon2_set(w);
		myradon2_lop(adj, false, np, nx, &cmm[iw*np], &cdd[iw*nx]);
		if(adj&&inv) myradon2_inv(&cmm[iw*np], &cmm[iw*np], eps);
	}


	if(adj){// m(tau,p)=sum_{i=0}^{nx} d(t=tau+p*x_i,x_i)
		matrix_transpose(cmm, np, nw);
		for(ip=0; ip<np; ip++) // loop over slopes
		{			
			memcpy(tmpc, &cmm[ip*nw], nw*sizeof(sf_complex));
		 	fftwf_execute(ifft1); // IFFT: cmm-->mm
			//for(iw=0; iw<nt; iw++) mm[ip][iw]=tmpr[iw]/nfft;
			for(iw=0; iw<nt; iw++) xxmm[iw+ip*nt]+=tmpr[iw]/nfft;
		}
	}else{// d(t,h)=sum_{i=0}^{np} m(tau=t-p_i*h,p_i)
		matrix_transpose(cdd, nx, nw);
		for(ix=0; ix<nx; ix++) // loop over offsets
		{
			memcpy(tmpc, &cdd[ix*nw], nw*sizeof(sf_complex));
		 	fftwf_execute(ifft1);// IFFT: cmm-->mm
			//for(iw=0; iw<nt; iw++) dd[ix][iw]=tmpr[iw]/nfft;
			for(iw=0; iw<nt; iw++) yydd[iw+ix*nt]+=tmpr[iw]/nfft;
		}
	}
}

void sf_radon2_close()
/*< free the allocated variables >*/
{
	free(xxtmp);
	free(*dd); free(dd);
	free(*mm); free(mm);
	free(cdd); 
	free(cmm); 
	fftwf_free(tmpr);
	fftwf_free(tmpc);
	fftwf_destroy_plan(fft1);
    	fftwf_destroy_plan(ifft1);
}

