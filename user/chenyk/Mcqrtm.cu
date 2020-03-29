/* Madagascar version of the 2-D QRTM (based on the original codes from Yufeng Wang)
*/
/*
  		Copyright (C) 
  		Yangkang Chen
		Nov 21, 2019  
		Zhejiang University
  
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
Wang, Y., H. Zhou, H. Chen, and Y. Chen, 2017, Adaptive stabilization for Q-compensated
reverse time migration, Geophysics, 83, S15-S32.
Wang, Y., H. Zhou, X. Zhao, Q. Zhang, P. Zhao, X. Yu, and Y. Chen, 2019, CuQ-RTM: A CUDA-based code package for stable and efficient Q-compensated RTM, Geophysics, 84, F1-F15.
Wang, Y., H. Zhou, X. Zhao, Q. Zhang, and Y. Chen, 2019, Q-compensated viscoelastic reverse time migration using mode-dependent adaptive stabilization scheme, Geophysics, 84, S301–S315.
*/

#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>
#include <string.h>
#include<time.h>
#include <cuda.h>

#include <cuda_runtime.h>
#define BLOCKSIZE 16
#define LTABLE 8
#define NTABLE 513
#define PI 3.1415926535898

extern "C" {
#include <rsf.h>
}

#include "mpi.h"



#define BLOCK_SIZE 16

#define OUTPUT_SNAP 0


//#include "headmulti.h"
#include "cufft.h"

//#define BATCH 834

//#include "Myfunctions.h"
//using namespace std;

extern "C"
struct Source
{
	int s_iz, s_ix, r_iz, *r_ix, r_n;
};

extern "C"
struct MultiGPU
{
	// cufft handle for forward and backward propagation
	cufftHandle PLAN_FORWARD;
	cufftHandle PLAN_BACKWARD;

	// host pagelock memory (variables needs to cudaMemcpyDeviceToHost)
	cufftComplex *u0, *u1, *u2;
	float *seismogram_obs;
	float *seismogram_dir;
	float *seismogram_syn;
	float *seismogram_rms;
	float *image_sources, *image_receivers;
	float *image_cor, *image_nor;

	// device global memory (variables needs to cudaMemcpyHostToDevice)
	int *d_r_ix;
	float *d_ricker;
	float *d_vp, *d_Gamma;
	cufftComplex *d_u0, *d_u1, *d_u2;
	cufftComplex *d_u0_inv, *d_u1_inv, *d_u2_inv;
	int *d_t_cp;
	float *d_u_cp;
	float *d_kx, *d_kz;
	float *d_kfilter, *d_kstabilization;	
	cufftComplex *d_uk, *d_uk0;
	cufftComplex *d_uk_inv, *d_uk0_inv;
	cufftComplex *d_Lap_uk, *d_amp_uk, *d_pha_uk;
	cufftComplex *d_Lap, *d_amp_Lap, *d_pha_Lap;
	float *d_seismogram;
	float *d_seismogram_rms;
	float *d_borders_up,*d_borders_bottom;
	float *d_borders_left,*d_borders_right;
	float *d_u2_final0, *d_u2_final1;
	float *d_image_sources,*d_image_receivers;
	float *d_image_cor, *d_image_nor;
};


// define multistream to prepare for streaming execution
struct Multistream
{
	cudaStream_t stream,stream_back;
};


//==========================================================
//  This subroutine is used for initializating wavefield variables
//  =========================================================
__global__ void cuda_kernel_initialization
(
	int ntx, int ntz, cufftComplex *u0, cufftComplex *u1, cufftComplex *u2, 
	cufftComplex *uk0, cufftComplex *uk, cufftComplex *Lap, cufftComplex *amp_Lap, cufftComplex *pha_Lap
)
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;
	int iz=by*BLOCK_SIZE+ty;
	int ix=bx*BLOCK_SIZE+tx;
	int ip=iz*ntx+ix;	

	if(iz>=0&&iz<=ntz-1&&ix>=0&&ix<=ntx-1)
	{
		u0[ip].x=0.0; u0[ip].y=0.0;
		u1[ip].x=0.0; u1[ip].y=0.0;
		u2[ip].x=0.0; u2[ip].y=0.0;
		uk0[ip].x=0.0; uk0[ip].y=0.0; 
		uk[ip].x=0.0; uk[ip].y=0.0; 
		Lap[ip].x=0.0; Lap[ip].y=0.0; 
		amp_Lap[ip].x=0.0; amp_Lap[ip].y=0.0; 
		pha_Lap[ip].x=0.0; pha_Lap[ip].y=0.0; 
	}
	__syncthreads();	
}


//==========================================================
//  This subroutine is used for updating wavefield variables
//  =========================================================
__global__ void cuda_kernel_update
(
	int ntx, int ntz, cufftComplex *u0, cufftComplex *u1, cufftComplex *u2
)
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;
	int iz=by*BLOCK_SIZE+ty;
	int ix=bx*BLOCK_SIZE+tx;
	int ip=iz*ntx+ix;	

	if(iz>=0&&iz<=ntz-1&&ix>=0&&ix<=ntx-1)
	{
		u0[ip].x=u1[ip].x;
		u0[ip].y=u1[ip].y;
		u1[ip].x=u2[ip].x;
		u1[ip].y=u2[ip].y;
	}
	__syncthreads();
}


//==========================================================
//  This subroutine is used for initializating image variables
//  =========================================================
__global__ void cuda_kernel_initialization_images
(
	int ntx, int ntz, float *image_cor, float *image_nor, float *image_sources, float *image_receivers
)
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;
	int iz=by*BLOCK_SIZE+ty;
	int ix=bx*BLOCK_SIZE+tx;
	int ip=iz*ntx+ix;	

	if(iz>=0&&iz<=ntz-1&&ix>=0&&ix<=ntx-1)
	{
		image_cor[ip]=0;
		image_nor[ip]=0;
		image_sources[ip]=0;
		image_receivers[ip]=0;
	}

	__syncthreads();
}


//=========================================================
//  This subroutine is used for adaptive stablization scheme 
//	For more details please refer to Eq(10) in our paper
//  =======================================================
__global__ void cuda_kernel_AdaSta
(
	int it, int ntx, int ntz, float dx, float dz, float dt, float *kstabilization,
	float *vp, float *Gamma, float avervp, float averGamma, float Omega0,
	float *kx, float *kz, cufftComplex *uk, float sigma, int Order
)
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;
	int iz=by*BLOCK_SIZE+ty;
	int ix=bx*BLOCK_SIZE+tx;
	int ip=iz*ntx+ix;

	float tau, ksi;
	if(iz>=0&&iz<=ntz-1&&ix>=0&&ix<=ntx-1)
	{
		tau = -powf(avervp,2*averGamma-1)*powf(Omega0,-2*averGamma)*sin(averGamma*PI);
		ksi = -tau*powf(avervp*cos(averGamma*PI/2),2)*powf((powf(kx[ix],2)+powf(kz[iz],2)), averGamma+0.5)/2;
		// adaptive stabilization coefficient
		kstabilization[ip] = (1+sigma*exp(Order*ksi*dt*it))/(1+sigma*exp(Order*ksi*dt*(it+1)));
		uk[ip].x *= kstabilization[ip]/(ntx*ntz);
		uk[ip].y *= kstabilization[ip]/(ntx*ntz);				
	}
	__syncthreads();
}


//=========================================================
//  This subroutine is used for low-pass filtering scheme 
//  =======================================================
__global__ void cuda_kernel_filter2d
(
	int ntx, int ntz, float dx, float dz, float *kfilter, float *kx, float *kz, cufftComplex *uk, 
	float kx_cut, float kz_cut, float taper_ratio
)
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;
	int iz=by*BLOCK_SIZE+ty;
	int ix=bx*BLOCK_SIZE+tx;
	int ip=iz*ntx+ix;

	float xc=kx_cut*ntx*dx;
	float zc=kz_cut*ntz*dz;
	float xs=(1-taper_ratio)*xc;
	float zs=(1-taper_ratio)*zc;
	int nxh=ntx/2;
	int nzh=ntz/2;	

	if(iz>=0&&iz<=ntz-1&&ix>=0&&ix<=ntx-1)
	{
		kfilter[ip] = 0;

		// filtering at x-direction
		if(ix>=0&&ix<xs&&iz>=0&&iz<ntz)
		{
			kfilter[ip]=1.0;
		}
		else if(ix>=xs&&ix<xc&&iz>=0&&iz<ntz)
		{
			kfilter[ip]=cos(PI/2.0*(ix-xs)/(xc-xs));		// cosin window
		}
		else if(ix>=xc&&ix<=nxh&&iz>=0&&iz<ntz)
		{
			kfilter[ip]=0.0;
		}
		else if(ix>=nxh&&ix<ntx-xc&&iz>=0&&iz<ntz)
		{
			kfilter[ip]=0.0;
		}
		else if(ix>=ntx-xc&&ix<ntx-xs&&iz>=0&&iz<ntz)
		{
			kfilter[ip]=sin(PI/2.0*(ix-(ntx-xc))/(xc-xs));
		}
		else if(ix>=ntx-xs&&ix<ntx&&iz>=0&&iz<ntz)
		{
			kfilter[ip]=1.0;
		}

		// filtering at z-direction
		if(iz>=0&&iz<zs&&ix>=0&&ix<ntx)
		{
			kfilter[ip]*=1.0;
		}
		else if(iz>=zs&&iz<zc&&ix>=0&&ix<ntx)
		{
			kfilter[ip]*=cos(PI/2.0*(iz-zs)/(zc-zs));		// cosin window
		}
		else if(iz>=zc&&iz<=nzh&&ix>=0&&ix<ntx)
		{
			kfilter[ip]*=0.0;
		}	
		else if(iz>=nzh&&iz<ntz-zc&&ix>=0&&ix<ntx)
		{
			kfilter[ip]*=0.0;
		}
		else if(iz>=ntz-zc&&iz<ntz-zs&&ix>=0&&ix<ntx)
		{
			kfilter[ip]*=sin(PI/2.0*(iz-(ntz-zc))/(zc-zs));
		}
		else if(iz>=ntz-zs&&iz<ntz&&ix>=0&&ix<ntx)
		{
			kfilter[ip]*=1.0;
		}
	}

	if(iz>=0&&iz<=ntz-1&&ix>=0&&ix<=ntx-1)
	{
		uk[ip].x *= kfilter[ip];
		uk[ip].y *= kfilter[ip];		
	}
	__syncthreads();
}


//==========================================================
//  This subroutine is used for defining k
// =========================================================
__global__ void cuda_kernel_k_define
(
	int ntx, int ntz, float dx, float dz, float *kx, float *kz
)
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;
	int iz=by*BLOCK_SIZE+ty;
	int ix=bx*BLOCK_SIZE+tx;

	int nxh=ntx/2;
	int nzh=ntz/2;
	float dkx=1.0/(ntx*dx);
	float dkz=1.0/(ntz*dz);

	if(ix>=0&&ix<=nxh)
	{
		kx[ix]=2*PI*ix*dkx;
	}
	if(ix>=nxh&&ix<ntx)
	{
		kx[ix]=kx[ntx-ix];
	}
	if(iz>=0&&iz<=nzh)
	{
		kz[iz]=2*PI*iz*dkz;
	}
	if(iz>=nzh&&iz<ntz)
	{
		kz[iz]=kz[ntz-iz];
	}
	__syncthreads();
}


//==========================================================
//  This subroutine is used for calculating forward wavefileds in k-space
//  ========================================================
__global__ void cuda_kernel_visco_PSM_2d_forward_k_space
(
	float beta1, float beta2,
	int it, int nt, int ntx, int ntz, float dx, float dz, float dt, 
	float *vp, float *Gamma, float averGamma, float f0, float Omega0, 
	float *kx, float *kz, 
	cufftComplex *uk, cufftComplex *uk0, 
	cufftComplex *Lap_uk, cufftComplex *amp_uk, cufftComplex *pha_uk
)
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;
	int iz=by*BLOCK_SIZE+ty;
	int ix=bx*BLOCK_SIZE+tx;
	int ip=iz*ntx+ix;

	if(iz>=0&&iz<=ntz-1&&ix>=0&&ix<=ntx-1)
	{
		Lap_uk[ip].x=-(powf(kx[ix],2)+powf(kz[iz],2))*uk[ip].x;
		Lap_uk[ip].y=-(powf(kx[ix],2)+powf(kz[iz],2))*uk[ip].y;
		if(beta1!=0)
		{
			pha_uk[ip].x=powf((powf(kx[ix],2)+powf(kz[iz],2)), averGamma+1)*uk[ip].x;
			pha_uk[ip].y=powf((powf(kx[ix],2)+powf(kz[iz],2)), averGamma+1)*uk[ip].y;
		}
		if(beta2!=0)
		{
			amp_uk[ip].x=powf((powf(kx[ix],2)+powf(kz[iz],2)), averGamma+0.5)*(uk[ip].x-uk0[ip].x)/dt;
			amp_uk[ip].y=powf((powf(kx[ix],2)+powf(kz[iz],2)), averGamma+0.5)*(uk[ip].y-uk0[ip].y)/dt;
		}					
		uk0[ip].x=uk[ip].x;
		uk0[ip].y=uk[ip].y;
	}
	__syncthreads();
}


//==========================================================
//  This subroutine is used for calculating forward wavefileds in x-space
//  ========================================================
__global__ void cuda_kernel_visco_PSM_2d_forward_x_space
(
	float beta1, float beta2,
	int it, int nt, int ntx, int ntz, int nx, int nz, int L, float dx, float dz, float dt, 
	float *vp, float *Gamma, float averGamma, float f0, float Omega0, 
	float *seismogram, int *r_ix, int r_iz, int rnmax, float *ricker, int s_ix, int s_iz,
	cufftComplex *u0, cufftComplex *u1, cufftComplex *u2,
	cufftComplex *Lap, cufftComplex *amp_Lap, cufftComplex *pha_Lap,
	float *borders_up, float *borders_bottom, float *borders_left, float *borders_right,
	float *u2_final0, float *u2_final1,
	int Sto_Rec, int vp_type
)
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;
	int iz=by*BLOCK_SIZE+ty;
	int ix=bx*BLOCK_SIZE+tx;
	int ip=iz*ntx+ix;

	int icp;
	float eta, tau;

	if(iz>=0&&iz<=ntz-1&&ix>=0&&ix<=ntx-1)
	{		
		eta= -powf(vp[ip],2*Gamma[ip])*powf(Omega0,-2*Gamma[ip])*cos(Gamma[ip]*PI);
		tau= -powf(vp[ip],2*Gamma[ip]-1)*powf(Omega0,-2*Gamma[ip])*sin(Gamma[ip]*PI);

		// scale fft by dividing (ntx*ntz)
		Lap[ip].x=Lap[ip].x/(ntx*ntz);
		pha_Lap[ip].x=pha_Lap[ip].x/(ntx*ntz);
		amp_Lap[ip].x=amp_Lap[ip].x/(ntx*ntz);
		u2[ip].x=powf(vp[ip]*cos(Gamma[ip]*PI/2),2)*powf(dt,2)
			*(
				Lap[ip].x
				+beta1*(eta*pha_Lap[ip].x-Lap[ip].x)
				+beta2*tau*amp_Lap[ip].x
			)
			+2*u1[ip].x-u0[ip].x;
	}
	// add Ricker source
	if(iz==s_iz&&ix==s_ix)
	{
		u2[ip].x+=ricker[it];
	}	
	// record Seismogram
	if(ix>=r_ix[0]&&ix<=r_ix[rnmax-1]&&iz==r_iz)
		seismogram[it*rnmax+ix-r_ix[0]]=u2[ip].x;

	// store borders and final two-step wavefileds for wavefield reconstruction
	if(Sto_Rec==0&&vp_type==2)
	{
		if(ix>=L&&ix<=ntx-L-1&&iz==L)
		{
			borders_up[it*nx+ix-L]=u2[ip].x;
		}
		if(ix>=L&&ix<=ntx-L-1&&iz==ntz-L-1)
		{
			borders_bottom[it*nx+ix-L]=u2[ip].x;
		}
		if(iz>=L&&iz<=ntz-L-1&&ix==L)
		{
			borders_left[it*nz+iz-L]=u2[ip].x;
		}
		if(iz>=L&&iz<=ntz-L-1&&ix==ntx-L-1)
		{
			borders_right[it*nz+iz-L]=u2[ip].x;
		}
		if(it==nt-1)
		{
			if(iz>=0&&iz<=ntz-1&&ix>=0&&ix<=ntx-1)
			{				
				u2_final0[ip]=u2[ip].x;
				u2_final1[ip]=u1[ip].x;				
			}			
		}
	}
	__syncthreads();
}


//==========================================================
//  This subroutine is used for writing checkpoints
//  ========================================================
__global__ void cuda_kernel_checkpoints_Out
(
	int it, int nt, int ntx, int ntz, int nx, int nz, int L, float dx, float dz, float dt, 
	cufftComplex *u1, cufftComplex *u2,
	float *u_cp, int N_cp, int *t_cp
)
{

	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;
	int iz=by*BLOCK_SIZE+ty;
	int ix=bx*BLOCK_SIZE+tx;
	int ip=iz*ntx+ix;

	int icp;

	for(icp=0;icp<N_cp;icp++)
	{
		if(icp%2==1&&it==t_cp[icp])
		{
			if(iz>=0&&iz<=ntz-1&&ix>=0&&ix<=ntx-1)
			{
				u_cp[icp*ntx*ntz+ip]=u2[ip].x;
				u_cp[(icp-1)*ntx*ntz+ip]=u1[ip].x;
			}
		}
	}
	__syncthreads();
}


//==========================================================
//  This two subroutines are used for initializing Final two wavefileds
//  =========================================================
__global__ void cuda_kernel_initialization_Finals
(
	int ntx, int ntz, cufftComplex *u0, cufftComplex *u1, float *u2_final0, float *u2_final1
)
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;
	int iz=by*BLOCK_SIZE+ty;
	int ix=bx*BLOCK_SIZE+tx;
	int ip=iz*ntx+ix;	
	if(iz>=0&&iz<=ntz-1&&ix>=0&&ix<=ntx-1)
	{
		u0[ip].x=u2_final0[ip];
		u1[ip].x=u2_final1[ip];
	}
	__syncthreads();	
}


/*==========================================================
  This subroutine is used for calculating reconstructed wavefileds in k-space
  ===========================================================*/

__global__ void cuda_kernel_visco_PSM_2d_reconstruction_k_space
(
	float beta1, float beta2,
	int it, int nt, int ntx, int ntz, int nx, int nz, int L, float dx, float dz, float dt, 
	float *vp, float *Gamma, float averGamma, float f0, float Omega0, 
	float *kx, float *kz, 
	cufftComplex *uk, cufftComplex *uk0, 
	cufftComplex *Lap_uk, cufftComplex *amp_uk, cufftComplex *pha_uk
)
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;
	int iz=by*BLOCK_SIZE+ty;
	int ix=bx*BLOCK_SIZE+tx;
	int ip=iz*ntx+ix;

	if(iz>=0&&iz<=ntz-1&&ix>=0&&ix<=ntx-1)
	{
		Lap_uk[ip].x=-(powf(kx[ix],2)+powf(kz[iz],2))*uk[ip].x;
		Lap_uk[ip].y=-(powf(kx[ix],2)+powf(kz[iz],2))*uk[ip].y;
		if(beta1!=0)
		{
			pha_uk[ip].x=powf((powf(kx[ix],2)+powf(kz[iz],2)), averGamma+1)*uk[ip].x;
			pha_uk[ip].y=powf((powf(kx[ix],2)+powf(kz[iz],2)), averGamma+1)*uk[ip].y;
		}
		if(beta2!=0)
		{
			amp_uk[ip].x=powf((powf(kx[ix],2)+powf(kz[iz],2)), averGamma+0.5)*(uk[ip].x-uk0[ip].x)/dt;
			amp_uk[ip].y=powf((powf(kx[ix],2)+powf(kz[iz],2)), averGamma+0.5)*(uk[ip].y-uk0[ip].y)/dt;
		}					
		uk0[ip].x=uk[ip].x;
		uk0[ip].y=uk[ip].y;			
	}
	__syncthreads();
}


//==========================================================
//  This subroutine is used for calculating reconstructed wavefileds in x-space
//  =========================================================
__global__ void cuda_kernel_visco_PSM_2d_reconstruction_x_space
(
	float beta1, float beta2,
	int it, int nt, int ntx, int ntz, int nx, int nz, int L, float dx, float dz, float dt, 
	float *vp, float *Gamma, float averGamma, float f0, float Omega0, 
	float *ricker, int s_ix, int s_iz,
	cufftComplex *u0, cufftComplex *u1, cufftComplex *u2,
	cufftComplex *Lap, cufftComplex *amp_Lap, cufftComplex *pha_Lap,
	float *borders_up, float *borders_bottom, float *borders_left, float *borders_right
)
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;
	int iz=by*BLOCK_SIZE+ty;
	int ix=bx*BLOCK_SIZE+tx;
	int ip=iz*ntx+ix;

	int icp;
	float eta, tau;

	if(iz>=0&&iz<=ntz-1&&ix>=0&&ix<=ntx-1)
	{		
		eta= -powf(vp[ip],2*Gamma[ip])*powf(Omega0,-2*Gamma[ip])*cos(Gamma[ip]*PI);
		tau= -powf(vp[ip],2*Gamma[ip]-1)*powf(Omega0,-2*Gamma[ip])*sin(Gamma[ip]*PI);

		// scale fft by dividing (ntx*ntz)
		Lap[ip].x=Lap[ip].x/(ntx*ntz);
		pha_Lap[ip].x=pha_Lap[ip].x/(ntx*ntz);
		amp_Lap[ip].x=amp_Lap[ip].x/(ntx*ntz);

		u2[ip].x=powf(vp[ip]*cos(Gamma[ip]*PI/2),2)*powf(dt,2)
			*(
				Lap[ip].x
				+beta1*(eta*pha_Lap[ip].x-Lap[ip].x)
				+beta2*tau*amp_Lap[ip].x
			)
			+2*u1[ip].x-u0[ip].x;
	}

	// add borders 
	if(ix>=L&&ix<=ntx-L-1&&iz==L)
	{
		u2[ip].x=borders_up[it*nx+ix-L];
	}
	if(ix>=L&&ix<=ntx-L-1&&iz==ntz-L-1)
	{
		u2[ip].x=borders_bottom[it*nx+ix-L];
	}
	if(iz>=L&&iz<=ntz-L-1&&ix==L)
	{
		u2[ip].x=borders_left[it*nz+iz-L];
	}
	if(iz>=L&&iz<=ntz-L-1&&ix==ntx-L-1)
	{
		u2[ip].x=borders_right[it*nz+iz-L];
	}
	__syncthreads();
}


//==========================================================
//  This subroutine is used for reading checkpoints
//  =========================================================
__global__ void cuda_kernel_checkpoints_In
(
	int it, int nt, int ntx, int ntz, int nx, int nz, int L, float dx, float dz, float dt, 
	cufftComplex *u1, cufftComplex *u2,
	float *u_cp, int N_cp, int *t_cp
)
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;
	int iz=by*BLOCK_SIZE+ty;
	int ix=bx*BLOCK_SIZE+tx;
	int ip=iz*ntx+ix;

	int icp;

	for(icp=0;icp<N_cp;icp++)
	{
		if(icp%2==0&&it==t_cp[icp])
		{
			if(iz>=0&&iz<=ntz-1&&ix>=0&&ix<=ntx-1)
			{
				u2[ip].x=u_cp[icp*ntx*ntz+ip];
				u1[ip].x=u_cp[(icp+1)*ntx*ntz+ip];
			}
		}
	}
	__syncthreads();
}


//==========================================================
//  This subroutine is used for calculating backward wavefileds in k-space
//  =========================================================
__global__ void cuda_kernel_visco_PSM_2d_backward_k_space
(
	float beta1, float beta2,
	int it, int nt, int ntx, int ntz, float dx, float dz, float dt, 
	float *vp, float *Gamma, float averGamma, float f0, float Omega0, 
	float *kx, float *kz, 
	cufftComplex *uk, cufftComplex *uk0, 
	cufftComplex *Lap_uk, cufftComplex *amp_uk, cufftComplex *pha_uk
)
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;
	int iz=by*BLOCK_SIZE+ty;
	int ix=bx*BLOCK_SIZE+tx;
	int ip=iz*ntx+ix;

	if(iz>=0&&iz<=ntz-1&&ix>=0&&ix<=ntx-1)
	{
		Lap_uk[ip].x=-(powf(kx[ix],2)+powf(kz[iz],2))*uk[ip].x;
		Lap_uk[ip].y=-(powf(kx[ix],2)+powf(kz[iz],2))*uk[ip].y;
		if(beta1!=0)
		{
			pha_uk[ip].x=powf((powf(kx[ix],2)+powf(kz[iz],2)), averGamma+1)*uk[ip].x;
			pha_uk[ip].y=powf((powf(kx[ix],2)+powf(kz[iz],2)), averGamma+1)*uk[ip].y;
		}
		if(beta2!=0)
		{
			amp_uk[ip].x=powf((powf(kx[ix],2)+powf(kz[iz],2)), averGamma+0.5)*(uk[ip].x-uk0[ip].x)/dt;
			amp_uk[ip].y=powf((powf(kx[ix],2)+powf(kz[iz],2)), averGamma+0.5)*(uk[ip].y-uk0[ip].y)/dt;
		}					
		uk0[ip].x=uk[ip].x;
		uk0[ip].y=uk[ip].y;			
	}
	__syncthreads();
}


//==========================================================
//  This subroutine is used for calculating backward wavefileds in x-space
//  ========================================================
__global__ void cuda_kernel_visco_PSM_2d_backward_x_space
(
	float beta1, float beta2,
	int it, int nt, int ntx, int ntz, float dx, float dz, float dt, 
	float *vp, float *Gamma, float averGamma, float f0, float Omega0, 
	float *seismogram_rms, int *r_ix, int r_iz, int s_ix, int rnmax, int nrx_obs,
	cufftComplex *u0, cufftComplex *u1, cufftComplex *u2,
	cufftComplex *Lap, cufftComplex *amp_Lap, cufftComplex *pha_Lap
)
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;
	int iz=by*BLOCK_SIZE+ty;
	int ix=bx*BLOCK_SIZE+tx;
	int ip=iz*ntx+ix;

	float eta, tau;

	if(iz>=0&&iz<=ntz-1&&ix>=0&&ix<=ntx-1)
	{
		eta= -powf(vp[ip],2*Gamma[ip])*powf(Omega0,-2*Gamma[ip])*cos(Gamma[ip]*PI);
		tau= -powf(vp[ip],2*Gamma[ip]-1)*powf(Omega0,-2*Gamma[ip])*sin(Gamma[ip]*PI);

		// scaling fft 
		Lap[ip].x=Lap[ip].x/(ntx*ntz);
		pha_Lap[ip].x=pha_Lap[ip].x/(ntx*ntz);
		amp_Lap[ip].x=amp_Lap[ip].x/(ntx*ntz);

		u2[ip].x=powf(vp[ip]*cos(Gamma[ip]*PI/2),2)*powf(dt,2)
			*(
				Lap[ip].x
				+beta1*(eta*pha_Lap[ip].x-Lap[ip].x)
				+beta2*tau*amp_Lap[ip].x
			)
			+2*u1[ip].x-u0[ip].x;	
	}

	// add seismogram as source
	int irx_min = s_ix-nrx_obs;
	int irx_max = s_ix+nrx_obs;
	if(irx_min<r_ix[0])
		irx_min = r_ix[0];
	if(irx_max>r_ix[rnmax-1])
		irx_max = r_ix[rnmax-1];
	if(ix>=irx_min&&ix<=irx_max&&iz==r_iz)
		u2[ip].x=seismogram_rms[it*rnmax+ix-r_ix[0]];
	__syncthreads();
}


//==========================================================
//  This subroutine is used for imaging
// ========================================================
__global__ void cuda_kernel_image
(
	int ntx, int ntz, int L,
	cufftComplex *u2_inv, cufftComplex *u2,
	float *image_cor, float *image_sources, float *image_receivers
)
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;
	int iz=by*BLOCK_SIZE+ty;
	int ix=bx*BLOCK_SIZE+tx;
	int ip=iz*ntx+ix;

	if(iz>=L&&iz<=ntz-L-1&&ix>=L&&ix<=ntx-L-1)
	{
		image_cor[ip]+=u2_inv[ip].x*u2[ip].x;
		image_sources[ip]+=u2_inv[ip].x*u2_inv[ip].x;
		image_receivers[ip]+=u2[ip].x*u2[ip].x;    
	}
	__syncthreads();
}

//==========================================================
//  This subroutine is used for absorbing boundary condition
//  ========================================================
__global__ void cuda_kernel_MTF_2nd
(
	int L, int ntx, int ntz, float dx, float dz, float dt, 
	float *vp, cufftComplex *u0, cufftComplex *u1, cufftComplex *u2
)
{
	int bx=blockIdx.x;
	int by=blockIdx.y;
	int tx=threadIdx.x;
	int ty=threadIdx.y;
	int iz=by*BLOCK_SIZE+ty;
	int ix=bx*BLOCK_SIZE+tx;
	int ip=iz*ntx+ix;

	int ipp=iz*ntx+(ntx-1-ix);
	int ippp=(ntz-1-iz)*ntx+ix;

	float alpha=1.0;
	float w, s, t1, t2, t3;

	// left ABC ...
	if(ix>=0&&ix<=L-1&&iz>=0&&iz<=ntz-1)
	{
		w=1-1.0*ix/L;
		s=alpha*vp[ip]*dt/dx;
		t1=(2-s)*(1-s)/2;
		t2=s*(2-s);
		t3=s*(s-1)/2;	

		u2[ip].x=w*
			(
				(1*2)*
				(
					t1*u1[ip].x+t2*u1[ip+1].x+t3*u1[ip+2].x
				)
				+(-1*1)*
				(
					t1*t1*u0[ip].x
					+2*t1*t2*u0[ip+1].x
					+(2*t1*t3+t2*t2)*u0[ip+2].x
					+2*t2*t3*u0[ip+3].x
					+t3*t3*u0[ip+4].x
				)
			)
			+(1-w)*u2[ip].x;								
	}

	// right ABC ...
	if(ix>=ntx-L&&ix<=ntx-1&&iz>=0&&iz<=ntz-1)
	{
		w=1-1.0*(ntx-1-ix)/L;

		s=alpha*vp[ip]*dt/dx;
		t1=(2-s)*(1-s)/2;
		t2=s*(2-s);
		t3=s*(s-1)/2;			

		u2[ip].x=w*
				(
					(1*2)*
					(
						t1*u1[ip].x
						+t2*u1[ip-1].x
						+t3*u1[ip-2].x
					)
					+(-1*1)*
					(
						t1*t1*u0[ip].x
						+2*t1*t2*u0[ip-1].x
						+(2*t1*t3+t2*t2)*u0[ip-2].x
						+2*t2*t3*u0[ip-3].x
						+t3*t3*u0[ip-4].x
					)
				)
				+(1-w)*u2[ip].x;								
	}


	// up ABC ...
	if(iz>=0&&iz<=L-1&&ix>=0&&ix<=ntx-1)
	{
		w=1-1.0*iz/L;
		s=alpha*vp[ip]*dt/dz;	
		t1=(2-s)*(1-s)/2;
		t2=s*(2-s);
		t3=s*(s-1)/2;

		u2[ip].x=w*
			(
				(1*2)*
				(
					t1*u1[ip].x
					+t2*u1[ip+ntx].x
					+t3*u1[ip+2*ntx].x
				)
				+(-1*1)*
				(
					t1*t1*u0[ip].x
					+2*t1*t2*u0[ip+ntx].x
					+(2*t1*t3+t2*t2)*u0[ip+2*ntx].x
					+2*t2*t3*u0[ip+3*ntx].x
					+t3*t3*u0[ip+4*ntx].x
				)
			)
			+(1-w)*u2[ip].x;			
	}

	// bottom ABC ...
	if(iz>=ntz-L&&iz<=ntz-1&&ix>=0&&ix<=ntx-1)
	{
		w=1-1.0*(ntz-1-iz)/L;
		s=alpha*vp[ip]*dt/dz;	
		t1=(2-s)*(1-s)/2;
		t2=s*(2-s);
		t3=s*(s-1)/2;

		u2[ip].x=w*
			(
				(1*2)*
				(
					t1*u1[ip].x
					+t2*u1[ip-ntx].x
					+t3*u1[ip-2*ntx].x
				)
				+(-1*1)*
				(
					t1*t1*u0[ip].x
					+2*t1*t2*u0[ip-ntx].x
					+(2*t1*t3+t2*t2)*u0[ip-2*ntx].x
					+2*t2*t3*u0[ip-3*ntx].x
					+t3*t3*u0[ip-4*ntx].x
				)
			)
			+(1-w)*u2[ip].x;		
	}
	__syncthreads();	
}







//==========================================================
//  This subroutine are used for forward modeling
//	For more details please refer to Eq(1) in our paper
//  =========================================================
extern "C"
void cuda_visco_PSM_2d_forward
(
	int beta1, int beta2,
	int nt, int ntx, int ntz, int ntp, int nx, int nz, int L, float dx, float dz, float dt,
	float *vp, float *Gamma, float avervp, float averGamma, float f0, float Omega0, float *ricker,
	int myid, int is, struct Source ss[], struct MultiGPU plan[], int GPU_N, int rnmax, int nrx_obs, int N_cp, int *t_cp,
	float kx_cut, float kz_cut, float sigma, int Order,	float taper_ratio, float *kfilter, float *kstabilization,
	int Sto_Rec, int vp_type, int Save_Not, int Filtertype
)
{
	int i, it, ix, iz, ip, icp;
	size_t size_model=sizeof(float)*ntp;
	FILE *fp;
	char filename[40];
	float *u2_real;
	u2_real = (float*)malloc(sizeof(float)*ntp);

	// define multistream  variable
	Multistream plans[GPU_N];

	// define streaming cufft handle (very important!!!)
	for(i=0;i<GPU_N;i++)
	{
		cudaSetDevice(i);
		cudaStreamCreate(&plans[i].stream);	
		cufftSetStream(plan[i].PLAN_FORWARD,plans[i].stream);
		cufftSetStream(plan[i].PLAN_BACKWARD,plans[i].stream);
	}	

	// block size 16*16; 
	dim3 dimBlock(BLOCK_SIZE,BLOCK_SIZE);
	// grid size ntx/16*ntz/16
	dim3 dimGrid((ntx+dimBlock.x-1)/dimBlock.x,(ntz+dimBlock.y-1)/dimBlock.y);

	// copy the vectors from the host to the device
	for(i=0;i<GPU_N;i++)
	{
		cudaSetDevice(i);
		cudaMemcpyAsync(plan[i].d_r_ix,ss[is+i].r_ix,sizeof(float)*rnmax,cudaMemcpyHostToDevice,plans[i].stream);
		cudaMemcpyAsync(plan[i].d_ricker,ricker,sizeof(float)*nt,cudaMemcpyHostToDevice,plans[i].stream);
		cudaMemcpyAsync(plan[i].d_vp,vp,size_model,cudaMemcpyHostToDevice,plans[i].stream);
		cudaMemcpyAsync(plan[i].d_Gamma,Gamma,size_model,cudaMemcpyHostToDevice,plans[i].stream);
		cudaMemcpyAsync(plan[i].d_t_cp,t_cp,N_cp*sizeof(int),cudaMemcpyHostToDevice,plans[i].stream);
	}

	// initializing wavefield variables and define k variables
	for(i=0;i<GPU_N;i++)
	{
		cudaSetDevice(i);
		cuda_kernel_initialization<<<dimGrid,dimBlock,0,plans[i].stream>>>
			(ntx, ntz, plan[i].d_u0, plan[i].d_u1, plan[i].d_u2, plan[i].d_uk0, plan[i].d_uk, 
			plan[i].d_Lap, plan[i].d_amp_Lap, plan[i].d_pha_Lap);
		cuda_kernel_k_define<<<dimGrid,dimBlock,0,plans[i].stream>>>
			(ntx, ntz, dx, dz, plan[i].d_kx, plan[i].d_kz);
	}

	// forward time iteration
	for(it=0;it<nt;it++)  
	{
		for(i=0;i<GPU_N;i++)
		{
			cudaSetDevice(i);
			cufftExecC2C(plan[i].PLAN_FORWARD,plan[i].d_u1,plan[i].d_uk,CUFFT_FORWARD); //CUFFT_FORWARD

			cuda_kernel_visco_PSM_2d_forward_k_space<<<dimGrid,dimBlock,0,plans[i].stream>>>
				(
					beta1, beta2,
					it, nt, ntx, ntz, dx, dz, dt, 
					plan[i].d_vp, plan[i].d_Gamma, averGamma, f0, Omega0,
					plan[i].d_kx, plan[i].d_kz,  
					plan[i].d_uk, plan[i].d_uk0, 
					plan[i].d_Lap_uk, plan[i].d_amp_uk, plan[i].d_pha_uk
				);		

			cufftExecC2C(plan[i].PLAN_BACKWARD,plan[i].d_Lap_uk,plan[i].d_Lap,CUFFT_INVERSE); //CUFFT_INVERSE

			if(beta1!=0)
			{	
				cufftExecC2C(plan[i].PLAN_BACKWARD,plan[i].d_pha_uk,plan[i].d_pha_Lap,CUFFT_INVERSE); //CUFFT_INVERSE					
			}

			if(beta2!=0)
			{
				if(beta2<0&&Filtertype==0)		// low-pass filtering
				{
					cuda_kernel_filter2d<<<dimGrid,dimBlock,0,plans[i].stream>>>
						(ntx, ntz, dx, dz, plan[i].d_kfilter, plan[i].d_kx, plan[i].d_kz, plan[i].d_amp_uk, kx_cut, kz_cut, taper_ratio);
				}
							
				cufftExecC2C(plan[i].PLAN_BACKWARD,plan[i].d_amp_uk,plan[i].d_amp_Lap,CUFFT_INVERSE); //CUFFT_INVERSE
			}

			cuda_kernel_visco_PSM_2d_forward_x_space<<<dimGrid,dimBlock,0,plans[i].stream>>>
				(
					beta1, beta2,
					it, nt, ntx, ntz, nx, nz, L, dx, dz, dt, 
					plan[i].d_vp, plan[i].d_Gamma, averGamma, f0, Omega0,
					plan[i].d_seismogram, plan[i].d_r_ix, ss[is+i].r_iz, rnmax, plan[i].d_ricker, ss[is+i].s_ix, ss[is+i].s_iz,
					plan[i].d_u0, plan[i].d_u1, plan[i].d_u2,
					plan[i].d_Lap, plan[i].d_amp_Lap, plan[i].d_pha_Lap,
					plan[i].d_borders_up, plan[i].d_borders_bottom, plan[i].d_borders_left, plan[i].d_borders_right,
					plan[i].d_u2_final0, plan[i].d_u2_final1,
					Sto_Rec, vp_type				
				);
		
			if(beta2<0&&Filtertype==1)			// adaptive stabilization
			{
				cufftExecC2C(plan[i].PLAN_FORWARD,plan[i].d_u2,plan[i].d_uk,CUFFT_FORWARD); //CUFFT_FORWARD
				cuda_kernel_AdaSta<<<dimGrid,dimBlock,0,plans[i].stream>>>
				(
					it, ntx, ntz, dx, dz, dt, plan[i].d_kstabilization,
					plan[i].d_vp, plan[i].d_Gamma, avervp, averGamma, Omega0,
					plan[i].d_kx, plan[i].d_kz, plan[i].d_uk, sigma, Order
				);	
				cufftExecC2C(plan[i].PLAN_BACKWARD,plan[i].d_uk,plan[i].d_u2,CUFFT_INVERSE); //CUFFT_INVERSE	
			}
		
			// MTF absorbing boundary condition
			cuda_kernel_MTF_2nd<<<dimGrid,dimBlock,0,plans[i].stream>>>
				(L, ntx, ntz, dx, dz, dt, plan[i].d_vp, plan[i].d_u0, plan[i].d_u1, plan[i].d_u2);

			// record wavefields at checkpoints
			if(Sto_Rec==0&&vp_type==2)
			{
				cuda_kernel_checkpoints_Out<<<dimGrid,dimBlock,0,plans[i].stream>>>
					(
						it, nt, ntx, ntz, nx, nz, L, dx, dz, dt, 
						plan[i].d_u1, plan[i].d_u2,
						plan[i].d_u_cp, N_cp, plan[i].d_t_cp			
					);
			}

			// write wavefields at checkpoints and last two time steps
			if(Sto_Rec==1&&vp_type==2||Save_Not==1)
			{
				cudaMemcpyAsync(plan[i].u2,plan[i].d_u2,sizeof(cufftComplex)*ntp,cudaMemcpyDeviceToHost,plans[i].stream);

				sprintf(filename,"./output/GPU_%d_u2_%d.dat",i,it);     
				fp=fopen(filename,"wb");
				for(ix=0;ix<ntx-0;ix++)
				{
					for(iz=0;iz<ntz-0;iz++)
					{
						u2_real[iz*ntx+ix]=plan[i].u2[iz*ntx+ix].x;			
						fwrite(&u2_real[iz*ntx+ix],sizeof(float),1,fp);
					}
				}
				fclose(fp);			
			}


			// write stabilized k-space wavefiled at timestep of 800 (for test our stabilization performance)
			if(it==500)
			{
				cudaMemcpyAsync(kstabilization,plan[i].d_kstabilization,sizeof(float)*ntp,cudaMemcpyDeviceToHost,plans[i].stream);
				cudaMemcpyAsync(kfilter,plan[i].d_kfilter,sizeof(float)*ntp,cudaMemcpyDeviceToHost,plans[i].stream);
				sprintf(filename,"./output/kstabilization.dat");     
				fp=fopen(filename,"wb");
				for(ix=0;ix<ntx-0;ix++)
				{
					for(iz=0;iz<ntz-0;iz++)
					{	
						fwrite(&kstabilization[iz*ntx+ix],sizeof(float),1,fp);
					}
				}
				fclose(fp);					
				sprintf(filename,"./output/kfilter.dat");     
				fp=fopen(filename,"wb");
				for(ix=0;ix<ntx-0;ix++)
				{
					for(iz=0;iz<ntz-0;iz++)
					{	
						fwrite(&kfilter[iz*ntx+ix],sizeof(float),1,fp);
					}
				}
				fclose(fp);	
			}

			// updating wavefields
			cuda_kernel_update<<<dimGrid,dimBlock,0,plans[i].stream>>>
				(ntx, ntz, plan[i].d_u0, plan[i].d_u1, plan[i].d_u2);

			/*if(myid==0&&it%100==0)
			{
				sf_warning("shot %d forward %d has finished!", is+i+1, it);
			}*/

		}// GPU_N end	
	}// nt end

	sf_warning("shot %d forward has finished!", is+i);

	for(i=0;i<GPU_N;i++)
	{
		cudaSetDevice(i);

		// copy seismograms to host memory
		if(vp_type==0)	// homogeneous model
		{
			cudaMemcpyAsync(plan[i].seismogram_dir,plan[i].d_seismogram,
					sizeof(float)*ss[is+i].r_n*nt,cudaMemcpyDeviceToHost,plans[i].stream);
		}
		else if(vp_type==1)	// ture model
		{
			cudaMemcpyAsync(plan[i].seismogram_obs,plan[i].d_seismogram,
					sizeof(float)*ss[is+i].r_n*nt,cudaMemcpyDeviceToHost,plans[i].stream);
		}
		else if(vp_type==2)	// initial model
		{
			cudaMemcpyAsync(plan[i].seismogram_syn,plan[i].d_seismogram,
					sizeof(float)*ss[is+i].r_n*nt,cudaMemcpyDeviceToHost,plans[i].stream);
		}
	}

	for(i=0;i<GPU_N;i++)
	{
		cudaSetDevice(i);
		cudaDeviceSynchronize();
	}

	//free the memory of DEVICE
	for(i=0;i<GPU_N;i++)
	{
		cudaSetDevice(i);	
		cudaStreamDestroy(plans[i].stream);
	}
	free(u2_real);
}



//==========================================================
//  This subroutine are used for backward modeling
//	For more details please refer to Eq(5) in our paper
//  =========================================================
extern "C"
void cuda_visco_PSM_2d_backward
(
	int beta1, int beta2,
	int nt, int ntx, int ntz, int ntp, int nx, int nz, int L, float dx, float dz, float dt,
	float *vp, float *Gamma, float avervp, float averGamma, float f0, float Omega0, float *ricker,
	int myid, int is, struct Source ss[], struct MultiGPU plan[], int GPU_N, int rnmax, int nrx_obs, int N_cp, int *t_cp,
	float kx_cut, float kz_cut, float sigma, int Order,	float taper_ratio, float *kfilter, float *kstabilization,
	int Sto_Rec, int Save_Not, int Filtertype
)
{
	int i, it, ix, iz, ip;
	size_t size_model=sizeof(float)*ntp;
	FILE *fp;
	char filename[40];
	float *u2_real;
	u2_real = (float*)malloc(sizeof(float)*ntp);

	// define multistream  variable
	Multistream plans[GPU_N];

	// define streaming cufft handle (very important!!!)
	for(i=0;i<GPU_N;i++)
	{
		cudaSetDevice(i);
		cudaStreamCreate(&plans[i].stream);	
		cufftSetStream(plan[i].PLAN_FORWARD,plans[i].stream);
		cufftSetStream(plan[i].PLAN_BACKWARD,plans[i].stream);
	}	

	dim3 dimBlock(BLOCK_SIZE,BLOCK_SIZE);
	dim3 dimGrid((ntx+dimBlock.x-1)/dimBlock.x,(ntz+dimBlock.y-1)/dimBlock.y);


	for(i=0;i<GPU_N;i++)
	{
		cudaSetDevice(i);

		// initializating wavefields for reconstruction
		cuda_kernel_initialization<<<dimGrid,dimBlock,0,plans[i].stream>>>
			(ntx, ntz, plan[i].d_u0_inv, plan[i].d_u1_inv, plan[i].d_u2_inv, plan[i].d_uk0_inv, plan[i].d_uk_inv,
			plan[i].d_Lap, plan[i].d_amp_Lap, plan[i].d_pha_Lap);

		// initializating last two wavefields for backward propagation
		cuda_kernel_initialization_Finals<<<dimGrid,dimBlock,0,plans[i].stream>>>
			(ntx, ntz, plan[i].d_u0_inv, plan[i].d_u1_inv, plan[i].d_u2_final0, plan[i].d_u2_final1);

		// initialization wavefield for backward propagation
		cuda_kernel_initialization<<<dimGrid,dimBlock,0,plans[i].stream>>>
			(ntx, ntz, plan[i].d_u0, plan[i].d_u1, plan[i].d_u2, plan[i].d_uk0, plan[i].d_uk,
			plan[i].d_Lap, plan[i].d_amp_Lap, plan[i].d_pha_Lap);

		// initialization image variables for imaging
		cuda_kernel_initialization_images<<<dimGrid,dimBlock,0,plans[i].stream>>>
			(ntx, ntz, plan[i].d_image_cor, plan[i].d_image_nor, plan[i].d_image_sources, plan[i].d_image_receivers);			
	}

	// copy the vectors from the host to the device
	for(i=0;i<GPU_N;i++)
	{
		cudaSetDevice(i);
		cudaMemcpyAsync(plan[i].d_seismogram_rms,plan[i].seismogram_rms,
			sizeof(float)*rnmax*nt,cudaMemcpyHostToDevice,plans[i].stream);
	}

	int beta1_inv=beta1;
	int beta2_inv=-1*beta2;		// source wavefield is compensated, so reconstruction process is attenuated
	float kx_cut_inv=kx_cut;
	float kz_cut_inv=kz_cut;


	// backward time iteration
	for(it=nt-3;it>=0;it--)  
	{
		for(i=0;i<GPU_N;i++)
		{
			cudaSetDevice(i);

			// reconstructing source wavefield using checkpointing scheme
			if(Sto_Rec==0)
			{
				cufftExecC2C(plan[i].PLAN_FORWARD,plan[i].d_u1_inv,plan[i].d_uk_inv,CUFFT_FORWARD); //CUFFT_FORWARD
				cuda_kernel_visco_PSM_2d_reconstruction_k_space<<<dimGrid,dimBlock,0,plans[i].stream>>>
					(
						beta1_inv, beta2_inv,
						it, nt, ntx, ntz, nx, nz, L, dx, dz, dt, 
						plan[i].d_vp, plan[i].d_Gamma, averGamma, f0, Omega0,
						plan[i].d_kx, plan[i].d_kz,  
						plan[i].d_uk_inv, plan[i].d_uk0_inv, 
						plan[i].d_Lap_uk, plan[i].d_amp_uk, plan[i].d_pha_uk
					);	

				cufftExecC2C(plan[i].PLAN_BACKWARD,plan[i].d_Lap_uk,plan[i].d_Lap,CUFFT_INVERSE); //CUFFT_INVERSE

				if(beta1_inv!=0)
				{	
					cufftExecC2C(plan[i].PLAN_BACKWARD,plan[i].d_pha_uk,plan[i].d_pha_Lap,CUFFT_INVERSE); //CUFFT_INVERSE					
				}
		
				if(beta2_inv!=0)
				{
					if(beta2_inv<0&&Filtertype==0)				// low-pass filtering
					{
						cuda_kernel_filter2d<<<dimGrid,dimBlock,0,plans[i].stream>>>
							(ntx, ntz, dx, dz, plan[i].d_kfilter, plan[i].d_kx, plan[i].d_kz, plan[i].d_amp_uk, kx_cut_inv, kz_cut_inv, taper_ratio);
					}

					cufftExecC2C(plan[i].PLAN_BACKWARD,plan[i].d_amp_uk,plan[i].d_amp_Lap,CUFFT_INVERSE); //CUFFT_INVERSE
				}
				
				cuda_kernel_visco_PSM_2d_reconstruction_x_space<<<dimGrid,dimBlock,0,plans[i].stream>>>
					(
						beta1_inv, beta2_inv,
						it, nt, ntx, ntz, nx, nz, L, dx, dz, dt, 
						plan[i].d_vp, plan[i].d_Gamma, averGamma, f0, Omega0,
						plan[i].d_ricker, ss[is+i].s_ix, ss[is+i].s_iz,
						plan[i].d_u0_inv, plan[i].d_u1_inv, plan[i].d_u2_inv,
						plan[i].d_Lap, plan[i].d_amp_Lap, plan[i].d_pha_Lap,			
						plan[i].d_borders_up, plan[i].d_borders_bottom, plan[i].d_borders_left, plan[i].d_borders_right
					);	

				
				if(beta2_inv<0&&Filtertype==1)					// adaptive stabilization scheme
				{
					cufftExecC2C(plan[i].PLAN_FORWARD,plan[i].d_u2_inv,plan[i].d_uk_inv,CUFFT_FORWARD); //CUFFT_FORWARD
					cuda_kernel_AdaSta<<<dimGrid,dimBlock,0,plans[i].stream>>>
					(
						nt-it, ntx, ntz, dx, dz, dt, plan[i].d_kstabilization,
						plan[i].d_vp, plan[i].d_Gamma, avervp, averGamma, Omega0,
						plan[i].d_kx, plan[i].d_kz, plan[i].d_uk_inv, sigma, Order
					);	
					cufftExecC2C(plan[i].PLAN_BACKWARD,plan[i].d_uk_inv,plan[i].d_u2_inv,CUFFT_INVERSE); //CUFFT_INVERSE	
				}


				// MTF absorbing boundary condition
				cuda_kernel_MTF_2nd<<<dimGrid,dimBlock,0,plans[i].stream>>>
					(L, ntx, ntz, dx, dz, dt, plan[i].d_vp, plan[i].d_u0_inv, plan[i].d_u1_inv, plan[i].d_u2_inv);

				// read wavefields at checkpoints
				cuda_kernel_checkpoints_In<<<dimGrid,dimBlock,0,plans[i].stream>>>
					(
						it, nt, ntx, ntz, nx, nz, L, dx, dz, dt, 
						plan[i].d_u1_inv, plan[i].d_u2_inv,
						plan[i].d_u_cp, N_cp, plan[i].d_t_cp			
					);							

				// updating wavefields
				cuda_kernel_update<<<dimGrid,dimBlock,0,plans[i].stream>>>
					(ntx, ntz, plan[i].d_u0_inv, plan[i].d_u1_inv, plan[i].d_u2_inv);
			}

			// read source wavefields from disk 
			if(Sto_Rec==1)
			{
				sprintf(filename,"./output/GPU_%d_u2_%d.dat",i,it); 
				fp=fopen(filename,"rb");
				for(ix=0;ix<ntx-0;ix++)
				{
					for(iz=0;iz<ntz-0;iz++)
					{								
						fread(&u2_real[iz*ntx+ix],sizeof(float),1,fp);
						plan[i].u2[iz*ntx+ix].x=u2_real[iz*ntx+ix];	
						plan[i].u2[iz*ntx+ix].y=0.0;
					}
				}
				fclose(fp);
				cudaMemcpyAsync(plan[i].d_u2_inv,plan[i].u2,sizeof(cufftComplex)*ntp,cudaMemcpyHostToDevice,plans[i].stream);
			}




			// backward propagation for imaging
			cufftExecC2C(plan[i].PLAN_FORWARD,plan[i].d_u1,plan[i].d_uk,CUFFT_FORWARD); //CUFFT_FORWARD

			cuda_kernel_visco_PSM_2d_backward_k_space<<<dimGrid,dimBlock,0,plans[i].stream>>>
				(
					beta1, beta2,
					it, nt, ntx, ntz, dx, dz, dt, 
					plan[i].d_vp, plan[i].d_Gamma, averGamma, f0, Omega0,
					plan[i].d_kx, plan[i].d_kz,  
					plan[i].d_uk, plan[i].d_uk0, 
					plan[i].d_Lap_uk, plan[i].d_amp_uk, plan[i].d_pha_uk
				);

			cufftExecC2C(plan[i].PLAN_BACKWARD,plan[i].d_Lap_uk,plan[i].d_Lap,CUFFT_INVERSE); //CUFFT_INVERSE

			if(beta1!=0)
			{	
				cufftExecC2C(plan[i].PLAN_BACKWARD,plan[i].d_pha_uk,plan[i].d_pha_Lap,CUFFT_INVERSE); //CUFFT_INVERSE					
			}

			if(beta2!=0)
			{
				if (beta2<0&&Filtertype==0)					// low-pass filtering
				{
					cuda_kernel_filter2d<<<dimGrid,dimBlock,0,plans[i].stream>>>
						(ntx, ntz, dx, dz, plan[i].d_kfilter, plan[i].d_kx, plan[i].d_kz, plan[i].d_amp_uk, kx_cut, kz_cut, taper_ratio);
				}

				cufftExecC2C(plan[i].PLAN_BACKWARD,plan[i].d_amp_uk,plan[i].d_amp_Lap,CUFFT_INVERSE); //CUFFT_INVERSE
			}

			cuda_kernel_visco_PSM_2d_backward_x_space<<<dimGrid,dimBlock,0,plans[i].stream>>>
				(
					beta1, beta2,
					it, nt, ntx, ntz, dx, dz, dt, 
					plan[i].d_vp, plan[i].d_Gamma, averGamma, f0, Omega0,
					plan[i].d_seismogram_rms, plan[i].d_r_ix, ss[is+i].r_iz, ss[is+i].s_ix, rnmax, nrx_obs,
					plan[i].d_u0, plan[i].d_u1, plan[i].d_u2,
					plan[i].d_Lap, plan[i].d_amp_Lap, plan[i].d_pha_Lap
				);
				
			if(beta2<0&&Filtertype==1)						// adaptive stabilization scheme
			{
				cufftExecC2C(plan[i].PLAN_FORWARD,plan[i].d_u2,plan[i].d_uk,CUFFT_FORWARD); //CUFFT_FORWARD
				cuda_kernel_AdaSta<<<dimGrid,dimBlock,0,plans[i].stream>>>
				(
					nt-it, ntx, ntz, dx, dz, dt, plan[i].d_kstabilization,
					plan[i].d_vp, plan[i].d_Gamma, avervp, averGamma, Omega0,
					plan[i].d_kx, plan[i].d_kz, plan[i].d_uk, sigma, Order
				);	
				cufftExecC2C(plan[i].PLAN_BACKWARD,plan[i].d_uk,plan[i].d_u2,CUFFT_INVERSE); //CUFFT_INVERSE	
			}
			
			// MTF absorbing boundary condition
			cuda_kernel_MTF_2nd<<<dimGrid,dimBlock,0,plans[i].stream>>>
				(L, ntx, ntz, dx, dz, dt, plan[i].d_vp, plan[i].d_u0, plan[i].d_u1, plan[i].d_u2);
		
			// imaging (exclude duration of explosion)
			int it0 = int(2/(f0*dt));
			if (it>it0)
			{
				cuda_kernel_image<<<dimGrid,dimBlock,0,plans[i].stream>>>
					(
						ntx, ntz, L,
						plan[i].d_u2_inv, plan[i].d_u2,
						plan[i].d_image_cor, plan[i].d_image_sources, plan[i].d_image_receivers
					);			
			}


			// write backward wavefields and reconstructed wavefields to disk
			if(Save_Not==1)
			{
				cudaMemcpyAsync(plan[i].u1,plan[i].d_u2_inv,sizeof(cufftComplex)*ntp,cudaMemcpyDeviceToHost,plans[i].stream);
				cudaMemcpyAsync(plan[i].u2,plan[i].d_u2,sizeof(cufftComplex)*ntp,cudaMemcpyDeviceToHost,plans[i].stream);
				cudaStreamSynchronize(plans[i].stream);
				sprintf(filename,"./output/GPU_%d_u2_inv_%d.dat",i,it); 
				fp=fopen(filename,"wb");
				for(ix=0;ix<ntx-0;ix++)
				{
					for(iz=0;iz<ntz-0;iz++)
					{
						u2_real[iz*ntx+ix]=plan[i].u1[iz*ntx+ix].x;
						fwrite(&u2_real[iz*ntx+ix],sizeof(float),1,fp);
					}
				}
				fclose(fp);
				sprintf(filename,"./output/GPU_%d_u2_bak_%d.dat",i,it); 
				fp=fopen(filename,"wb");
				for(ix=0;ix<ntx-0;ix++)
				{
					for(iz=0;iz<ntz-0;iz++)
					{	
						u2_real[iz*ntx+ix]=plan[i].u2[iz*ntx+ix].x;				
						fwrite(&u2_real[iz*ntx+ix],sizeof(float),1,fp);
					}
				}
				fclose(fp);
			}

			// updating wavefields 
			cuda_kernel_update<<<dimGrid,dimBlock,0,plans[i].stream>>>
				(ntx, ntz, plan[i].d_u0, plan[i].d_u1, plan[i].d_u2);

			/*if(myid==0&&it%100==0)
			{
				sf_warning("shot %d reconstruction and backward %d has finished!", is+i+1, it);
			}*/
		}// GPU_N end	
	}// nt end

	
	sf_warning("shot %d reconstruction and backward has finished!", is+i);

	for(i=0;i<GPU_N;i++)
	{
		cudaSetDevice(i);
		// output images 
		cudaMemcpyAsync(plan[i].image_cor,plan[i].d_image_cor,sizeof(float)*ntp,cudaMemcpyDeviceToHost,plans[i].stream);
		cudaMemcpyAsync(plan[i].image_sources,plan[i].d_image_sources,sizeof(float)*ntp,cudaMemcpyDeviceToHost,plans[i].stream);
		cudaMemcpyAsync(plan[i].image_receivers,plan[i].d_image_receivers,sizeof(float)*ntp,cudaMemcpyDeviceToHost,plans[i].stream);
	}

	for(i=0;i<GPU_N;i++)
	{
		cudaSetDevice(i);
		cudaDeviceSynchronize();
	}

	//free the memory of DEVICE
	for(i=0;i<GPU_N;i++)
	{
		cudaSetDevice(i);
		cudaStreamDestroy(plans[i].stream);
	}
	free(u2_real);	
}


//==========================================================
//  This two subroutines are used for Laplace filteing
//  ========================================================
extern "C"
void Laplace_filtering
(
	float *image, int ntx, int ntz, float dx, float dz
)
{ 
	int ix,iz,ip,K,NX,NZ;

	K=(int)ceil(log(1.0*ntx)/log(2.0));
	NX=(int)pow(2.0,K);

	K=(int)ceil(log(1.0*ntz)/log(2.0));
	NZ=(int)pow(2.0,K);

	float dkx,dkz;
	float kx,kz;

	dkx=(float)1.0/((NX)*dx);
	dkz=(float)1.0/((NZ)*dz);

	int NTP=NX*NZ;

	cufftComplex *pp,*temp,*tempout;		

	cudaMallocHost((void **)&pp, sizeof(cufftComplex)*NX*NZ);
	cudaMalloc((void **)&temp,sizeof(cufftComplex)*NX*NZ);
	cudaMalloc((void **)&tempout,sizeof(cufftComplex)*NX*NZ);

	cufftHandle plan;
	cufftPlan2d(&plan,NX,NZ,CUFFT_C2C);

	for(ip=0;ip<NTP;ip++)
	{ 
		pp[ip].x=0.0;
		pp[ip].y=0.0; 
	} 

	for(ix=0;ix<ntx;ix++)
	{            
		for(iz=0;iz<ntz;iz++)
		{
			pp[ix*NZ+iz].x=image[iz*ntx+ix];
		}
	} 

	cudaMemcpy(temp,pp,sizeof(cufftComplex)*NX*NZ,cudaMemcpyHostToDevice);
	cufftExecC2C(plan,temp,tempout,CUFFT_FORWARD);
	cudaMemcpy(pp,tempout,sizeof(cufftComplex)*NX*NZ,cudaMemcpyDeviceToHost);

	for(ix=0;ix<NX;ix++)
	{            
		for(iz=0;iz<NZ;iz++)
		{
			if(ix<NX/2)
			{
				kx=2*PI*ix*dkx;
			}
			if(ix>NX/2)	
			{
				kx=2*PI*(NX-1-ix)*dkx;
			}

			if(iz<NZ/2)
			{
				kz=2*PI*iz*dkz;//2*PI*(NZ/2-1-iz)*dkz;//0.0;//
			}
			if(iz>NZ/2)
			{
				kz=2*PI*(NZ-1-iz)*dkz;//2*PI*(iz-NZ/2)*dkz;//0.0;//
			}

			ip=ix*NZ+iz;

			pp[ip].x=pp[ip].x*(kx*kx+kz*kz);
			pp[ip].y=pp[ip].y*(kx*kx+kz*kz);

		}
	} 

	cudaMemcpy(temp,pp,sizeof(cufftComplex)*NX*NZ,cudaMemcpyHostToDevice);
	cufftExecC2C(plan,temp,tempout,CUFFT_INVERSE);
	cudaMemcpy(pp,tempout,sizeof(cufftComplex)*NX*NZ,cudaMemcpyDeviceToHost);

	for(ix=0;ix<ntx;ix++)
	{            
		for(iz=0;iz<ntz;iz++)
		{
			image[iz*ntx+ix]=pp[ix*NZ+iz].x/(NX*NZ);
		}
	} 
	cudaFreeHost(pp);
	cudaFree(temp);
	cudaFree(tempout);
	cufftDestroy(plan);

	return;
}


//=========================================================
//  Allocate the memory for variables in device
//  =======================================================
extern "C"
void cuda_Device_malloc
(
	int ntx, int ntz, int ntp, int nx, int nz, int nt, 
	float dx, float dz, int L, int rnmax, int N_cp,
	struct MultiGPU plan[], int GPU_N
)
{
	int i;
	size_t size_model=sizeof(float)*ntp;

	for(i=0;i<GPU_N;i++)
	{
		cudaSetDevice(i);

		cufftPlan2d(&plan[i].PLAN_FORWARD,ntz,ntx,CUFFT_C2C);
		cufftPlan2d(&plan[i].PLAN_BACKWARD,ntz,ntx,CUFFT_C2C);

		cudaMallocHost((void **)&plan[i].u0, sizeof(cufftComplex)*ntp);		
		cudaMallocHost((void **)&plan[i].u1, sizeof(cufftComplex)*ntp);
		cudaMallocHost((void **)&plan[i].u2, sizeof(cufftComplex)*ntp);		

		cudaMallocHost((void **)&plan[i].seismogram_obs, sizeof(float)*nt*rnmax);
		cudaMallocHost((void **)&plan[i].seismogram_dir, sizeof(float)*nt*rnmax);
		cudaMallocHost((void **)&plan[i].seismogram_syn, sizeof(float)*nt*rnmax);
		cudaMallocHost((void **)&plan[i].seismogram_rms, sizeof(float)*nt*rnmax);
		
		cudaMallocHost((void **)&plan[i].image_sources, sizeof(float)*ntp);
		cudaMallocHost((void **)&plan[i].image_receivers, sizeof(float)*ntp);
		cudaMallocHost((void **)&plan[i].image_cor, sizeof(float)*ntp);
		cudaMallocHost((void **)&plan[i].image_nor, sizeof(float)*ntp);

		cudaMalloc((void**)&plan[i].d_r_ix,sizeof(int)*rnmax);
		cudaMalloc((void**)&plan[i].d_ricker,sizeof(float)*nt);        //ricker

		cudaMalloc((void**)&plan[i].d_vp,size_model);
		cudaMalloc((void**)&plan[i].d_Gamma,size_model);

		cudaMalloc((void**)&plan[i].d_u0,sizeof(cufftComplex)*ntp);
		cudaMalloc((void**)&plan[i].d_u1,sizeof(cufftComplex)*ntp);
		cudaMalloc((void**)&plan[i].d_u2,sizeof(cufftComplex)*ntp);

		cudaMalloc((void**)&plan[i].d_u0_inv,sizeof(cufftComplex)*ntp);
		cudaMalloc((void**)&plan[i].d_u1_inv,sizeof(cufftComplex)*ntp);
		cudaMalloc((void**)&plan[i].d_u2_inv,sizeof(cufftComplex)*ntp);

		cudaMalloc((void**)&plan[i].d_t_cp,sizeof(int)*N_cp);
		cudaMalloc((void**)&plan[i].d_u_cp,size_model*N_cp);		//checkpoints

		cudaMalloc((void**)&plan[i].d_kx,sizeof(float)*ntx);
		cudaMalloc((void**)&plan[i].d_kz,sizeof(float)*ntz);
		cudaMalloc((void**)&plan[i].d_kfilter,sizeof(float)*ntp);
		cudaMalloc((void**)&plan[i].d_kstabilization,sizeof(float)*ntp);

		cudaMalloc((void **)&plan[i].d_uk,sizeof(cufftComplex)*ntp);
		cudaMalloc((void **)&plan[i].d_uk0,sizeof(cufftComplex)*ntp);

		cudaMalloc((void **)&plan[i].d_uk_inv,sizeof(cufftComplex)*ntp);
		cudaMalloc((void **)&plan[i].d_uk0_inv,sizeof(cufftComplex)*ntp);

		cudaMalloc((void **)&plan[i].d_Lap_uk,sizeof(cufftComplex)*ntp);
		cudaMalloc((void **)&plan[i].d_amp_uk,sizeof(cufftComplex)*ntp);
		cudaMalloc((void **)&plan[i].d_pha_uk,sizeof(cufftComplex)*ntp);

		cudaMalloc((void **)&plan[i].d_Lap,sizeof(cufftComplex)*ntp);
		cudaMalloc((void **)&plan[i].d_amp_Lap,sizeof(cufftComplex)*ntp);
		cudaMalloc((void **)&plan[i].d_pha_Lap,sizeof(cufftComplex)*ntp);
							
		cudaMalloc((void**)&plan[i].d_seismogram,sizeof(float)*nt*rnmax);
		cudaMalloc((void**)&plan[i].d_seismogram_rms,sizeof(float)*nt*rnmax);

		cudaMalloc((void**)&plan[i].d_borders_up,sizeof(float)*nt*nx);
		cudaMalloc((void**)&plan[i].d_borders_bottom,sizeof(float)*nt*nx);
		cudaMalloc((void**)&plan[i].d_borders_left,sizeof(float)*nt*nz);
		cudaMalloc((void**)&plan[i].d_borders_right,sizeof(float)*nt*nz);

		cudaMalloc((void**)&plan[i].d_u2_final0,size_model);
		cudaMalloc((void**)&plan[i].d_u2_final1,size_model);

		cudaMalloc((void**)&plan[i].d_image_sources,size_model);
		cudaMalloc((void**)&plan[i].d_image_receivers,size_model);
		cudaMalloc((void**)&plan[i].d_image_cor,size_model);
		cudaMalloc((void**)&plan[i].d_image_nor,size_model);
	}
}


//=========================================================
//  Free the memory for variables in device
//  =======================================================
extern "C"
void cuda_Device_free
(
	int ntx, int ntz, int ntp, int nx, int nz, int nt, 
	float dx, float dz, int L, int rnmax, int N_cp,
	struct MultiGPU plan[], int GPU_N
)
{
	int i;
	 

	for(i=0;i<GPU_N;i++)
	{
		cudaSetDevice(i);

		cufftDestroy(plan[i].PLAN_FORWARD);
		cufftDestroy(plan[i].PLAN_BACKWARD);

		cudaFreeHost(plan[i].u0);
		cudaFreeHost(plan[i].u1);
		cudaFreeHost(plan[i].u2); 

		cudaFreeHost(plan[i].seismogram_obs);
		cudaFreeHost(plan[i].seismogram_dir);
		cudaFreeHost(plan[i].seismogram_syn); 
		cudaFreeHost(plan[i].seismogram_rms);
 
		cudaFreeHost(plan[i].image_cor);
		cudaFreeHost(plan[i].image_nor);
		cudaFreeHost(plan[i].image_sources);
		cudaFreeHost(plan[i].image_receivers);

		cudaFree(plan[i].d_r_ix);
		cudaFree(plan[i].d_ricker);

		cudaFree(plan[i].d_vp);
		cudaFree(plan[i].d_Gamma);

		cudaFree(plan[i].d_u0);
		cudaFree(plan[i].d_u1);
		cudaFree(plan[i].d_u2);

		cudaFree(plan[i].d_u0_inv);
		cudaFree(plan[i].d_u1_inv);
		cudaFree(plan[i].d_u2_inv);

		cudaFree(plan[i].d_t_cp);
		cudaFree(plan[i].d_u_cp);

		cudaFree(plan[i].d_kx);
		cudaFree(plan[i].d_kz);
		cudaFree(plan[i].d_kfilter);
		cudaFree(plan[i].d_kstabilization);

		cudaFree(plan[i].d_uk);
		cudaFree(plan[i].d_uk0);

		cudaFree(plan[i].d_uk_inv);
		cudaFree(plan[i].d_uk0_inv);

		cudaFree(plan[i].d_Lap_uk);
		cudaFree(plan[i].d_amp_uk);
		cudaFree(plan[i].d_pha_uk);

		cudaFree(plan[i].d_Lap);
		cudaFree(plan[i].d_amp_Lap);
		cudaFree(plan[i].d_pha_Lap);

		cudaFree(plan[i].d_seismogram);
		cudaFree(plan[i].d_seismogram_rms);

		cudaFree(plan[i].d_borders_up);
		cudaFree(plan[i].d_borders_bottom);
		cudaFree(plan[i].d_borders_left);
		cudaFree(plan[i].d_borders_right);

		cudaFree(plan[i].d_u2_final0);
		cudaFree(plan[i].d_u2_final1);

		cudaFree(plan[i].d_image_sources);
		cudaFree(plan[i].d_image_receivers);
		cudaFree(plan[i].d_image_cor);
		cudaFree(plan[i].d_image_nor);
	}
}


//=========================================================
//  Initializating the memory for variables in device
//  =======================================================
extern "C"
void cuda_Host_initialization
(
	int ntx, int ntz, int ntp, int nx, int nz, int nt, 
	float dx, float dz, int L, int rnmax, int N_cp,
	struct MultiGPU plan[], int GPU_N
)
{
	int i;
	for(i=0;i<GPU_N;i++)
	{
		cudaSetDevice(i);
		memset(plan[i].u0, 0, ntx*ntz*sizeof(float));
		memset(plan[i].u1, 0, ntx*ntz*sizeof(float));
		memset(plan[i].u2, 0, ntx*ntz*sizeof(float));
		memset(plan[i].seismogram_obs, 0, nt*rnmax*sizeof(float));
		memset(plan[i].seismogram_dir, 0, nt*rnmax*sizeof(float));
		memset(plan[i].seismogram_syn, 0, nt*rnmax*sizeof(float));
		memset(plan[i].seismogram_rms, 0, nt*rnmax*sizeof(float));
 		memset(plan[i].image_cor, 0, ntx*ntz*sizeof(float));
 		memset(plan[i].image_nor, 0, ntx*ntz*sizeof(float));
 		memset(plan[i].image_sources, 0, ntx*ntz*sizeof(float));
 		memset(plan[i].image_receivers, 0, ntx*ntz*sizeof(float));
	}
}

extern "C"
void getdevice(int *GPU_N)
{	
	cudaGetDeviceCount(GPU_N);	
}


//==========================================================
//  This subroutine is used for calculating the ricker wave
//  ========================================================
void ricker_wave
(
	float *ricker, int nt, float f0, float t0, float dt, int flag
)
{
	float pi=3.1415927;
	int   it;
	float temp,max=0.0;
	FILE *fp;
	if(flag==1)
	{
		for(it=0;it<nt;it++)
		{
			temp=pi*f0*(it*dt-t0);
			temp=temp*temp;
			ricker[it]=(1.0-2.0*temp)*exp(-temp);
		}
		fp=fopen("./output/ricker.dat","wb");    
		for(it=0;it<nt;it++)
		{
			fwrite(&ricker[it],sizeof(float),1,fp);
		}    
		fclose(fp);
	}
	if(flag==2)
	{
		for(it=0;it<nt;it++)
		{
			temp=pi*f0*(it*dt-t0);
			temp=temp*temp;         
			ricker[it]=(it*dt-t0)*exp(-temp);

			if(max<fabs(ricker[it]))
			{
				max=fabs(ricker[it]);
			}
		}
		for(it=0;it<nt;it++)
		{
			ricker[it]=ricker[it]/max;
		}
		fp=fopen("./output/ricker_integration.dat","wb");    
		for(it=0;it<nt;it++)
		{
			fwrite(&ricker[it],sizeof(float),1,fp);
		}    
		fclose(fp);
	}
	if(flag==3)
	{	
		for(it=0;it<nt;it++)
		{
			temp=pi*f0*(it*dt-t0);
			ricker[it]=(4*powf(pi*f0,4)*powf((it*dt-t0),3)-6*powf(pi*f0,2)*(it*dt-t0))*exp(-powf(temp,2));  

			if(max<fabs(ricker[it]))
			{
				max=fabs(ricker[it]);
			}
		}
		for(it=0;it<nt;it++)
		{
			ricker[it]=ricker[it]/max;
		}
		fp=fopen("./output/ricker_derivative.dat","wb");    
		for(it=0;it<nt;it++)
		{
			fwrite(&ricker[it],sizeof(float),1,fp);
		}    
		fclose(fp);
	}
	return;
}


/*==========================================================
  This subroutine is used for initializing the true model...
  ===========================================================*/
void get_acc_model
(
	float *vp0, float *Qp0, float *vp, float *Qp, int ntp, int ntx, int ntz, int L
)
{
	int ip,ipp,iz,ix,ip0; 
	FILE *fp;

	for(ix=L;ix<ntx-L;ix++)
	{
		for(iz=L;iz<ntz-L;iz++)
		{
			ip=iz*ntx+ix;
			ip0=(iz-L)*(ntx-2*L)+(ix-L);
			vp[ip]=vp0[ip0];
		}
	}

	for(iz=0;iz<=L-1;iz++)
	{
		for(ix=0;ix<=L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=L*ntx+L;

			vp[ip]=vp[ipp];
		}
		for(ix=L;ix<=ntx-L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=L*ntx+ix;

			vp[ip]=vp[ipp];
		}
		for(ix=ntx-L;ix<ntx;ix++)
		{
			ip=iz*ntx+ix;
			ipp=L*ntx+ntx-L-1;

			vp[ip]=vp[ipp];
		}
	}
	for(iz=L;iz<=ntz-L-1;iz++)
	{
		for(ix=0;ix<=L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=iz*ntx+L;

			vp[ip]=vp[ipp];
		}
		for(ix=ntx-L;ix<ntx;ix++)
		{
			ip=iz*ntx+ix;
			ipp=iz*ntx+ntx-L-1;

			vp[ip]=vp[ipp];
		}
	}

	for(iz=ntz-L;iz<ntz;iz++)
	{
		for(ix=0;ix<=L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=(ntz-L-1)*ntx+L;

			vp[ip]=vp[ipp];
		}
		for(ix=L;ix<=ntx-L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=(ntz-L-1)*ntx+ix;

			vp[ip]=vp[ipp];
		}
		for(ix=ntx-L;ix<ntx;ix++)
		{
			ip=iz*ntx+ix;
			ipp=(ntz-L-1)*ntx+ntx-L-1;
			vp[ip]=vp[ipp];
		}
	}

	for(ix=L;ix<ntx-L;ix++)
	{
		for(iz=L;iz<ntz-L;iz++)
		{
			ip=iz*ntx+ix;
			ip0=(iz-L)*(ntx-2*L)+(ix-L);
			Qp[ip]=Qp0[ip0];
		}
	}


	

	for(iz=0;iz<=L-1;iz++)
	{
		for(ix=0;ix<=L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=L*ntx+L;
			Qp[ip]=Qp[ipp];
		}
		for(ix=L;ix<=ntx-L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=L*ntx+ix;
			Qp[ip]=Qp[ipp];
		}
		for(ix=ntx-L;ix<ntx;ix++)
		{
			ip=iz*ntx+ix;
			ipp=L*ntx+ntx-L-1;
			Qp[ip]=Qp[ipp];
		}
	}
	for(iz=L;iz<=ntz-L-1;iz++)
	{
		for(ix=0;ix<=L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=iz*ntx+L;
			Qp[ip]=Qp[ipp];
		}
		for(ix=ntx-L;ix<ntx;ix++)
		{
			ip=iz*ntx+ix;
			ipp=iz*ntx+ntx-L-1;
			Qp[ip]=Qp[ipp];
		}
	}

	for(iz=ntz-L;iz<ntz;iz++)
	{
		for(ix=0;ix<=L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=(ntz-L-1)*ntx+L;
			Qp[ip]=Qp[ipp];
		}
		for(ix=L;ix<=ntx-L-1;ix++)
		{
			ip=iz*ntx+ix;
			ipp=(ntz-L-1)*ntx+ix;
			Qp[ip]=Qp[ipp];
		}
		for(ix=ntx-L;ix<ntx;ix++)
		{
			ip=iz*ntx+ix;
			ipp=(ntz-L-1)*ntx+ntx-L-1;
			Qp[ip]=Qp[ipp];
		}
	}
	return;
}


// ==========================================================
//  This subroutine is used for initializing the homogeneous model...
//  =========================================================
void get_homo_model
(
	float *vp, int ntp, int ntx, int ntz, int L
)
{
	int ip,ipp,iz,ix;  
	FILE *fp;
	for(ix=0;ix<ntx;ix++)
	{
		for(iz=0;iz<ntz;iz++)
		{
			ip=iz*ntx+ix;
			if(iz>L+1)
			{
				ipp=(L+1)*ntx+ix;
				vp[ip]=vp[ipp];
			}
		}
	}
	return;
}


//==========================================================
//  This subroutine is used for initializing the initial model...
//  ========================================================
void get_ini_model
(
	float *vp, int ntp, int ntx, int ntz, int span
)
{
	int ix, ixw, ixx;
	int iz, izw, izz;

	float *s_a;
	s_a=(float*)malloc(sizeof(float)*ntp);

	float *sx_a;
	sx_a=(float*)malloc(sizeof(float)*ntp);

	for(iz=0;iz<ntz;iz++)
	{
		for(ix=0; ix<ntx; ix++)
		{
			s_a[iz*ntx+ix]=0.0;
			for(ixw=ix-span;ixw<=ix+span;ixw++)
			{
				if(ixw<0)
					ixx=0;
				else if(ixw>ntx-1)
					ixx=ntx-1;
				else
					ixx=ixw;
				s_a[iz*ntx+ix]+=vp[iz*ntx+ixx]/(2*span+1);
			}		
		}	
	}

	for(iz=0;iz<ntz;iz++)
	{
		for(ix=0; ix<ntx; ix++)
		{		
			sx_a[iz*ntx+ix]=s_a[iz*ntx+ix];		
		}	
	}

	for(iz=0;iz<ntz;iz++)
	{
		for(ix=0; ix<ntx; ix++)
		{
			s_a[iz*ntx+ix]=0.0;
			for(izw=iz-span;izw<=iz+span;izw++)
			{
				if(izw<0)
					izz=0;
				else if(izw>ntz-1)
					izz=ntz-1;
				else
					izz=izw;
				s_a[iz*ntx+ix]+=sx_a[izz*ntx+ix]/(2*span+1);
			}		
		}	
	}

	for(iz=0;iz<ntz;iz++)
	{
		for(ix=0; ix<ntx; ix++)
		{
			vp[iz*ntx+ix]=s_a[iz*ntx+ix];		
		}	
	}

	free(s_a);	
	free(sx_a);		
}


//==========================================================
//  This subroutine is used for cutting the direct wave
//  =======================================================

void cut_dir
(
	float *seismogram_obs, float *seismogram_rms, 
	int rnmax, int nt, int is, float dx, float dz, float dt, 
	int r_iz, int s_ix, int s_iz, float t0, float *vp
)
{
	int it, ix;
	float removeline[rnmax];

	for(it=0;it<nt; it++)
	{
		for(ix=0; ix<rnmax; ix++)
		{
			removeline[ix]=(sqrt(powf((ix-s_ix)*dx,2)+powf((r_iz-s_iz)*dz,2))/vp[1*rnmax+ix]+4.0*t0)/dt;

			if(it<removeline[ix])
				seismogram_rms[it*rnmax+ix]=0.0;
			else
				seismogram_rms[it*rnmax+ix]=seismogram_obs[it*rnmax+ix];
		}	
	}
}

//==========================================================
//  This subroutine is used for Laplace filtering
//  ========================================================
void Laplace_FD_filtering
(
	float *image, int ntx, int ntz, float dx, float dz
)
{ 
	int ix,iz,ip;
	float diff1, diff2;
	float *tmp;
	tmp = (float*)malloc(sizeof(float)*ntx*ntz);
	memset(tmp, 0, ntx*ntz*sizeof(float));
	for(iz=1;iz<ntz-1;iz++)
	{
		for(ix=1;ix<ntx-1;ix++)
		{
			ip=iz*ntx+ix;
			diff1=(image[ip+ntx]-2.0*image[ip]+image[ip-ntx])/(dz*dz);
			diff2=(image[ip+1]-2.0*image[ip]+image[ip-1])/(dx*dx);	
			tmp[ip]=diff1+diff2;          
		}
	}

	for(iz=0;iz<=ntz-1;iz++)
	{
		for(ix=0;ix<=ntx-1;ix++)
		{
			ip=iz*ntx+ix;
			image[ip]=tmp[ip];          
		}
	}	
	free(tmp);
	return;
}


int main(int argc,char *argv[])
{
	sf_init(argc,argv);

	/*NEW PROGRAM ENDS HERE*/
	sf_file in, inq, out, out2;
        /* set I/O files */
        in=sf_input("in");        /*  velocity  */
	inq=sf_input("q");	  /* Q */

        out=sf_output("out");      /*  Image with correlation IC  */
        out2=sf_output("out2");      /* Image with normalization IC  */

	/*NEW PROGRAM ENDS HERE*/

	int myid,numprocs,namelen;	
	MPI_Comm comm=MPI_COMM_WORLD;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(comm,&myid);
	MPI_Comm_size(comm,&numprocs);
	MPI_Get_processor_name(processor_name,&namelen);
	if(myid==0)
		sf_warning("Number of MPI thread is %d",numprocs);


	//=========================================================
	//  Parameters of the time of the system...
	//  =======================================================

	clock_t start, end;
	float runtime=0.0;
	int RTMtype;

	//=========================================================
	//  Parameters of Cartesian coordinate...
	//  =======================================================

	int it,ix,iz,ip,ipp;
	int nz;				// vertical samples
	int nx;				// horizontal samples
	int L;				// Layers of MTF absorbing boundary condition
	int ntz;			// total vertical samples 
	int ntx;			// total horizontal samples 
	int nxh;			// half of vertical samples
	int nzh;			// half of horizontal samples
	int ntp;		// total samples of the model
	int np;			// total samples of the simulating domain
	float dz;			// vertical interval
	float dx;			// horizontal interval
	FILE *fp;				// file pointer
	char filename[40];		// filename char



	/*=========================================================
	  Parameters of ricker wave...
	  ========================================================*/

	float f0;			// dominate frequency of the ricker
	int nt;			// time samples
	float dt;			// time interval
	float *ricker;			// ricker wavelet
	ricker=(float *) malloc(nt*sizeof(float));



        if (!sf_getint("type",&RTMtype)) RTMtype=3;	
	/*RTM type parameter, 0 for acoustic RTM, 1 for viscoacoustic RTM, 2 for QRTM using low-pass filtering, 3 for QRTM using adaptive stabilization scheme*/

        if (!sf_getint("nz",&nz)) nz=234;	
	/*NZ*/

        if (!sf_getint("nx",&nx)) nx=663;	
	/*NX*/

        if (!sf_getfloat("dz",&dz)) dz=10;	
	/*DZ*/

        if (!sf_getfloat("dx",&dx)) dx=10;	
	/*DX*/

        if (!sf_getint("L",&L)) L=20;	
	/*L*/

        if (!sf_getint("nt",&nt)) nt=2001;	
	/*nt*/

        if (!sf_getfloat("dt",&dt)) dt=0.001;	
	/*dt*/

        if (!sf_getfloat("f0",&f0)) f0=20;	
	/*dominant frequency of the Ricker wavelet*/

	ntz=nz+2*L;
	ntx=nx+2*L;
	nxh=ntx/2;
	nzh=ntz/2;
	ntp=ntz*ntx;
	np=nx*nz;
	

	float t0=1/f0;			// ricker delay time
	float Omega0=20*PI*f0;	// reference frequency

	//=========================================================
	//  Flags (in this package we set flags=0 as defult)
	//  =======================================================
	//int RTMtype=3;		// RTMtype=0 for acoustic RTM
						// RTMtype=1 for viscoacoustic RTM without compensation
						// RTMtype=2 for QRTM using low-pass filtering
						// RTMtype=3 for QRTM using adaptive stabilization scheme

	int Save_Not=0;		// Save_Not=1 for save forward, reconstruction and backward wavefileds
						// Save_Not=0 for don't save forward, reconstruction and backward wavefileds

	int Sto_Rec=0;		// Sto_Rec=1 for wavefiled storage
						// Sto_Rec=0 for wavefiled reconstruction

	int Cut_Sub=0;		// Cut_Sub=1 for Cut direct wave from records
						// Cut_Sub=0 for Substract direct wave by modeling

	int Filtertype;		// Filtertype=0 for Empirical time-invaiant low-pass filtering
						// Filtertype=1 for Adaptive stabilization
						// Filtertype=2 for Non-stabilization

	int beta1, beta2;				// beta1=0 or 1 for non-disperssive or disperssive; beta2=0 or 1 for non-attenuating or attenuating
	int beta1_c, beta2_c;			// beta1_c=beta1_c for correcting phase distortion; beta2_c=-1*beta2 for compensating amplitude loss

	if(RTMtype==0)			// RTMtype=0 for acoustic RTM
	{
		Filtertype=2;		// Non-stabilization
		beta1=0;			// Non-dispersion
		beta2=0;			// Non-attenuation
		beta1_c=0;			// Non-dispersion correction
		beta2_c=0;			// Non-attenuation compensation
		sf_warning("Acoustic RTM is selected!");
	}

	else if(RTMtype==1)		// RTMtype=1 for viscoacoustic RTM without compensation
	{
		Filtertype=2;		// Non-stabilization
		beta1=1;			// dispersion
		beta2=1;			// attenuation
		beta1_c=0;			// Non-dispersion correction
		beta2_c=0;			// Non-attenuation compensation
		sf_warning("Viscocoustic RTM without compensation is selected!");
	}
	else if(RTMtype==2)		// RTMtype=2 for QRTM using low-pass filtering
	{
		Filtertype=0;		// low-pass filtering
		beta1=1;			// dispersion
		beta2=1;			// attenuation
		beta1_c=1;			// dispersion correction
		beta2_c=-1;			// attenuation compensation
		sf_warning("Q-RTM using is low-pass filtering selected!");
	}
	else if(RTMtype==3)		// RTMtype=3 for QRTM using adaptive stabilization scheme
	{
		Filtertype=1;		// adaptive stabilization scheme
		beta1=1;			// dispersion
		beta2=1;			// attenuation
		beta1_c=1;			// dispersion correction
		beta2_c=-1;			// attenuation compensation
		sf_warning("Q-RTM using adaptive stabilization is selected!");
	}
	else
	{
		sf_warning("Please define RTM type!");
	}	
											
	int Ckptstype=0;	// Ckptstype=0 for Ave-Distribution
						// Ckptstype=1 for Log-Distribution
						// Ckptstype=2 for Hyb-Distribution 

	int vp_type=0;		// vp_type=0 for homogeneous model
						// vp_type=1 for ture model
						// vp_type=2 for initial model			


	if(myid==0)
	{
		ricker_wave(ricker,nt,f0,t0,dt,3); // where parameter 3 means ricker_derivative
		sf_warning("Ricker wave is done!");
	}
	MPI_Bcast(ricker,nt,MPI_FLOAT,0,comm);	// broadcast ricker to each node


	//=========================================================
	//  Parameters of Sources and Receivers...
	//  =======================================================

  	int is, irx,irz;
	int nsid,modsr,prcs;
	int iss,eachsid,offsets;

	int ds;				// shots interval
	int s0;				// position of the first shot is L+s0
	int ns;	// total shot number is 64 

	int rnmax=0;					// maxium receiver numbers
	int nrx_obs=100; 				// 100 receivers every side of shot for backward propagation

        if (!sf_getint("ds",&ds)) ds=10;	
	/*ds*/
	
        if (!sf_getint("s0",&s0)) s0=15;	
	/* position of the first shot is L+s0*/	

        if (!sf_getint("ns",&ns)) ns=64; /*maximum: (ntx-2*L-2*s0)/ds+1*/
	/*total shot number*/


	if (myid==0)
	{
		sf_warning("type=%d",RTMtype);
		sf_warning("ds=%d, s0=%d, ns=%d",ds,s0,ns);
		sf_warning("dt=%g, nt=%d",dt,nt);
		sf_warning("f0=%g",f0);
	}




	struct Source ss[ns];			// struct pointer for source variables
	for(is=0;is<ns;is++)
	{
		ss[is].s_ix=is*ds+L+s0;	// receiver horizontal position
		ss[is].s_iz=L+5;			// shot vertical position
		ss[is].r_iz=L+5;			// receiver vertical position
		ss[is].r_n=nx; 				// one shot to all receivers
		ss[is].r_ix=(int*)malloc(sizeof(int)*ss[is].r_n);
		for(ip=0;ip<ss[is].r_n;ip++)
		{
			ss[is].r_ix[ip]=L+ip;	// shot horzitonal position
		}	
		if(rnmax<ss[is].r_n)
			rnmax=ss[is].r_n;		// maxium receiver numbers
	}
	if(myid==0)
	{
		sf_warning("The total shot number is %d",ns);
		sf_warning("The maximum trace number for source is %d",rnmax);
	}


	//=========================================================
	//  Parameters of checkpoints...
	//  =======================================================

	int check_steps;			// checkpoints interval
	int N_cp;					// total number of checkpoints
	int *t_cp;					// time point for checkpoints

	if(Ckptstype==0)	// for Ave-Distribution
	{
		check_steps=200;							// every 200 timestep has a checkpoint
		N_cp=2*(int)(nt/check_steps);				// total number of checkpoints (in pair)
		t_cp= (int *) malloc(N_cp*sizeof(int));		// checkpoints pointer
		for(int icp=0;icp<N_cp;icp++)
		{
			if(icp%2==0)
				t_cp[icp]=check_steps*(icp/2+1);	// checkpoints at even timestep
			else
				t_cp[icp]=t_cp[icp-1]+1;			// checkpoints at odd timestep
			if(myid==0)
			{
				sf_warning("checkpoints time is %d", t_cp[icp]);
			}
		}			
	}		

	else if(Ckptstype==1)	// Log-Distribution
	{
		check_steps=20;											// every 20*2^n timestep has a checkpoint
		N_cp=2*(int)((log(nt/check_steps)/log(2)+1));

		if(pow(2,N_cp/2-1)+pow(2,N_cp/2-2) < nt/check_steps)
			N_cp+=2;

		t_cp= (int *) malloc(N_cp*sizeof(int));

		for(int icp=0;icp<N_cp;icp++)
		{
			if(icp%2==0)
			{
				t_cp[icp]=check_steps*(int)(pow(2,icp/2));		// checkpoints at even timestep
				if(t_cp[icp]> nt)
					t_cp[icp]-=t_cp[icp-2]/2;
			}
			else
				t_cp[icp]=t_cp[icp-1]+1;						// checkpoints at odd timestep

			if(myid==0)
			{
				sf_warning("checkpoints time is %d", t_cp[icp]);
			}		
		}		
	}

	else if(Ckptstype==2)	// Hyb-Distribution
	{
		check_steps=200;
		float min_steps=20;
		int min_N_cp;
		min_N_cp=2*(int)((log(check_steps/min_steps)/log(2)+1));
		if(pow(2,min_N_cp/2-1)+pow(2,min_N_cp/2-2) < check_steps/min_steps)
			min_N_cp+=2;
		N_cp= min_N_cp+2*(int)(nt/check_steps);
		t_cp= (int *) malloc(N_cp*sizeof(int*));

		for(int icp=0;icp<min_N_cp;icp++)
		{
			if(icp%2==0)
			{
				t_cp[icp]=min_steps*(int)(pow(2,icp/2));
				if(t_cp[icp]> check_steps)
					t_cp[icp]-=t_cp[icp-2]/2;
			}
			else
				t_cp[icp]=t_cp[icp-1]+1;

			if(myid==0)
			{
				sf_warning("checkpoints time is %d", t_cp[icp]);
			}	
		}
		for(int icp=min_N_cp;icp<N_cp;icp++)
		{
			if(icp%2==0)
				t_cp[icp]=check_steps*((icp-min_N_cp)/2+1);
			else
				t_cp[icp]=t_cp[icp-1]+1;

			if(myid==0)
			{
				sf_warning("checkpoints time is %d", t_cp[icp]);
			}	
		}				
	}

	if(myid==0)
	{
		sf_warning("%d Checkpoints are done!", N_cp);
	}


	//=========================================================
	//  Parameters of GPU...(we assume each node has the same number of GPUs)
	//  =======================================================

	int i,GPU_N;						// GPU_N stands for the total number of GPUs per node
	getdevice(&GPU_N);					// Obtain the number of GPUs per node
	sf_warning("The available Device number is %d on %s",GPU_N,processor_name);
	struct MultiGPU plan[GPU_N];		// struct pointer for MultiGPU variables

	nsid=ns/(GPU_N*numprocs);			// shots number per GPU card
	modsr=ns%(GPU_N*numprocs);			// residual shots after average shot distribution
	prcs=modsr/GPU_N;					// which GPUs at each node will have one more shot
	if(myid<prcs)						
	{
		eachsid=nsid+1;					// if thread ID less than prcs, the corresponding GUPs have one more shot
		offsets=myid*(nsid+1)*GPU_N;	// the offset of the shots
	}
	else
	{
		eachsid=nsid;												// the rest GUPs have nsid shots
		offsets=prcs*(nsid+1)*GPU_N+(myid-prcs)*nsid*GPU_N;			// the offset of the shots (including previous shots)
	}
	

	//=========================================================
	//  Parameters of model...
	//  =======================================================

	float *vp, *Qp;						// velocity and Q models
	float *Gamma, averGamma;			// Q-related parameter Gamma and its average
	float vp_max,Qp_max;				
	float vp_min,Qp_min;
	float avervp;

	vp = (float*)malloc(sizeof(float)*ntp);
	Qp = (float*)malloc(sizeof(float)*ntp);
	Gamma = (float*)malloc(sizeof(float)*ntp);	

	float *vp0,*Qp0;
	vp0=sf_floatalloc(nx*nz);
	Qp0=sf_floatalloc(nx*nz);

	if(myid==0)
	{
		sf_floatread(vp0,nx*nz,in);
		sf_floatread(Qp0,nx*nz,inq);

		get_acc_model(vp0,Qp0,vp,Qp,ntp,ntx,ntz,L);
		
		/*fp=fopen("./output/acc_vp.dat","wb");
		for(ix=L;ix<=ntx-L-1;ix++)
		{
			for(iz=L;iz<=ntz-L-1;iz++)
			{
				fwrite(&vp[iz*ntx+ix],sizeof(float),1,fp);

			}
		}
		fclose(fp);
		fp=fopen("./output/acc_Qp.dat","wb");
		for(ix=L;ix<=ntx-L-1;ix++)
		{
			for(iz=L;iz<=ntz-L-1;iz++)
			{
				fwrite(&Qp[iz*ntx+ix],sizeof(float),1,fp);

			}
		}
		fclose(fp);*/

		vp_max=0.0;
		Qp_max=0.0;
		vp_min=5000.0;
		Qp_min=5000.0;

		for(ip=0;ip<ntp;ip++)
		{     
			if(vp[ip]>=vp_max)
			{
				vp_max=vp[ip];
			}
			if(Qp[ip]>=Qp_max)
			{
				Qp_max=fabs(Qp[ip]);
			}
			if(vp[ip]<=vp_min)
			{
				vp_min=vp[ip];
			}
			if(Qp[ip]<=Qp_min)
			{
				Qp_min=fabs(Qp[ip]);
			}
		}

		sf_warning("vp_max = %f",vp_max); 
		sf_warning("Qp_max = %f",Qp_max);
		sf_warning("vp_min = %f",vp_min); 
		sf_warning("Qp_min = %f",Qp_min);

		averGamma=0.0;
		avervp=0.0;

		for(iz=0;iz<=ntz-1;iz++)  
			for(ix=0;ix<=ntx-1;ix++)
			{
				Qp[iz*ntx+ix] /= 1;
				Gamma[iz*ntx+ix]=atan(1/Qp[iz*ntx+ix])/PI;
				averGamma+=Gamma[iz*ntx+ix]/(ntx*ntz);
				avervp+=vp[iz*ntx+ix]/(ntx*ntz);		
			}

		sf_warning("The true model is done!"); 
	}

	MPI_Bcast(vp, ntp, MPI_FLOAT, 0, comm);	
	MPI_Bcast(Qp, ntp, MPI_FLOAT, 0, comm);
	MPI_Bcast(&vp_max, 1, MPI_FLOAT, 0, comm);
	MPI_Bcast(&Qp_max, 1, MPI_FLOAT, 0, comm);
	MPI_Bcast(&vp_min, 1, MPI_FLOAT, 0, comm);
	MPI_Bcast(&Qp_min, 1, MPI_FLOAT, 0, comm);	
	MPI_Bcast(Gamma, ntp, MPI_FLOAT, 0, comm);
	MPI_Bcast(&averGamma, 1, MPI_FLOAT, 0, comm);
	MPI_Bcast(&avervp, 1, MPI_FLOAT, 0, comm);

	float *Inner_image_cor, *Inner_image_nor, *Final_image_cor, *Final_image_nor;

	Inner_image_cor=(float*)malloc(sizeof(float)*np);
	Inner_image_nor=(float*)malloc(sizeof(float)*np);
	Final_image_cor=(float*)malloc(sizeof(float)*np);
	Final_image_nor=(float*)malloc(sizeof(float)*np);
	memset(Inner_image_cor,0,np*sizeof(float));
	memset(Inner_image_nor,0,np*sizeof(float));
	memset(Final_image_cor,0,np*sizeof(float));
	memset(Final_image_nor,0,np*sizeof(float));


	// The following two functions are responsible for alloc and initialization struct varibale plan for GPU_N
	cuda_Device_malloc(ntx, ntz, ntp, nx, nz, nt, dx, dz, L, rnmax, N_cp, plan, GPU_N);
	cuda_Host_initialization(ntx, ntz, ntp, nx, nz, nt, dx, dz, L, rnmax, N_cp, plan, GPU_N);





	//=========================================================
	//  Filtering and Stabilization parameters
	//  =======================================================
	// empirical cutoff wavenumber for low-pass filtering
	float kx_cut=3.0*2*PI*f0/vp_max;
	float kz_cut=3.0*2*PI*f0/vp_max;
	float kx_cut_inv=3.0*2*PI*f0/vp_max;
	float kz_cut_inv=3.0*2*PI*f0/vp_max;
	float taper_ratio=0.2;				// taper ratio for turkey window filter 

	// parameter for adaptive satbiliztion scheme
	float sigma=2.5e-3;
	int Order=1;

	// define stabilized k-space field to verify two stabilization schemes
	float *kfilter, *kstabilization;
	kfilter = (float*)malloc(sizeof(float)*ntp);
	kstabilization = (float*)malloc(sizeof(float)*ntp);

	MPI_Barrier(comm);	// MPI barrier to ensure that all variables have been well-defined





	//=======================================================
	//  Calculate the Observed seismograms...
	//  (this process is designed for synthetic example)
	//========================================================

	for(iss=0; iss<eachsid; iss++)	// each GPU card compute eachside shots 
	{		
		is=offsets+iss*GPU_N;		// current shot index
		vp_type = 1; 				// ture model
		Save_Not=0;					// Don't save any wavefield snapshots
		sf_warning("Shot number is %d/%d",is+1,ns);
		// viscoacoustic wave equation modeling using ture model to obtain observed seismogram
		cuda_visco_PSM_2d_forward
		(
			beta1, beta2,
			nt, ntx, ntz, ntp, nx, nz, L, dx, dz, dt,
			vp, Gamma, avervp, averGamma, f0, Omega0, ricker,
			myid, is, ss, plan, GPU_N, rnmax, nrx_obs, N_cp, t_cp,
			kx_cut, kz_cut, sigma, Order, taper_ratio, kfilter, kstabilization,
			Sto_Rec, vp_type, Save_Not, Filtertype
		);

		// write observed seismogram to disk
		for(i=0;i<GPU_N;i++)
		{
			sprintf(filename,"./output/%dsource_seismogram_obs.dat",is+i+1);
			fp=fopen(filename,"wb");
			for(ix=0;ix<ss[is+i].r_n;ix++)
			{
				for(it=0;it<nt;it++)
				{
					fwrite(&plan[i].seismogram_obs[it*ss[is+i].r_n+ix],sizeof(float),1,fp);
				}
			}
			fclose(fp);
		}

		// remove direct wave from observed seismogram by mutting and obtain the residual seismogram
		if(Cut_Sub==1)
		{
			for(i=0;i<GPU_N;i++)
			{
				cut_dir(plan[i].seismogram_obs, plan[i].seismogram_rms, rnmax, nt, is, dx, dz, dt, 
					ss[is+i].r_iz, ss[is+i].s_ix, ss[is+i].s_iz, t0, vp);
				sprintf(filename,"./output/%dsource_seismogram_rms.dat",is+i+1);
				fp=fopen(filename,"wb");
				for(ix=0;ix<ss[is+i].r_n;ix++)
				{
					for(it=0;it<nt;it++)
					{
						fwrite(&plan[i].seismogram_rms[it*ss[is+i].r_n+ix],sizeof(float),1,fp);
					}
				}
				fclose(fp);
			}				
		}
	}

	if(myid==0)
	{
		sf_warning("seismogram_obs is obtained!");
	}
	if(Cut_Sub==1&&myid==0)
	{
		sf_warning("seismogram_rms is obtained!");
	}


	//=======================================================
	//  Calculate the direct and rms seismograms...
	//  (this process is designed for removing direct wave from 
	//  observed seismogram by subtracting direct wave, which is
	//  simulated by homogeneous velocity and Q model)
	//========================================================

	if(Cut_Sub==0)
	{
		if(myid==0)
		{
			get_homo_model(vp,ntp,ntx,ntz,L);
			get_homo_model(Gamma,ntp,ntx,ntz,L);
			fp=fopen("./output/homo_vp.dat","wb");
			for(ix=L;ix<=ntx-L-1;ix++)
			{
				for(iz=L;iz<=ntz-L-1;iz++)
				{
					fwrite(&vp[iz*ntx+ix],sizeof(float),1,fp);
				}
			}
			fclose(fp);
			sf_warning("The homogeneous model is done!"); 			
		}
		MPI_Bcast(vp, ntp, MPI_FLOAT, 0, comm);
		MPI_Bcast(Gamma, ntp, MPI_FLOAT, 0, comm);

		for(int iss=0; iss<eachsid; iss++)  
		{
			is=offsets+iss*GPU_N;
			vp_type = 0; // homogeneous model
			Save_Not=0;

			//// viscoacoustic wave equation modeling using homogeneous model to obtain direct seismogram
			cuda_visco_PSM_2d_forward
			(
				beta1, beta2,
				nt, ntx, ntz, ntp, nx, nz, L, dx, dz, dt,
				vp, Gamma, avervp, averGamma, f0, Omega0, ricker,
				myid, is, ss, plan, GPU_N, rnmax, nrx_obs, N_cp, t_cp,
				kx_cut, kz_cut, sigma, Order, taper_ratio, kfilter, kstabilization,
				Sto_Rec, vp_type, Save_Not, Filtertype
			);

			// write direct seismogram to disk and calculate residual seismogram 
			for(i=0;i<GPU_N;i++)
			{
				sprintf(filename,"./output/%dsource_seismogram_dir.dat",is+i+1);
				fp=fopen(filename,"wb");
				for(ix=0;ix<ss[is+i].r_n;ix++)
				{
					for(it=0;it<nt;it++)
					{
						fwrite(&plan[i].seismogram_dir[it*ss[is+i].r_n+ix],sizeof(float),1,fp);
					}
				}
				fclose(fp);

				sprintf(filename,"./output/%dsource_seismogram_obs.dat",is+i+1);
				fp=fopen(filename,"rb");
				for(ix=0;ix<ss[is+i].r_n;ix++)
				{
					for(it=0;it<nt;it++)
					{
						fread(&plan[i].seismogram_obs[it*ss[is+i].r_n+ix],sizeof(float),1,fp);
					}
				}
				fclose(fp);

				for(ip=0;ip<ss[is+i].r_n*nt;ip++)
				{
					plan[i].seismogram_rms[ip]=plan[i].seismogram_obs[ip]-plan[i].seismogram_dir[ip];
				}

				sprintf(filename,"./output/%dsource_seismogram_rms.dat",is+i+1);
				fp=fopen(filename,"wb");
				for(ix=0;ix<ss[is+i].r_n;ix++)
				{
					for(it=0;it<nt;it++)
					{
						fwrite(&plan[i].seismogram_rms[it*ss[is+i].r_n+ix],sizeof(float),1,fp);
					}
				}
				fclose(fp);
			}
		}
		if(myid==0)
		{
			sf_warning("seismogram_dir is obtained!");
			sf_warning("seismogram_rms is obtained!");
		}
	}



	//=======================================================
	//  Construct the forward wavefields and Back-propagate
	//  the RMS seismograms, Meanwhile the images are computed... 
	//========================================================

	if(myid==0)
	{
		sf_warning("====================");
		sf_warning("    RTM BEGIN");
		sf_warning("====================");

		start=clock();				// start clock
	}

	// obtain the initial model for forward and backward propagation


	if(myid==0)
	{

		// obtain the ture velocity and Q model
		get_acc_model(vp0,Qp0,vp,Qp,ntp,ntx,ntz,L);

		for(iz=0;iz<=ntz-1;iz++)  
			for(ix=0;ix<=ntx-1;ix++)
			{
				Qp[iz*ntx+ix] /= 1;
				Gamma[iz*ntx+ix]=atan(1/Qp[iz*ntx+ix])/PI;	
			}
		sf_warning("The true model is done!"); 

		// obtain the initial model
		get_ini_model(vp,ntp,ntx,ntz,20);
		fp=fopen("./output/ini_vp.dat","wb");
		for(ix=L;ix<=ntx-L-1;ix++)
		{
			for(iz=L;iz<=ntz-L-1;iz++)
			{
				fwrite(&vp[iz*ntx+ix],sizeof(float),1,fp);
			}
		}
		fclose(fp);
		sf_warning("The initial model is done!"); 
	}
	MPI_Bcast(vp, ntp, MPI_FLOAT, 0, comm);
	MPI_Bcast(Gamma, ntp, MPI_FLOAT, 0, comm);



	for(int iss=0; iss<eachsid; iss++)  
	{
		is=offsets+iss*GPU_N;
		vp_type = 2;				// initial model
		Save_Not=0;

		// compensated source wavefiled propagation using initial velocity
		cuda_visco_PSM_2d_forward
		(
			beta1_c, beta2_c,
			nt, ntx, ntz, ntp, nx, nz, L, dx, dz, dt,
			vp, Gamma, avervp, averGamma, f0, Omega0, ricker,
			myid, is, ss, plan, GPU_N, rnmax, nrx_obs, N_cp, t_cp,
			kx_cut, kz_cut, sigma, Order, taper_ratio, kfilter, kstabilization,
			Sto_Rec, vp_type, Save_Not, Filtertype
		);

		// write the synthetic seismogram to disk
		for(i=0;i<GPU_N;i++)
		{
			sprintf(filename,"./output/%dsource_seismogram_syn.dat",is+i+1);
			fp=fopen(filename,"wb");
			for(ix=0;ix<ss[is+i].r_n;ix++)
			{
				for(it=0;it<nt;it++)
				{
					fwrite(&plan[i].seismogram_syn[it*ss[is+i].r_n+ix],sizeof(float),1,fp);
				}
			}
			fclose(fp);
		}  

		// read the residual seismogram from disk
		for(i=0;i<GPU_N;i++)
		{
			sprintf(filename,"./output/%dsource_seismogram_rms.dat",is+i+1);
			fp=fopen(filename,"rb");
			for(ix=0;ix<ss[is+i].r_n;ix++)
			{
				for(it=0;it<nt;it++)
				{
					fread(&plan[i].seismogram_rms[it*ss[is+i].r_n+ix],sizeof(float),1,fp);
				}
			}
			fclose(fp);
		}

		// compensated receiver wavefiled propagation and imaging using initial velocity
		cuda_visco_PSM_2d_backward
		(
			beta1_c, beta2_c,
			nt, ntx, ntz, ntp, nx, nz, L, dx, dz, dt,
			vp, Gamma, avervp, averGamma, f0, Omega0, ricker,
			myid, is, ss, plan, GPU_N, rnmax, nrx_obs, N_cp, t_cp,
			kx_cut, kz_cut, sigma, Order, taper_ratio, kfilter, kstabilization,
			Sto_Rec, Save_Not, Filtertype
		);

		// output images
		for(i=0;i<GPU_N;i++)
		{
			// calculate normalized image
			float image_sources_max=0.0;
			for(ip=0;ip<ntp;ip++)
			{
				if(image_sources_max<fabs(plan[i].image_sources[ip]))
				{
					image_sources_max=fabs(plan[i].image_sources[ip]);
				}
			}
			for(ip=0;ip<ntp;ip++)		
			{
				plan[i].image_nor[ip]=plan[i].image_cor[ip]
						/(plan[i].image_sources[ip]+1.0e-5*image_sources_max);
			}

			// Laplace filtering
			Laplace_FD_filtering(plan[i].image_cor,ntx,ntz,dx,dz);
			Laplace_FD_filtering(plan[i].image_nor,ntx,ntz,dx,dz);

			// write all images to disk
			sprintf(filename,"./output/image_sources%d.dat",is+i+1);
			fp=fopen(filename,"wb");
			fwrite(&plan[i].image_sources[0],sizeof(float),ntp,fp);
			fclose(fp);
			sprintf(filename,"./output/image_receivers%d.dat",is+i+1);
			fp=fopen(filename,"wb");
			fwrite(&plan[i].image_receivers[0],sizeof(float),ntp,fp);
			fclose(fp);
			sprintf(filename,"./output/image_cor%d.dat",is+i+1);
			fp=fopen(filename,"wb");
			fwrite(&plan[i].image_cor[0],sizeof(float),ntp,fp);
			fclose(fp);
			sprintf(filename,"./output/image_nor%d.dat",is+i+1);
			fp=fopen(filename,"wb");
			fwrite(&plan[i].image_nor[0],sizeof(float),ntp,fp);
			fclose(fp);

			// calculate inner images at simulation domain (exclude absorbing boundary) 
			for(iz=L;iz<=ntz-L-1;iz++)
			{
				for(ix=L;ix<=ntx-L-1;ix++)
				{
					ip=iz*ntx+ix;
					ipp=(ix-L)*nz+iz-L;
					Inner_image_cor[ipp]+=plan[i].image_cor[ip]*vp[ip]*vp[ip];
					Inner_image_nor[ipp]+=plan[i].image_nor[ip]*vp[ip]*vp[ip];
				}
			}
		}
	}//end is (shotnumbers)

	// MPI barrier to ensure all shots have been calculted before stacking
	MPI_Barrier(comm);

	// stacking all shot's images
	MPI_Allreduce(Inner_image_cor,Final_image_cor,np,MPI_FLOAT,MPI_SUM,comm);
	MPI_Allreduce(Inner_image_nor,Final_image_nor,np,MPI_FLOAT,MPI_SUM,comm);

	//==========================================================
	//  Output the final images,...
	//===========================================================
	if(myid==0)
	{
		/*sprintf(filename,"./output/Final_image_cor_type%d.dat",RTMtype);
		fp=fopen(filename,"wb");
		fwrite(&Final_image_cor[0],sizeof(float),np,fp);
		fclose(fp);


		sprintf(filename,"./output/Final_image_nor_type%d.dat",RTMtype);
		fp=fopen(filename,"wb");
		fwrite(&Final_image_nor[0],sizeof(float),np,fp);
		fclose(fp);*/
		sf_putint(out,"n1",nz);
		sf_putint(out,"n2",nx);
		sf_putfloat(out,"d1",dz);
		sf_putfloat(out,"d2",dx);

		sf_putint(out2,"n1",nz);
		sf_putint(out2,"n2",nx);
		sf_putfloat(out2,"d1",dz);
		sf_putfloat(out2,"d2",dx);

		sf_floatwrite(&Final_image_cor[0], np, out);		
		sf_floatwrite(&Final_image_nor[0], np, out2);

	}

	MPI_Barrier(comm);

	if(myid==0)
	{
		sf_warning("====================");
		sf_warning("      THE END");
		sf_warning("====================");

		end=clock();
		sf_warning("The cost of the run time is %f seconds",
				(double)(end-start)/CLOCKS_PER_SEC);
	}

	//==========================================================
	//  Free the variables...
	//==========================================================

	// the following function is responsible for free struct variable plan for GPU_N
	cuda_Device_free(ntx, ntz, ntp, nx, nz, nt, dx, dz, L, rnmax, N_cp, plan, GPU_N);

	for(is=0;is<ns;is++)
	{
		free(ss[is].r_ix);
	} 
	free(ricker);
	free(vp); free(Qp); free(Gamma);
	free(t_cp);
	free(kfilter); 
	free(kstabilization);
	free(Inner_image_cor);
	free(Inner_image_nor);
	free(Final_image_cor);
	free(Final_image_nor);

	// MPI end
	MPI_Barrier(comm);
	MPI_Finalize();


	return 0;
}


