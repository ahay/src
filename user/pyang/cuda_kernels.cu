/*
  Copyright (C) 2013  Xi'an Jiaotong University (Pengliang Yang)
    Email: ypl.2100@gmail.com	

  Acknowledgement: This code is written with the help of Baoli Wang.


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

// set the positions of sources and geophones in whole domain
__global__ void cuda_set_sg(int *sxz, int sxbeg, int szbeg, int jsx, int jsz, int ns, int npml, int nnz)
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<ns) sxz[id]=nnz*(sxbeg+id*jsx+npml)+(szbeg+id*jsz+npml);
}

// generate ricker wavelet with time deley
__global__ void cuda_ricker_wavelet(float *wlt, float fm, float dt, int nt)
{
	int it=threadIdx.x+blockDim.x*blockIdx.x;
    	float tmp = PI*fm*fabsf(it*dt-1.0/fm);	//delay the wavelet to exhibit all waveform
    	tmp *=tmp;
    	if (it<nt) wlt[it]= (1.0-2.0*tmp)*expf(-tmp);	// ricker wavelet at time: t=nt*dt
}

// add==true, add (inject) the source; add==false, subtract the source
__global__ void cuda_add_source(float *p, float *source, int *Sxz, int ns, bool add)
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<ns)
	{
		if (add)	p[Sxz[id]]+=source[id];
		else 		p[Sxz[id]]-=source[id];
	}	
}


// record the seismogram at time kt
__global__ void cuda_record(float*p, float *seis_kt, int *Gxz, int ng)
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<ng) seis_kt[id]=p[Gxz[id]];
}

// mute the direct arrival according to the given velocity vmute
__global__ void cuda_mute(float *seis_kt, int gzbeg, int szbeg, int gxbeg, int sxc, int jgx, int kt, int ntd, float vmute, float dt, float dz, float dx, int ng)
{
    	int id=threadIdx.x+blockDim.x*blockIdx.x;
	float a=dx*abs(gxbeg+id*jgx-sxc);
	float b=dz*(gzbeg-szbeg);
	float t0=sqrtf(a*a+b*b)/vmute;
	int ktt=int(t0/dt)+ntd;// ntd is manually added to obtain the best muting effect.
    	if (id<ng && kt<ktt) seis_kt[id]=0.0;
}

// initialize the PML coefficients along x direction
__global__ void cuda_init_abcx(float *vel, float *bx1, float *bx2, float dx, float dz, float dt, int npml, int nnz, int nnx)
{
	// bx1: left and right PML ABC coefficients, decay p (px,pz) along x direction
	// bx2: left and right PML ABC coefficients, decay v (vx,vz) along x direction
	// only 2 blocks used horizontally, blockIdx.x=0, 1

	// id: position in top or bottom PML zone itself
	// blockIdx.x==0, left PML zone; blockIdx.x==1, right PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=blockIdx.y*npml+threadIdx.y;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);
	int ik=i1+nnz*i2;

	float Rc=1.0e-5f;
	float d=npml*MAX(dx,dz);
	float d0=-3.0f*vel[id]*logf(Rc)/d/2.0f;
	float tmp1, tmp2;

	if (i2<npml) 	// left PML zone
	{
		tmp1=(float)(npml-i2);
		tmp2=tmp1-0.5f;
	
	}
	else		// right PML zone
	{
		tmp1=i2-npml+0.5f;
		tmp2=(tmp1+0.5f);	
	}
	tmp1=tmp1/npml;
	tmp2=tmp2/npml;
	tmp1=tmp1*tmp1;
	tmp2=tmp2*tmp2;
	bx1[ik]=expf(-d0*tmp1*dt);
	bx2[ik]=expf(-d0*tmp2*dt);
}


// initialize the PML coefficients along z-axis
__global__ void cuda_init_abcz(float *vel, float *bz1, float *bz2, float dx, float dz, float dt, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone

	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int id=nnz*i2+(blockIdx.x*(nnz-npml)+threadIdx.x);
	int ik=i1+2*npml*i2;

	float Rc=1.0e-5f;
	float d=npml*MAX(dx,dz);
	float d0=-3.0f*vel[id]*logf(Rc)/d/2.0f;
	float tmp1, tmp2;

	if (i1<npml) 	// top PML zone
	{	
		tmp1=(float)(npml-i1);
		tmp2=tmp1-0.5f;	
	}
	else		// bottom PML zone
	{
		tmp1=i1-npml+0.5f;
		tmp2=(float)(tmp1+0.5f);	
	}
	tmp1=tmp1/npml;
	tmp2=tmp2/npml;
	tmp1=tmp1*tmp1;
	tmp2=tmp2*tmp2;
	bz1[ik]=expf(-d0*tmp1*dt);
	bz2[ik]=expf(-d0*tmp2*dt);
}


__global__ void cuda_init_bell(float *bell)
{
	int i1=threadIdx.x;
	int i2=threadIdx.y;
	int id=i1+i2*(2*nbell+1);
	float s = 0.5*nbell;
	bell[id]=expf(-((i1-nbell)*(i1-nbell)+(i2-nbell)*(i2-nbell))/s);
}

// inject Gaussian bell smoothed wavelet
// lauch configuration: <<<dim3(ns,1), dim3(2*nbell+1,2*nbell+1)>>>
// add==true, add (inject) the wavelet; add==false, subtract the wavelet
__global__ void cuda_add_bellwlt(float *p, float *bell, float *wlt, int *Sxz, int ns, int npml, int nnz, int nnx, bool add)
{
	int i1=threadIdx.x;
	int i2=threadIdx.y;
	int is=blockIdx.x;// source wavelet index
	
    	if (is<ns)
	{
		if (add)	p[Sxz[is]+(i1-nbell)+(i2-nbell)*nnz]+=bell[i1+i2*(2*nbell+1)]*wlt[is];
		else 		p[Sxz[is]+(i1-nbell)+(i2-nbell)*nnz]-=bell[i1+i2*(2*nbell+1)]*wlt[is];
	}	
}

//================================================ NJ=2 ========================================================
__global__ void cuda_forward_v_2(float *p, float *vx, float *vz, float _dx, float _dz, int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+i2*nnz;

	__shared__ float s_p[Block_Size1+1][Block_Size2+1];
    	s_p[threadIdx.x][threadIdx.y]=p[id];
	if (threadIdx.x>blockDim.x-2)
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+1][threadIdx.y]=p[id+1];
		else				s_p[threadIdx.x+1][threadIdx.y]=0.0f;
	} 
	if (threadIdx.y>blockDim.y-2)
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x][threadIdx.y+1]=p[id+nnz];
		else				s_p[threadIdx.x][threadIdx.y+1]=0.0f;
	}
	__syncthreads();

	float diff1=(s_p[threadIdx.x+1][threadIdx.y]-s_p[threadIdx.x][threadIdx.y]);// .x-->1st dim--> i1
	float diff2=(s_p[threadIdx.x][threadIdx.y+1]-s_p[threadIdx.x][threadIdx.y]);// .y-->2nd dim--> i2
	vz[id]=_dz*diff1;
	vx[id]=_dx*diff2;
}

__global__ void cuda_PML_vz_2(float *p, float *convpz, float *bz, float *vz, float _dz,	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0,1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=(blockIdx.x*(nnz-npml)+threadIdx.x)+nnz*i2;

	__shared__ float s_p[33][Block_Size2];	
    	s_p[threadIdx.x][threadIdx.y]=p[id]; 
	if (threadIdx.x>30)
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+1][threadIdx.y]=p[id+1];
		else				s_p[threadIdx.x+1][threadIdx.y]=0.0f;
	}
	__syncthreads();

	float diff1=(s_p[threadIdx.x+1][threadIdx.y]-s_p[threadIdx.x][threadIdx.y]);
	convpz[ik]=bz[ik]*convpz[ik]+(bz[ik]-1.0f)*_dz*diff1;	
	vz[id]+=convpz[ik];
}

__global__ void cuda_PML_vx_2(float *p,  float *convpx, float *bx, float *vx, float _dx, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+nnz*i2;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+	threadIdx.y);

	__shared__ float s_p[Block_Size1][33];
    	s_p[threadIdx.x][threadIdx.y]=p[id]; 
	if (threadIdx.y>30)
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x][threadIdx.y+1]=p[id+nnz];
		else				s_p[threadIdx.x][threadIdx.y+1]=0.0f;
	}
	__syncthreads();

	float diff2=(s_p[threadIdx.x][threadIdx.y+1]-s_p[threadIdx.x][threadIdx.y]);
	convpx[ik]=bx[ik]*convpx[ik]+(bx[ik]-1.0f)*_dx*diff2;	
	vx[id]+=convpx[ik];
}

__global__ void cuda_forward_up_2(float *up, float *vx, float *vz, float _dx, float _dz, int npml, int nnz, int nnx)
{
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int id=i1+i2*nnz;

	__shared__ float s_v1[Block_Size1+1][Block_Size2];
	__shared__ float s_v2[Block_Size1][Block_Size2+1];
    	s_v1[threadIdx.x+1][threadIdx.y]=vz[id]; 
    	s_v2[threadIdx.x][threadIdx.y+1]=vx[id]; 
	if (threadIdx.x<1)
	{
		if (blockIdx.x) 		s_v1[threadIdx.x][threadIdx.y]=vz[id-1];
		else				s_v1[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y<1)
	{
		if (blockIdx.y)			s_v2[threadIdx.x][threadIdx.y]=vx[id-nnz];
		else 				s_v2[threadIdx.x][threadIdx.y]=0.0f;
	}
	__syncthreads();

	float diff1=(s_v1[threadIdx.x+1][threadIdx.y]-s_v1[threadIdx.x][threadIdx.y]);
	float diff2=(s_v2[threadIdx.x][threadIdx.y+1]-s_v2[threadIdx.x][threadIdx.y]);
	up[id]=_dz*diff1+_dx*diff2;
}



__global__ void cuda_PML_upz_2(float *up,  float *convvz, float *bz, float *vz, float _dz, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=(blockIdx.x*(nnz-npml)+	threadIdx.x)+	nnz*	i2;

	__shared__ float s_v1[33][Block_Size2];
    	s_v1[threadIdx.x+1][threadIdx.y]=vz[id]; 
	if (threadIdx.x<1)
	{
		if (blockIdx.x)			s_v1[threadIdx.x][threadIdx.y]=vz[id-1];
		else 				s_v1[threadIdx.x][threadIdx.y]=0.0f;
	}
	__syncthreads();

	float diff1=(s_v1[threadIdx.x+1][threadIdx.y]-s_v1[threadIdx.x][threadIdx.y]);
	convvz[ik]=bz[ik]*convvz[ik]+(bz[ik]-1.0f)*_dz*diff1;
	up[id]+=convvz[ik];
}

__global__ void cuda_PML_upx_2(float *up,  float *convvx, float *bx, float *vx, float _dx, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+nnz*i2;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);

	__shared__ float s_v2[Block_Size1][33];
    	s_v2[threadIdx.x][threadIdx.y+1]=vx[id]; 
	if (threadIdx.y<1)
	{
		if (blockIdx.y) 		s_v2[threadIdx.x][threadIdx.y]=vx[id-nnz];
		else				s_v2[threadIdx.x][threadIdx.y]=0.0f;
	}	
	__syncthreads();

	float diff2=(s_v2[threadIdx.x][threadIdx.y+1]-s_v2[threadIdx.x][threadIdx.y]);
	convvx[ik]=bx[ik]*convvx[ik]+(bx[ik]-1.0f)*_dx*diff2;
	up[id]+=convvx[ik];
}


//================================================ NJ=4 ========================================================
__global__ void cuda_forward_v_4(float *p, float *vx, float *vz, float _dx, float _dz, int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+i2*nnz;

	__shared__ float s_p[Block_Size1+3][Block_Size2+3];
    	s_p[threadIdx.x+1][threadIdx.y+1]=p[id]; 
	if (threadIdx.x<1)
	{
		if (blockIdx.x)			s_p[threadIdx.x][threadIdx.y+1]=p[id-1];
		else 				s_p[threadIdx.x][threadIdx.y+1]=0.0f;
	}
	if (threadIdx.x>blockDim.x-3)
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+3][threadIdx.y+1]=p[id+2];
		else				s_p[threadIdx.x+3][threadIdx.y+1]=0.0f;
	}
	if (threadIdx.y<1)
	{
		if (blockIdx.y) 		s_p[threadIdx.x+1][threadIdx.y]=p[id-nnz];
		else				s_p[threadIdx.x+1][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>blockDim.y-3)
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x+1][threadIdx.y+3]=p[id+2*nnz];
		else				s_p[threadIdx.x+1][threadIdx.y+3]=0.0f;
	}
	__syncthreads();   
  
	float diff1=1.125f*(s_p[threadIdx.x+2][threadIdx.y+1]-s_p[threadIdx.x+1][threadIdx.y+1])
			-0.041666666666667f*(s_p[threadIdx.x+3][threadIdx.y+1]-s_p[threadIdx.x][threadIdx.y+1]);
	float diff2=1.125f*(s_p[threadIdx.x+1][threadIdx.y+2]-s_p[threadIdx.x+1][threadIdx.y+1])
			-0.041666666666667f*(s_p[threadIdx.x+1][threadIdx.y+3]-s_p[threadIdx.x+1][threadIdx.y]);
	vz[id]=_dz*diff1;
	vx[id]=_dx*diff2;
}

__global__ void cuda_PML_vz_4(float *p,  float *convpz, float *bz, float *vz, float _dz, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=blockIdx.x*(nnz-npml)+threadIdx.x+nnz*i2;

	__shared__ float s_p[35][Block_Size2];
    	s_p[threadIdx.x+1][threadIdx.y]=p[id]; 
	if (threadIdx.x<1)
	{
		if (blockIdx.x)			s_p[threadIdx.x][threadIdx.y]=p[id-1];
		else 				s_p[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>29)
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+3][threadIdx.y]=p[id+2];
		else				s_p[threadIdx.x+3][threadIdx.y]=0.0f;
	}
	__syncthreads();


	float diff1=1.125f*(s_p[threadIdx.x+2][threadIdx.y]-s_p[threadIdx.x+1][threadIdx.y])
			-0.041666666666667f*(s_p[threadIdx.x+3][threadIdx.y]-s_p[threadIdx.x][threadIdx.y]);
	convpz[ik]=bz[ik]*convpz[ik]+(bz[ik]-1.0f)*_dz*diff1;	
	vz[id]+=convpz[ik];
}
__global__ void cuda_PML_vx_4(float *p,  float *convpx, float *bx, float *vx, float _dx, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+i2*nnz;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);

	__shared__ float s_p[Block_Size1][35];// npml+3=35; Block_SizeX=32; Block_SizeY=8;
    	s_p[threadIdx.x][threadIdx.y+1]=p[id]; 
	if (threadIdx.y<1)
	{
		if (blockIdx.y) 		s_p[threadIdx.x][threadIdx.y]=p[id-nnz];
		else				s_p[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>29)
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x][threadIdx.y+3]=p[id+2*nnz];
		else				s_p[threadIdx.x][threadIdx.y+3]=0.0f;
	}
	__syncthreads();

	float diff2=1.125f*(s_p[threadIdx.x][threadIdx.y+2]-s_p[threadIdx.x][threadIdx.y+1])
			-0.041666666666667f*(s_p[threadIdx.x][threadIdx.y+3]-s_p[threadIdx.x][threadIdx.y]);
	convpx[ik]=bx[ik]*convpx[ik]+(bx[ik]-1.0f)*_dx*diff2;	
	vx[id]+=convpx[ik];
}

__global__ void cuda_forward_up_4(float *up, float *vx, float *vz, float _dx, float _dz, int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+i2*nnz;

	__shared__ float s_vx[Block_Size1][Block_Size2+3];
	__shared__ float s_vz[Block_Size1+3][Block_Size2];
    	s_vx[threadIdx.x][threadIdx.y+2]=vx[id]; 
    	s_vz[threadIdx.x+2][threadIdx.y]=vz[id]; 

	if (threadIdx.x<2)
	{
		if (blockIdx.x)			s_vz[threadIdx.x][threadIdx.y]=vz[id-2];
		else 				s_vz[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>blockDim.x-2)
	{
		if (blockIdx.x<gridDim.x-1)	s_vz[threadIdx.x+3][threadIdx.y]=vz[id+1];
		else				s_vz[threadIdx.x+3][threadIdx.y]=0.0f;
	}
	if (threadIdx.y<2)
	{
		if (blockIdx.y) 		s_vx[threadIdx.x][threadIdx.y]=vx[id-2*nnz];
		else				s_vx[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>blockDim.y-2)
	{
		if (blockIdx.y<gridDim.y-1)	s_vx[threadIdx.x][threadIdx.y+3]=vx[id+nnz];
		else				s_vx[threadIdx.x][threadIdx.y+3]=0.0f;
	}
	__syncthreads();

	float diff2=1.125f*(s_vx[threadIdx.x][threadIdx.y+2]-s_vx[threadIdx.x][threadIdx.y+1])
			 -0.041666666666667f*(s_vx[threadIdx.x][threadIdx.y+3]-s_vx[threadIdx.x][threadIdx.y]);
	float diff1=1.125f*(s_vz[threadIdx.x+2][threadIdx.y]-s_vz[threadIdx.x+1][threadIdx.y])+
			 -0.041666666666667f*(s_vz[threadIdx.x+3][threadIdx.y]-s_vz[threadIdx.x][threadIdx.y]);
	up[id]=_dz*diff1+_dx*diff2;
}

__global__ void cuda_PML_upz_4(float *up,  float *convvz, float *bz, float *vz, float _dz, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=blockIdx.x*(nnz-npml)+threadIdx.x+nnz*i2;

	__shared__ float s_vz[35][Block_Size2];
    	s_vz[threadIdx.x+2][threadIdx.y]=vz[id]; 
	if (threadIdx.x<2)
	{
		if (blockIdx.x)			s_vz[threadIdx.x][threadIdx.y]=vz[id-2];
		else 				s_vz[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>30)
	{
		if (blockIdx.x<gridDim.x-1)	s_vz[threadIdx.x+3][threadIdx.y]=vz[id+1];
		else				s_vz[threadIdx.x+3][threadIdx.y]=0.0f;
	}
	__syncthreads();

	float diff1=1.125f*(s_vz[threadIdx.x+2][threadIdx.y]-s_vz[threadIdx.x+1][threadIdx.y])+
			-0.041666666666667f*(s_vz[threadIdx.x+3][threadIdx.y]-s_vz[threadIdx.x][threadIdx.y]);
	convvz[ik]=bz[ik]*convvz[ik]+(bz[ik]-1.0f)*_dz*diff1;
	up[id]+=convvz[ik];
}
__global__ void cuda_PML_upx_4(float *up, float *convvx, float *bx, float *vx, float _dx, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+i2*nnz;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);

	__shared__ float s_vx[Block_Size1][35];
    	s_vx[threadIdx.x][threadIdx.y+2]=vx[id]; 
	if (threadIdx.y<2)
	{
		if (blockIdx.y) 		s_vx[threadIdx.x][threadIdx.y]=vx[id-2*nnz];
		else				s_vx[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>30)
	{
		if (blockIdx.y<gridDim.y-1)	s_vx[threadIdx.x][threadIdx.y+3]=vx[id+nnz];
		else				s_vx[threadIdx.x][threadIdx.y+3]=0.0f;
	}
	__syncthreads();

	float diff2=1.125f*(s_vx[threadIdx.x][threadIdx.y+2]-s_vx[threadIdx.x][threadIdx.y+1])
			-0.041666666666667f*(s_vx[threadIdx.x][threadIdx.y+3]-s_vx[threadIdx.x][threadIdx.y]);
	convvx[ik]=bx[ik]*convvx[ik]+(bx[ik]-1.0f)*_dx*diff2;
	up[id]+=convvx[ik];
}


//================================================ NJ=6 ========================================================
__global__ void cuda_forward_v_6(float *p, float *vx, float *vz, float _dx, float _dz, int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+i2*nnz;

	__shared__ float s_p[Block_Size1+5][Block_Size2+5];
    	s_p[threadIdx.x+2][threadIdx.y+2]=p[id]; 
	if (threadIdx.x<2)
	{
		if (blockIdx.x)			s_p[threadIdx.x][threadIdx.y+2]=p[id-2];
		else 				s_p[threadIdx.x][threadIdx.y+2]=0.0f;
	}
	if (threadIdx.x>blockDim.x-4)
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+5][threadIdx.y+2]=p[id+3];
		else				s_p[threadIdx.x+5][threadIdx.y+2]=0.0f;
	}
	if (threadIdx.y<2)
	{
		if (blockIdx.y) 		s_p[threadIdx.x+2][threadIdx.y]=p[id-2*nnz];
		else				s_p[threadIdx.x+2][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>blockDim.y-4)
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x+2][threadIdx.y+5]=p[id+3*nnz];
		else				s_p[threadIdx.x+2][threadIdx.y+5]=0.0f;
	}
	__syncthreads();

	float diff1=1.171875f*(s_p[threadIdx.x+3][threadIdx.y+2]-s_p[threadIdx.x+2][threadIdx.y+2])
			-0.065104166666667f*(s_p[threadIdx.x+4][threadIdx.y+2]-s_p[threadIdx.x+1][threadIdx.y+2])
			+0.0046875f*(s_p[threadIdx.x+5][threadIdx.y+2]-s_p[threadIdx.x][threadIdx.y+2]);
	float diff2=1.171875f*(s_p[threadIdx.x+2][threadIdx.y+3]-s_p[threadIdx.x+2][threadIdx.y+2])
			-0.065104166666667f*(s_p[threadIdx.x+2][threadIdx.y+4]-s_p[threadIdx.x+2][threadIdx.y+1])
			+0.0046875f*(s_p[threadIdx.x+2][threadIdx.y+5]-s_p[threadIdx.x+2][threadIdx.y]);
	vz[id]=_dz*diff1;
	vx[id]=_dx*diff2;
}

__global__ void cuda_PML_vz_6(float *p,  float *convpz, float *bz, float *vz, float _dz, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=blockIdx.x*(nnz-npml)+threadIdx.x+nnz*i2;

	__shared__ float s_p[37][Block_Size2];
    	s_p[threadIdx.x+2][threadIdx.y]=p[id]; 
	if (threadIdx.x<2)
	{
		if (blockIdx.x)			s_p[threadIdx.x][threadIdx.y]=p[id-2];
		else 				s_p[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>28)
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+5][threadIdx.y]=p[id+3];
		else				s_p[threadIdx.x+5][threadIdx.y]=0.0f;
	}
	__syncthreads();

	float diff1=1.171875f*(s_p[threadIdx.x+3][threadIdx.y]-s_p[threadIdx.x+2][threadIdx.y])
			-0.065104166666667f*(s_p[threadIdx.x+4][threadIdx.y]-s_p[threadIdx.x+1][threadIdx.y])
			+0.0046875f*(s_p[threadIdx.x+5][threadIdx.y]-s_p[threadIdx.x][threadIdx.y]);
	convpz[ik]=bz[ik]*convpz[ik]+(bz[ik]-1.0f)*_dz*diff1;	
	vz[id]+=convpz[ik];
}
__global__ void cuda_PML_vx_6(float *p,  float *convpx, float *bx, float *vx, float _dx, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+i2*nnz;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);

	__shared__ float s_p[Block_Size1][37];
    	s_p[threadIdx.x][threadIdx.y+2]=p[id]; 
	if (threadIdx.y<2)
	{
		if (blockIdx.y) 		s_p[threadIdx.x][threadIdx.y]=p[id-2*nnz];
		else				s_p[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>28)
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x][threadIdx.y+5]=p[id+3*nnz];
		else				s_p[threadIdx.x][threadIdx.y+5]=0.0f;
	}
	__syncthreads();

	float diff2=1.171875f*(s_p[threadIdx.x][threadIdx.y+3]-s_p[threadIdx.x][threadIdx.y+2])
			-0.065104166666667f*(s_p[threadIdx.x][threadIdx.y+4]-s_p[threadIdx.x][threadIdx.y+1])
			+0.0046875f*(s_p[threadIdx.x][threadIdx.y+5]-s_p[threadIdx.x][threadIdx.y]);
	convpx[ik]=bx[ik]*convpx[ik]+(bx[ik]-1.0f)*_dx*diff2;	
	vx[id]+=convpx[ik];
}
__global__ void cuda_forward_up_6(float *up, float *vx, float *vz, float _dx, float _dz, int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+i2*nnz;

	__shared__ float s_vx[Block_Size1][Block_Size2+5];
	__shared__ float s_vz[Block_Size1+5][Block_Size2];
    	s_vx[threadIdx.x][threadIdx.y+3]=vx[id]; 
    	s_vz[threadIdx.x+3][threadIdx.y]=vz[id]; 

	if (threadIdx.x<3)
	{
		if (blockIdx.x)			s_vz[threadIdx.x][threadIdx.y]=vz[id-3];
		else 				s_vz[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>blockDim.x-3)
	{
		if (blockIdx.x<gridDim.x-1)	s_vz[threadIdx.x+5][threadIdx.y]=vz[id+2];
		else				s_vz[threadIdx.x+5][threadIdx.y]=0.0f;
	}
	if (threadIdx.y<3)
	{
		if (blockIdx.y) 		s_vx[threadIdx.x][threadIdx.y]=vx[id-3*nnz];
		else				s_vx[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>blockDim.y-3)
	{
		if (blockIdx.y<gridDim.y-1)	s_vx[threadIdx.x][threadIdx.y+5]=vx[id+2*nnz];
		else				s_vx[threadIdx.x][threadIdx.y+5]=0.0f;
	}
	__syncthreads();


	float diff2=1.171875f*(s_vx[threadIdx.x][threadIdx.y+3]-s_vx[threadIdx.x][threadIdx.y+2])
			-0.065104166666667f*(s_vx[threadIdx.x][threadIdx.y+4]-s_vx[threadIdx.x][threadIdx.y+1])+
			0.0046875f*(s_vx[threadIdx.x][threadIdx.y+5]-s_vx[threadIdx.x][threadIdx.y]);
	float diff1=1.171875f*(s_vz[threadIdx.x+3][threadIdx.y]-s_vz[threadIdx.x+2][threadIdx.y])+
			-0.065104166666667f*(s_vz[threadIdx.x+4][threadIdx.y]-s_vz[threadIdx.x+1][threadIdx.y])+
			0.0046875f*(s_vz[threadIdx.x+5][threadIdx.y]-s_vz[threadIdx.x][threadIdx.y]);
	up[id]=(_dz*diff1+_dx*diff2);
}

__global__ void cuda_PML_upz_6(float *up,  float *convvz, float *bz, float *vz, float _dz, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=blockIdx.x*(nnz-npml)+threadIdx.x+nnz*i2;

	__shared__ float s_vz[37][Block_Size2];// npml+5=37; Block_SizeX=32; Block_SizeY=8;
    	s_vz[threadIdx.x+3][threadIdx.y]=vz[id]; 
	if (threadIdx.x<3)
	{
		if (blockIdx.x)			s_vz[threadIdx.x][threadIdx.y]=vz[id-3];
		else 				s_vz[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>29)
	{
		if (blockIdx.x<gridDim.x-1)	s_vz[threadIdx.x+5][threadIdx.y]=vz[id+2];
		else				s_vz[threadIdx.x+5][threadIdx.y]=0.0f;
	}
	__syncthreads();

	float diff1=1.171875f*(s_vz[threadIdx.x+3][threadIdx.y]-s_vz[threadIdx.x+2][threadIdx.y])
			-0.065104166666667f*(s_vz[threadIdx.x+4][threadIdx.y]-s_vz[threadIdx.x+1][threadIdx.y])
			+0.0046875f*(s_vz[threadIdx.x+5][threadIdx.y]-s_vz[threadIdx.x][threadIdx.y]);
	convvz[ik]=bz[ik]*convvz[ik]+(bz[ik]-1.0f)*_dz*diff1;
	up[id]+=convvz[ik];
}
__global__ void cuda_PML_upx_6(float *up, float *convvx, float *bx, float *vx, float _dx, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+i2*nnz;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);

	__shared__ float s_vx[Block_Size1][37];
    	s_vx[threadIdx.x][threadIdx.y+3]=vx[id]; 
	if (threadIdx.y<3)
	{
		if (blockIdx.y) 		s_vx[threadIdx.x][threadIdx.y]=vx[id-3*nnz];
		else				s_vx[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>29)
	{
		if (blockIdx.y<gridDim.y-1)	s_vx[threadIdx.x][threadIdx.y+5]=vx[id+2*nnz];
		else				s_vx[threadIdx.x][threadIdx.y+5]=0.0f;
	}
	__syncthreads();

	float diff2=1.171875f*(s_vx[threadIdx.x][threadIdx.y+3]-s_vx[threadIdx.x][threadIdx.y+2])
			-0.065104166666667f*(s_vx[threadIdx.x][threadIdx.y+4]-s_vx[threadIdx.x][threadIdx.y+1])
			+0.0046875f*(s_vx[threadIdx.x][threadIdx.y+5]-s_vx[threadIdx.x][threadIdx.y]);
	convvx[ik]=bx[ik]*convvx[ik]+(bx[ik]-1.0f)*_dx*diff2;
	up[id]+=convvx[ik];
}


//================================================ NJ=8 ========================================================
__global__ void cuda_forward_v_8(float *p, float *vx, float *vz, float _dx, float _dz, int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+nnz*i2;

	__shared__ float s_p[Block_Size1+7][Block_Size2+7];
    	s_p[threadIdx.x+3][threadIdx.y+3]=p[id]; 
	if (threadIdx.x<3)
	{
		if (blockIdx.x)			s_p[threadIdx.x][threadIdx.y+3]=p[id-3];
		else 				s_p[threadIdx.x][threadIdx.y+3]=0.0f;
	}
	if (threadIdx.x>blockDim.x-5)
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+7][threadIdx.y+3]=p[id+4];
		else				s_p[threadIdx.x+7][threadIdx.y+3]=0.0f;
	}
	if (threadIdx.y<3)
	{
		if (blockIdx.y) 		s_p[threadIdx.x+3][threadIdx.y]=p[id-3*nnz];
		else				s_p[threadIdx.x+3][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>blockDim.y-5)
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x+3][threadIdx.y+7]=p[id+4*nnz];
		else				s_p[threadIdx.x+3][threadIdx.y+7]=0.0f;
	}
	__syncthreads();

	float diff1=1.1962890625000f*(s_p[threadIdx.x+4][threadIdx.y+3]-s_p[threadIdx.x+3][threadIdx.y+3])
			-0.0797526041667f*(s_p[threadIdx.x+5][threadIdx.y+3]-s_p[threadIdx.x+2][threadIdx.y+3])
			+0.0095703125000f*(s_p[threadIdx.x+6][threadIdx.y+3]-s_p[threadIdx.x+1][threadIdx.y+3])
			-0.0006975446429f*(s_p[threadIdx.x+7][threadIdx.y+3]-s_p[threadIdx.x][threadIdx.y+3]);
	float diff2=1.1962890625000f*(s_p[threadIdx.x+3][threadIdx.y+4]-s_p[threadIdx.x+3][threadIdx.y+3])
			-0.0797526041667f*(s_p[threadIdx.x+3][threadIdx.y+5]-s_p[threadIdx.x+3][threadIdx.y+2])
			+0.0095703125000f*(s_p[threadIdx.x+3][threadIdx.y+6]-s_p[threadIdx.x+3][threadIdx.y+1])
			-0.0006975446429f*(s_p[threadIdx.x+3][threadIdx.y+7]-s_p[threadIdx.x+3][threadIdx.y]);

	vz[id]=_dz*diff1;
	vx[id]=_dx*diff2;
}


__global__ void cuda_PML_vz_8(float *p, float *convpz, float *bz, float *vz, float _dz, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=blockIdx.x*(nnz-npml)+threadIdx.x+nnz*i2;

	__shared__ float s_p[39][Block_Size2];
    	s_p[threadIdx.x+3][threadIdx.y]=p[id]; 
	if (threadIdx.x<3)
	{
		if (blockIdx.x)			s_p[threadIdx.x][threadIdx.y]=p[id-3];
		else 				s_p[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>27)
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+7][threadIdx.y]=p[id+4];
		else				s_p[threadIdx.x+7][threadIdx.y]=0.0f;
	}
	__syncthreads();

	float diff1=1.1962890625000f*(s_p[threadIdx.x+4][threadIdx.y]-s_p[threadIdx.x+3][threadIdx.y])
			-0.0797526041667f*(s_p[threadIdx.x+5][threadIdx.y]-s_p[threadIdx.x+2][threadIdx.y])
			+0.0095703125000f*(s_p[threadIdx.x+6][threadIdx.y]-s_p[threadIdx.x+1][threadIdx.y])
			-0.0006975446429f*(s_p[threadIdx.x+7][threadIdx.y]-s_p[threadIdx.x][threadIdx.y]);
	convpz[ik]=bz[ik]*convpz[ik]+(bz[ik]-1.0f)*_dz*diff1;	
	vz[id]+=convpz[ik];
}
__global__ void cuda_PML_vx_8(float *p, float *convpx, float *bx, float *vx, float _dx, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+i2*nnz;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);

	__shared__ float s_p[Block_Size1][39];// npml+7=39; Block_SizeX=32; Block_SizeY=8;
    	s_p[threadIdx.x][threadIdx.y+3]=p[id]; 
	if (threadIdx.y<3)
	{
		if (blockIdx.y) 		s_p[threadIdx.x][threadIdx.y]=p[id-3*nnz];
		else				s_p[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>27)//32-4
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x][threadIdx.y+7]=p[id+4*nnz];
		else				s_p[threadIdx.x][threadIdx.y+7]=0.0f;
	}
	__syncthreads();

	float diff2=1.1962890625000f*(s_p[threadIdx.x][threadIdx.y+4]-s_p[threadIdx.x][threadIdx.y+3])
			-0.0797526041667f*(s_p[threadIdx.x][threadIdx.y+5]-s_p[threadIdx.x][threadIdx.y+2])+
			0.0095703125000f*(s_p[threadIdx.x][threadIdx.y+6]-s_p[threadIdx.x][threadIdx.y+1])+
			-0.0006975446429f*(s_p[threadIdx.x][threadIdx.y+7]-s_p[threadIdx.x][threadIdx.y]);
	convpx[ik]=bx[ik]*convpx[ik]+(bx[ik]-1.0f)*_dx*diff2;	
	vx[id]+=convpx[ik];
}

__global__ void cuda_forward_up_8(float *up, float *vx, float *vz, float _dx, float _dz, int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+i2*nnz;

	__shared__ float s_vx[Block_Size1][Block_Size2+7];
	__shared__ float s_vz[Block_Size1+7][Block_Size2];
    	s_vx[threadIdx.x][threadIdx.y+4]=vx[id]; 
    	s_vz[threadIdx.x+4][threadIdx.y]=vz[id]; 

	if (threadIdx.x<4)
	{
		if (blockIdx.x)			s_vz[threadIdx.x][threadIdx.y]=vz[id-4];
		else 				s_vz[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>blockDim.x-4)
	{
		if (blockIdx.x<gridDim.x-1)	s_vz[threadIdx.x+7][threadIdx.y]=vz[id+3];
		else				s_vz[threadIdx.x+7][threadIdx.y]=0.0f;
	}
	if (threadIdx.y<4)
	{
		if (blockIdx.y) 		s_vx[threadIdx.x][threadIdx.y]=vx[id-4*nnz];
		else				s_vx[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>blockDim.y-4)
	{
		if (blockIdx.y<gridDim.y-1)	s_vx[threadIdx.x][threadIdx.y+7]=vx[id+3*nnz];
		else				s_vx[threadIdx.x][threadIdx.y+7]=0.0f;
	}
	__syncthreads();


	float diff2=1.1962890625000f*(s_vx[threadIdx.x][threadIdx.y+4]-s_vx[threadIdx.x][threadIdx.y+3])
			-0.0797526041667f*(s_vx[threadIdx.x][threadIdx.y+5]-s_vx[threadIdx.x][threadIdx.y+2])+
			0.0095703125000f*(s_vx[threadIdx.x][threadIdx.y+6]-s_vx[threadIdx.x][threadIdx.y+1])+
			-0.0006975446429f*(s_vx[threadIdx.x][threadIdx.y+7]-s_vx[threadIdx.x][threadIdx.y]);
	float diff1=1.1962890625000f*(s_vz[threadIdx.x+4][threadIdx.y]-s_vz[threadIdx.x+3][threadIdx.y])+
			-0.0797526041667f*(s_vz[threadIdx.x+5][threadIdx.y]-s_vz[threadIdx.x+2][threadIdx.y])+
			0.0095703125000f*(s_vz[threadIdx.x+6][threadIdx.y]-s_vz[threadIdx.x+1][threadIdx.y])+
			-0.0006975446429f*(s_vz[threadIdx.x+7][threadIdx.y]-s_vz[threadIdx.x][threadIdx.y]);
	up[id]=_dz*diff1+_dx*diff2;
}
__global__ void cuda_PML_upz_8(float *up, float *convvz, float *bz, float *vz, float _dz, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=blockIdx.x*(nnz-npml)+threadIdx.x+nnz*i2;

	__shared__ float s_vz[39][Block_Size2];
    	s_vz[threadIdx.x+4][threadIdx.y]=vz[id]; 
	if (threadIdx.x<4)
	{
		if (blockIdx.x)			s_vz[threadIdx.x][threadIdx.y]=vz[id-4];
		else 				s_vz[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>28)
	{
		if (blockIdx.x<gridDim.x-1)	s_vz[threadIdx.x+7][threadIdx.y]=vz[id+3];
		else				s_vz[threadIdx.x+7][threadIdx.y]=0.0f;
	}
	__syncthreads();

	float diff1=1.1962890625000f*(s_vz[threadIdx.x+4][threadIdx.y]-s_vz[threadIdx.x+3][threadIdx.y])
			-0.0797526041667f*(s_vz[threadIdx.x+5][threadIdx.y]-s_vz[threadIdx.x+2][threadIdx.y])
			+0.0095703125000f*(s_vz[threadIdx.x+6][threadIdx.y]-s_vz[threadIdx.x+1][threadIdx.y])
			-0.0006975446429f*(s_vz[threadIdx.x+7][threadIdx.y]-s_vz[threadIdx.x][threadIdx.y]);
	convvz[ik]=bz[ik]*convvz[ik]+(bz[ik]-1.0f)*_dz*diff1;
	up[id]+=convvz[ik];
}
__global__ void cuda_PML_upx_8(float *up, float *convvx, float *bx, float *vx, float _dx, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+i2*nnz;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);

	__shared__ float s_vx[Block_Size1][39];
    	s_vx[threadIdx.x][threadIdx.y+4]=vx[id]; 
	if (threadIdx.y<4)
	{
		if (blockIdx.y) 		s_vx[threadIdx.x][threadIdx.y]=vx[id-4*nnz];
		else				s_vx[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>28)
	{
		if (blockIdx.y<gridDim.y-1)	s_vx[threadIdx.x][threadIdx.y+7]=vx[id+3*nnz];
		else				s_vx[threadIdx.x][threadIdx.y+7]=0.0f;
	}
	__syncthreads();

	float diff2=1.1962890625000f*(s_vx[threadIdx.x][threadIdx.y+4]-s_vx[threadIdx.x][threadIdx.y+3])
			-0.0797526041667f*(s_vx[threadIdx.x][threadIdx.y+5]-s_vx[threadIdx.x][threadIdx.y+2])
			+0.0095703125000f*(s_vx[threadIdx.x][threadIdx.y+6]-s_vx[threadIdx.x][threadIdx.y+1])
			-0.0006975446429f*(s_vx[threadIdx.x][threadIdx.y+7]-s_vx[threadIdx.x][threadIdx.y]);
	convvx[ik]=bx[ik]*convvx[ik]+(bx[ik]-1.0f)*_dx*diff2;
	up[id]+=convvx[ik];
}


//================================================ NJ=10 ========================================================
__global__ void cuda_forward_v_10(float *p, float *vx, float *vz, float _dx, float _dz,	int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+i2*nnz;

	__shared__ float s_p[Block_Size1+9][Block_Size2+9];
    	s_p[threadIdx.x+4][threadIdx.y+4]=p[id]; 
	if (threadIdx.x<4)
	{
		if (blockIdx.x)			s_p[threadIdx.x][threadIdx.y+4]=p[id-4];
		else 				s_p[threadIdx.x][threadIdx.y+4]=0.0f;
	}
	if (threadIdx.x>blockDim.x-6)
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+9][threadIdx.y+4]=p[id+5];
		else				s_p[threadIdx.x+9][threadIdx.y+4]=0.0f;
	}
	if (threadIdx.y<4)
	{
		if (blockIdx.y) 		s_p[threadIdx.x+4][threadIdx.y]=p[id-4*nnz];
		else				s_p[threadIdx.x+4][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>blockDim.y-6)
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x+4][threadIdx.y+9]=p[id+5*nnz];
		else				s_p[threadIdx.x+4][threadIdx.y+9]=0.0f;
	}
	__syncthreads();

	float diff1=1.211242675781250f*(s_p[threadIdx.x+5][threadIdx.y+4]-s_p[threadIdx.x+4][threadIdx.y+4])
			-0.089721679687500f*(s_p[threadIdx.x+6][threadIdx.y+4]-s_p[threadIdx.x+3][threadIdx.y+4])
			+0.013842773437500f*(s_p[threadIdx.x+7][threadIdx.y+4]-s_p[threadIdx.x+2][threadIdx.y+4])
			-0.001765659877232f*(s_p[threadIdx.x+8][threadIdx.y+4]-s_p[threadIdx.x+1][threadIdx.y+4])
			+0.000118679470486f*(s_p[threadIdx.x+9][threadIdx.y+4]-s_p[threadIdx.x][threadIdx.y+4]);
	float diff2= 1.211242675781250f*(s_p[threadIdx.x+4][threadIdx.y+5]-s_p[threadIdx.x+4][threadIdx.y+4])
			-0.089721679687500f*(s_p[threadIdx.x+4][threadIdx.y+6]-s_p[threadIdx.x+4][threadIdx.y+3])
			+0.013842773437500f*(s_p[threadIdx.x+4][threadIdx.y+7]-s_p[threadIdx.x+4][threadIdx.y+2])
			-0.001765659877232f*(s_p[threadIdx.x+4][threadIdx.y+8]-s_p[threadIdx.x+4][threadIdx.y+1])
			+0.000118679470486f*(s_p[threadIdx.x+4][threadIdx.y+9]-s_p[threadIdx.x+4][threadIdx.y]);
	vz[id]=_dz*diff1;
	vx[id]=_dx*diff2;
}

__global__ void cuda_PML_vz_10(float *p, float *convpz, float *bz, float *vz, float _dz, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=blockIdx.x*(nnz-npml)+threadIdx.x+nnz*i2;

	__shared__ float s_p[41][Block_Size2];
    	s_p[threadIdx.x+4][threadIdx.y]=p[id]; 
	if (threadIdx.x<4)
	{
		if (blockIdx.x)			s_p[threadIdx.x][threadIdx.y]=p[id-4];
		else 				s_p[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>26)
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+9][threadIdx.y]=p[id+5];
		else				s_p[threadIdx.x+9][threadIdx.y]=0.0f;
	}
	__syncthreads();

	float diff1=1.1962890625000f*(s_p[threadIdx.x+5][threadIdx.y]-s_p[threadIdx.x+4][threadIdx.y])
			-0.089721679687500f*(s_p[threadIdx.x+6][threadIdx.y]-s_p[threadIdx.x+3][threadIdx.y])
			+0.013842773437500f*(s_p[threadIdx.x+7][threadIdx.y]-s_p[threadIdx.x+2][threadIdx.y])
			-0.001765659877232f*(s_p[threadIdx.x+8][threadIdx.y]-s_p[threadIdx.x+1][threadIdx.y])
			+0.000118679470486f*(s_p[threadIdx.x+9][threadIdx.y]-s_p[threadIdx.x][threadIdx.y]);
	convpz[ik]=bz[ik]*convpz[ik]+(bz[ik]-1.0f)*_dz*diff1;	
	vz[id]+=convpz[ik];
}
__global__ void cuda_PML_vx_10(float *p, float *convpx, float *bx, float *vx, float _dx, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+i2*nnz;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);
	__shared__ float s_p[Block_Size1][41];
    	s_p[threadIdx.x][threadIdx.y+4]=p[id]; 
	if (threadIdx.y<4)
	{
		if (blockIdx.y) 		s_p[threadIdx.x][threadIdx.y]=p[id-4*nnz];
		else				s_p[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>26)
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x][threadIdx.y+9]=p[id+5*nnz];
		else				s_p[threadIdx.x][threadIdx.y+9]=0.0f;
	}
	__syncthreads();

	float diff2=1.1962890625000f*(s_p[threadIdx.x][threadIdx.y+5]-s_p[threadIdx.x][threadIdx.y+4])
			-0.089721679687500f*(s_p[threadIdx.x][threadIdx.y+6]-s_p[threadIdx.x][threadIdx.y+3])
			+0.013842773437500f*(s_p[threadIdx.x][threadIdx.y+7]-s_p[threadIdx.x][threadIdx.y+2])
			-0.001765659877232f*(s_p[threadIdx.x][threadIdx.y+8]-s_p[threadIdx.x][threadIdx.y+1])
			+0.000118679470486f*(s_p[threadIdx.x][threadIdx.y+9]-s_p[threadIdx.x][threadIdx.y]);
	convpx[ik]=bx[ik]*convpx[ik]+(bx[ik]-1.0f)*_dx*diff2;	
	vx[id]+=convpx[ik];
}

__global__ void cuda_forward_up_10(float *up, float *vx, float *vz, float _dx, float _dz, int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+i2*nnz;

	__shared__ float s_vx[Block_Size1][Block_Size2+9];
	__shared__ float s_vz[Block_Size1+9][Block_Size2];
    	s_vx[threadIdx.x][threadIdx.y+5]=vx[id]; 
    	s_vz[threadIdx.x+5][threadIdx.y]=vz[id]; 

	if (threadIdx.x<5)
	{
		if (blockIdx.x)			s_vz[threadIdx.x][threadIdx.y]=vz[id-5];
		else 				s_vz[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>blockDim.x-5)
	{
		if (blockIdx.x<gridDim.x-1)	s_vz[threadIdx.x+9][threadIdx.y]=vz[id+4];
		else				s_vz[threadIdx.x+9][threadIdx.y]=0.0f;
	}
	if (threadIdx.y<5)
	{
		if (blockIdx.y) 		s_vx[threadIdx.x][threadIdx.y]=vx[id-5*nnz];
		else				s_vx[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>blockDim.y-5)
	{
		if (blockIdx.y<gridDim.y-1)	s_vx[threadIdx.x][threadIdx.y+9]=vx[id+4*nnz];
		else				s_vx[threadIdx.x][threadIdx.y+9]=0.0f;
	}
	__syncthreads();	

	float diff1=1.1962890625000f*(s_vz[threadIdx.x+5][threadIdx.y]-s_vz[threadIdx.x+4][threadIdx.y])
			-0.089721679687500f*(s_vz[threadIdx.x+6][threadIdx.y]-s_vz[threadIdx.x+3][threadIdx.y])
			+0.013842773437500f*(s_vz[threadIdx.x+7][threadIdx.y]-s_vz[threadIdx.x+2][threadIdx.y])
			-0.001765659877232f*(s_vz[threadIdx.x+8][threadIdx.y]-s_vz[threadIdx.x+1][threadIdx.y])
			+0.000118679470486f*(s_vz[threadIdx.x+9][threadIdx.y]-s_vz[threadIdx.x][threadIdx.y]);
	float diff2=1.1962890625000f*(s_vx[threadIdx.x][threadIdx.y+5]-s_vx[threadIdx.x][threadIdx.y+4])
			-0.089721679687500f*(s_vx[threadIdx.x][threadIdx.y+6]-s_vx[threadIdx.x][threadIdx.y+3])
			+0.013842773437500f*(s_vx[threadIdx.x][threadIdx.y+7]-s_vx[threadIdx.x][threadIdx.y+2])
			-0.001765659877232f*(s_vx[threadIdx.x][threadIdx.y+8]-s_vx[threadIdx.x][threadIdx.y+1])
			+0.000118679470486f*(s_vx[threadIdx.x][threadIdx.y+9]-s_vx[threadIdx.x][threadIdx.y]);
	up[id]=(_dz*diff1+_dx*diff2);
}
__global__ void cuda_PML_upz_10(float *up, float *convvz, float *bz, float *vz, float _dz, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=blockIdx.x*(nnz-npml)+threadIdx.x+nnz*i2;

	__shared__ float s_vz[41][Block_Size2];
    	s_vz[threadIdx.x+5][threadIdx.y]=vz[id]; 
	if (threadIdx.x<5)
	{
		if (blockIdx.x)			s_vz[threadIdx.x][threadIdx.y]=vz[id-5];
		else 				s_vz[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>27)
	{
		if (blockIdx.x<gridDim.x-1)	s_vz[threadIdx.x+9][threadIdx.y]=vz[id+4];
		else				s_vz[threadIdx.x+9][threadIdx.y]=0.0f;
	}
	__syncthreads();


	float diff1=1.1962890625000f*(s_vz[threadIdx.x+5][threadIdx.y]-s_vz[threadIdx.x+4][threadIdx.y])	 
			-0.089721679687500f*(s_vz[threadIdx.x+6][threadIdx.y]-s_vz[threadIdx.x+3][threadIdx.y])
			+0.013842773437500f*(s_vz[threadIdx.x+7][threadIdx.y]-s_vz[threadIdx.x+2][threadIdx.y])
			-0.001765659877232f*(s_vz[threadIdx.x+8][threadIdx.y]-s_vz[threadIdx.x+1][threadIdx.y])
			+0.000118679470486f*(s_vz[threadIdx.x+9][threadIdx.y]-s_vz[threadIdx.x][threadIdx.y]);
	convvz[ik]=bz[ik]*convvz[ik]+(bz[ik]-1.0f)*_dz*diff1;
	up[id]+=convvz[ik];
}
__global__ void cuda_PML_upx_10(float *up, float *convvx, float *bx, float *vx, float _dx, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+i2*nnz;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);

	__shared__ float s_vx[Block_Size1][41];
    	s_vx[threadIdx.x][threadIdx.y+5]=vx[id]; 
	if (threadIdx.y<5)
	{
		if (blockIdx.y) 		s_vx[threadIdx.x][threadIdx.y]=vx[id-5*nnz];
		else				s_vx[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>27)
	{
		if (blockIdx.y<gridDim.y-1)	s_vx[threadIdx.x][threadIdx.y+9]=vx[id+4*nnz];
		else				s_vx[threadIdx.x][threadIdx.y+9]=0.0f;
	}
	__syncthreads();

	float diff2=1.1962890625000f*(s_vx[threadIdx.x][threadIdx.y+5]-s_vx[threadIdx.x][threadIdx.y+4])
			 -0.089721679687500f*(s_vx[threadIdx.x][threadIdx.y+6]-s_vx[threadIdx.x][threadIdx.y+3])+
			0.013842773437500f*(s_vx[threadIdx.x][threadIdx.y+7]-s_vx[threadIdx.x][threadIdx.y+2])
			 -0.001765659877232f*(s_vx[threadIdx.x][threadIdx.y+8]-s_vx[threadIdx.x][threadIdx.y+1])+
			0.000118679470486f*(s_vx[threadIdx.x][threadIdx.y+9]-s_vx[threadIdx.x][threadIdx.y]);
	convvx[ik]=bx[ik]*convvx[ik]+(bx[ik]-1.0f)*_dx*diff2;
	up[id]+=convvx[ik];
}

// update wavefield p
__global__ void cuda_step_forward(float *vel, float *up, float *p0, float *p1, float dt, int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+i2*nnz;
	float c=dt*vel[id];c=c*c;
	p0[id]=2*p1[id]-p0[id]+c*up[id];
}


//====================================== read and write/save the boundary ==================================
// read and write the inner computation zone boundary coefficients from and into RAM along z direction
// read==flase, write/save boundary; read==true, read the boundary
__global__ void cuda_rw_boundarytb(float *boundarytb, float *p, int npml, int nnz, int nnx, int NJ, bool read)
{
	//int nx=nnx-2*npml;
	int nz=nnz-2*npml;
  	int i1=threadIdx.x+blockDim.x*blockIdx.x;
	int i2=threadIdx.y+blockDim.y*blockIdx.y;
	int id=i1+2*(NJ-1)*i2;
	int i1p=i1+npml-(NJ-1);
	int i2p=i2;
	int idp=i1p+nnz*i2p;

	if (i1<NJ-1 && i2<nnx)
	{
		if(read)
		{
			p[idp]=boundarytb[id];	
			p[idp+nz+NJ-1]=boundarytb[id+NJ-1];
		}
		else	
		{
			boundarytb[id]=p[idp];	
			boundarytb[id+NJ-1]=p[idp+nz+NJ-1];
		}
	}
}

// read and write the inner computation zone boundary coefficients from and into RAM along x direction
// read==flase, write and save boundary; read==true, read the boundary
__global__ void cuda_rw_boundarylr(float *boundarylr, float *p, int npml, int nnz, int nnx, int NJ, bool read)
{
	int nx=nnx-2*npml;
	//int nz=nnz-2*npml;
	int i1=threadIdx.x+blockDim.x*blockIdx.x;
	int i2=threadIdx.y+blockDim.y*blockIdx.y;
	int id=i1+nnz*i2;
	int i1p=i1;
	int i2p=i2+npml-(NJ-1);
	int idp=i1p+nnz*i2p;

	if (i1<nnz && i2<NJ-1)
	{
		if (read)
		{
			p[idp]=boundarylr[id];
			p[idp+nnz*(nx+NJ-1)]=boundarylr[id+nnz*(NJ-1)];
		}
		else
		{
			boundarylr[id]=p[idp];
			boundarylr[id+nnz*(NJ-1)]=p[idp+nnz*(nx+NJ-1)];
		}
	}
}



//========================================== imaging condition ====================================================
__global__ void cuda_cross_correlate(float *Isg, float *Iss, float *sp,	float *gp, int npml, int nnz, int nnx)
{
	int i1=threadIdx.x+blockDim.x*blockIdx.x;
	int i2=threadIdx.y+blockDim.y*blockIdx.y;
	int id=i1+i2*nnz;

    	if(i1>=npml && i1<nnz-npml && i2>=npml && i2<nnx-npml)
    	{
		float ps=sp[id];
		float pg=gp[id];
		Isg[id] += ps*pg;
		Iss[id] += ps*ps;			
    	}
}


__global__ void cuda_imaging(float *Isg, float *Iss, float *I1, float *I2, int npml, int nnz, int nnx)
{
	int nz=nnz-2*npml;
	int nx=nnx-2*npml;
	int i1=threadIdx.x+blockDim.x*blockIdx.x;
	int i2=threadIdx.y+blockDim.y*blockIdx.y;
	int id=i1+i2*nnz;

    	if(i1>=npml && i1<nnz-npml && i2>=npml && i2<nnx-npml) 
	{
		I1[id]+=Isg[id];		// correlation imaging condition
		I2[id]+=Isg[id]/(Iss[id]+EPS);// normalized image
	}
	__syncthreads();

	if 	(i2>=0 	&& i2<npml && i1>=0 && i1<npml) 	// top left
	{	I1[id]=I1[nnz*npml+npml];		I2[id]=I2[nnz*npml+npml];		}
	else if (i2>=npml+nx && i2<nnx && i1>=0 && i1<npml) 	// top right
	{	I1[id]=I1[nnz*(npml+nx-1)+npml]; 	I2[id]=I2[nnz*(npml+nx-1)+npml]; 	}
	else if (i2>=0 && i2<npml && i1>=npml+nz && i1<nnz)	// bottom left
	{	I1[id]=I1[nnz*npml+(npml+nz-1)];	I2[id]=I2[nnz*npml+(npml+nz-1)];	}
	else if (i2>=npml+nx && i2<nnx && i1>=npml+nz && i1<nnz)// bottom right
	{	I1[id]=I1[nnz*(npml+nx-1)+(npml+nz-1)]; I2[id]=I2[nnz*(npml+nx-1)+(npml+nz-1)];	}
	else if (i2>=npml && i2<npml+nx && i1>=0 && i1<npml)	// top
	{	I1[id]=I1[nnz*i2+npml];			I2[id]=I2[nnz*i2+npml];			}
	else if (i2>=npml && i2<npml+nx && i1>=npml+nz && i1<nnz)// bottom
	{	I1[id]=I1[nnz*i2+(npml+nz-1)];		I2[id]=I2[nnz*i2+(npml+nz-1)];		}
	else if (i2>=0 && i2<npml && i1>=npml && i1<npml+nz)	// left
	{	I1[id]=I1[nnz*npml+i1];			I2[id]=I2[nnz*npml+i1];			}
	else if (i2>=npml+nx && i2<nnx && i1>=npml && i1<npml+nz)// right
	{	I1[id]=I1[nnz*(npml+nx-1)+i1];		I2[id]=I2[nnz*(npml+nx-1)+i1];		}
}

__global__ void cuda_laplace_filter(float *Img, float *laplace, float _dz, float _dx, int npml, int nnz, int nnx)
{
	int i1=threadIdx.x+blockDim.x*blockIdx.x;
	int i2=threadIdx.y+blockDim.y*blockIdx.y;
	int id=i1+i2*nnz;
	float diff1=0.0f;
	float diff2=0.0f;
    	if(i1>=npml && i1<nnz-npml && i2>=npml && i2<nnx-npml) 
	{
		diff1=Img[id+1]-2.0*Img[id]+Img[id-1];
		diff2=Img[id+nnz]-2.0*Img[id]+Img[id-nnz];
	}
	laplace[id]=_dz*_dz*diff1+_dx*_dx*diff2;
}
