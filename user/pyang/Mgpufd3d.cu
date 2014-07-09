/* GPU-based finite difference on 3-D grid
*/
/*
  Copyright (C) 2014  Xi'an Jiaotong University (Pengliang Yang)

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

  Reference:
	Micikevicius, Paulius. "3D finite difference computation on GPUs
	using CUDA." Proceedings of 2nd Workshop on General Purpose 
	Processing on Graphics Processing Units. ACM, 2009.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>

extern "C" {
#include <rsf.h>
}

#ifndef PI
#define PI 	3.141592653589793f
#endif
#define BDIM1 16// tile size in 1st-axis
#define BDIM2 16// tile size in 2nd-axis
#define radius 4// half of the order in space

__constant__ float coeff[5]={-205.0/72.0,8.0/5.0,-1.0/5.0,8.0/315.0,-1.0/560.0};

__global__ void cuda_ricker_wavelet(float *wlt, float fm, float dt, int nt)
/*< generate ricker wavelet with time deley >*/
{
	int it=threadIdx.x+blockDim.x*blockIdx.x;
    	if (it<nt){
	  float tmp = PI*fm*fabsf(it*dt-1.0/fm);//delay the wavelet to exhibit all waveform
	  tmp *=tmp;
	  wlt[it]= (1.0-2.0*tmp)*expf(-tmp);// ricker wavelet at time: t=nt*dt
	}
}



__global__ void cuda_set_sg(int *szxy, int szbeg, int sxbeg, int sybeg, int jsz, int jsx, int jsy, int ns, int n1, int n2, int n3)
/*< set the positions of sources and geophones in whole domain >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<ns) 
	  szxy[id]=(szbeg+id*jsz)+n1*(sxbeg+id*jsx)+n1*n2*(sybeg+id*jsy);
}


__global__ void cuda_add_source(bool add, float *p, float *source, int *szxy, int ns)
/*< add/subtract sources: length of source[]=ns, index stored in szxy[] >*/
{
  int id=threadIdx.x+blockIdx.x*blockDim.x;

  if(id<ns){
    if(add){
      p[szxy[id]]+=source[id];
    }else{
      p[szxy[id]]-=source[id];
    }
  }
}

__global__ void cuda_step_forward3d(float *p0, float *p1, float *vv, float _dz2, float _dx2, float _dy2, int n1, int n2, int n3)
/*< step forward: FD order=8 >*/
{
  __shared__ float s_data[BDIM2+2*radius][BDIM1+2*radius];

  int i1=threadIdx.x+blockIdx.x*blockDim.x;
  int i2=threadIdx.y+blockIdx.y*blockDim.y;
  int in_idx=i1+n1*i2;// index for reading input
  int out_idx=0;// index for writing output
  int stride=n1*n2;

  float infront1, infront2, infront3, infront4;// in front of current slice
  float behind1, behind2, behind3, behind4;// behind the current slice
  float current;

  int t1=threadIdx.x+radius;//1st thread index in shared memory tile
  int t2=threadIdx.y+radius;//2nd thread index in shared memory tile

  // fill the infronts and behinds
  behind3=p1[in_idx]; in_idx+=stride;
  behind2=p1[in_idx]; in_idx+=stride;
  behind1=p1[in_idx]; in_idx+=stride;
  current=p1[in_idx]; out_idx=in_idx; in_idx+=stride; 
  infront1=p1[in_idx]; in_idx+=stride;
  infront2=p1[in_idx]; in_idx+=stride;
  infront3=p1[in_idx]; in_idx+=stride;
  infront4=p1[in_idx]; in_idx+=stride;

  for(int i=radius; i<n3-radius; i++){
    //advance the slice (move the thread-front)
    behind4=behind3;
    behind3=behind2;
    behind2=behind1;
    behind1=current;
    current=infront1;
    infront1=infront2;
    infront2=infront3;
    infront3=infront4;
    infront4=p1[in_idx];

    in_idx+=stride;
    out_idx+=stride;
    __syncthreads();

    // update the slice in shared mem
    if (threadIdx.y<radius){
      if(blockIdx.y) s_data[t2-radius][t1]=p1[out_idx-radius*n1];
      else s_data[t2-radius][t1]=0.0f;
    }
    if (threadIdx.y>=blockDim.y-radius){
      if(blockIdx.y<gridDim.y-1) s_data[t2+radius][t1]=p1[out_idx+radius*n1];
      else s_data[t2+radius][t1]=0.0f;
    }
    if(threadIdx.x<radius){
      if(blockIdx.x) s_data[t2][t1-radius]=p1[out_idx-radius];
      else s_data[t2][t1-radius]=0.0;
    }
    if(threadIdx.x>=blockDim.x-radius){
      if(blockIdx.x<gridDim.x-1) s_data[t2][t1+radius]=p1[out_idx+radius];
      else s_data[t2][t1+radius]=0.0;
    }
    /*
    if(threadIdx.y<radius){
      s_data[threadIdx.y][t1]=(blockIdx.y>0)?p1[out_idx-radius*n1]:0.0;
      s_data[threadIdx.y+BDIM2+radius][t1]=(blockIdx.y<blockDim.y-1)?p1[out_idx+BDIM2*n1]:0.0;
    }
    if(threadIdx.x<radius){
      s_data[t2][threadIdx.x]=(blockIdx.x>0)?p1[out_idx-radius]:0.0;
      s_data[t2][threadIdx.x+BDIM1+radius]=(blockIdx.x<blockDim.x-1)?p1[out_idx+BDIM1]:0.0;
    }
    */
    //update the slice in shared mem
    s_data[t2][t1]=current;
    __syncthreads();

    // compute the output value
    float diff1, diff2, diff3;
    diff1=diff2=diff3=coeff[0]*current;
    diff1+=coeff[1]*(s_data[t2][t1-1]+s_data[t2][t1+1])+
      coeff[2]*(s_data[t2][t1-2]+s_data[t2][t1+2])+
      coeff[3]*(s_data[t2][t1-3]+s_data[t2][t1+3])+
      coeff[4]*(s_data[t2][t1-4]+s_data[t2][t1+4]);
    diff1*=_dz2;
    diff2+=coeff[1]*(s_data[t2-1][t1]+s_data[t2+1][t1])+
      coeff[2]*(s_data[t2-2][t1]+s_data[t2+2][t1])+
      coeff[3]*(s_data[t2-3][t1]+s_data[t2+3][t1])+
      coeff[4]*(s_data[t2-4][t1]+s_data[t2+4][t1]);
    diff2*=_dx2;
    diff3+=coeff[1]*(infront1+behind1)+
      coeff[2]*(infront2+behind2)+
      coeff[3]*(infront3+behind3)+
      coeff[4]*(infront4+behind4);
    diff3*=_dy2;
    p0[out_idx]=2.0f*current-p0[out_idx]+vv[out_idx]*(diff1+diff2+diff3);
    /*
    float div=coeff[0]*current*(_dz2+_dx2+_dy2);
    div+=coeff[1]*(_dy2*(infront1+behind1)+
		   _dx2*(s_data[t2-1][t1]+s_data[t2+1][t1])+
		   _dz2*(s_data[t2][t1-1]+s_data[t2][t1+1]));
    div+=coeff[2]*(_dy2*(infront2+behind2)+
		   _dx2*(s_data[t2-2][t1]+s_data[t2+2][t1])+
		   _dz2*(s_data[t2][t1-2]+s_data[t2][t1+2]));
    div+=coeff[3]*(_dy2*(infront3+behind3)+
		   _dx2*(s_data[t2-3][t1]+s_data[t2+3][t1])+
		   _dz2*(s_data[t2][t1-3]+s_data[t2][t1+3]));
    div+=coeff[4]*(_dy2*(infront4+behind4)+
		   _dx2*(s_data[t2-4][t1]+s_data[t2+4][t1])+
		   _dz2*(s_data[t2][t1-4]+s_data[t2][t1+4]));
    p0[out_idx]=2.0f*current-p0[out_idx]+div*vv[out_idx];
    */
  }
}

void velocity_transform(float *v0, float dt, int n1, int n2, int n3)
 /*< velocity transform: vv=v0*dt; vv<--vv^2 >*/
{
  int i1, i2, i3;
  float tmp;

  for(i3=0; i3<n3; i3++){
    for(i2=0; i2<n2; i2++){
      for(i1=0; i1<n1; i1++){
	tmp=v0[i1+n1*i2+n1*n2*i3]*dt;
	v0[i1+n1*i2+n1*n2*i3]=tmp*tmp;
      }
    }
  }  
}



int main(int argc, char* argv[])
{
	bool verb;
	int nz, nx, ny, ns, nt, kt, it, is, szbeg, sxbeg, sybeg, jsz, jsx, jsy;
	int *d_szxy;
	float dz, dx, dy, fm, dt, _dz2, _dx2, _dy2;
	float *v0, *d_wlt, *d_vv, *d_p0, *d_p1, *ptr;
	sf_file Fv, Fw;

    	sf_init(argc,argv);
	Fv=sf_input("in");
	Fw=sf_output("out");

    	if (!sf_getbool("verb",&verb)) verb=false; /* verbosity */
    	if (!sf_histint(Fv,"n1",&nz)) sf_error("No n1= in input");
    	if (!sf_histint(Fv,"n2",&nx)) sf_error("No n2= in input");
    	if (!sf_histint(Fv,"n3",&ny)) sf_error("No n3= in input");
    	if (!sf_histfloat(Fv,"d1",&dz)) sf_error("No d1= in input");
    	if (!sf_histfloat(Fv,"d2",&dx)) sf_error("No d2= in input");
    	if (!sf_histfloat(Fv,"d3",&dy)) sf_error("No d3= in input");
   	if (!sf_getint("nt",&nt))  sf_error("nt required");
    	if (!sf_getint("kt",&kt)) sf_error("kt required");
	/* record wavefield at time kt */
	if (kt>nt) sf_error("make sure kt<=nt");
   	if (!sf_getfloat("dt",&dt))  sf_error("dt required");
	/* time sampling interval */
   	if (!sf_getfloat("fm",&fm))  fm=20;
	/* dominant frequency of Ricker wavelet */
   	if (!sf_getint("ns",&ns))  ns=1;
	/* number of sources */
	if (!sf_getint("szbeg",&szbeg)) sf_error("No szbeg");
	/* source beginning of z-axis */
	if (!sf_getint("sxbeg",&sxbeg)) sf_error("No sxbeg");
	/* source beginning of x-axis */
	if (!sf_getint("sybeg",&sybeg)) sf_error("No sybeg");
	/* source beginning of y-axis */
	if (!sf_getint("jsz",&jsz)) sf_error("No jsz");
	/* source jump interval in z-axis */
	if (!sf_getint("jsx",&jsx)) sf_error("No jsx");
	/* source jump interval in x-axis */
	if (!sf_getint("jsy",&jsy)) sf_error("No jsy");
	/* source jump interval in y-axis */

	sf_putint(Fw,"n1",nz);
	sf_putint(Fw,"n2",nx);
	sf_putint(Fw,"n3",ny);

	_dz2=1.0/(dz*dz);
	_dx2=1.0/(dx*dx);
	_dy2=1.0/(dy*dy);
    	v0=(float*)malloc(nz*nx*ny*sizeof(float));
	sf_floatread(v0, nz*nx*ny, Fv);// read velocity model
		
    	cudaSetDevice(0);// initialize device, default device=0;
    	cudaError_t err = cudaGetLastError ();
    	if (cudaSuccess != err) 
	sf_warning("Cuda error: Failed to initialize device: %s", cudaGetErrorString(err));

	dim3 dimg((nz+BDIM1-1)/BDIM1, (nx+BDIM2-1)/BDIM2), dimb(BDIM1, BDIM2);

	/* allocate memory on device */
	cudaMalloc(&d_wlt, nt*sizeof(float));
	cudaMalloc(&d_vv, nz*nx*ny*sizeof(float));
	cudaMalloc(&d_p0, nz*nx*ny*sizeof(float));
	cudaMalloc(&d_p1, nz*nx*ny*sizeof(float));
	cudaMalloc(&d_szxy, ns*sizeof(int));

	cuda_ricker_wavelet<<<(nt+511)/512, 512>>>(d_wlt, fm, dt, nt);
	velocity_transform(v0, dt, nz, nx, ny);
	cudaMemcpy(d_vv, v0, nz*nx*ny*sizeof(float), cudaMemcpyHostToDevice);
	cuda_set_sg<<<1, ns>>>(d_szxy, szbeg, sxbeg, sybeg, jsz, jsx, jsy, ns, nz, nx, ny);

	for(is=0; is<ns; is++){
	  cudaMemset(d_p0, 0, nz*nx*ny*sizeof(float));
	  cudaMemset(d_p1, 0, nz*nx*ny*sizeof(float));
	  for(it=0; it<nt; it++){
	    cuda_add_source<<<1,1>>>(true, d_p1, &d_wlt[it], &d_szxy[is], 1);
	    cuda_step_forward3d<<<dimg,dimb>>>(d_p0, d_p1, d_vv, _dz2, _dx2, _dy2, nz, nx, ny);
	    ptr=d_p0; d_p0=d_p1; d_p1=ptr;

	    if(it==kt){
	      cudaMemcpy(v0, d_p0, nz*nx*ny*sizeof(float), cudaMemcpyDeviceToHost);
	      sf_floatwrite(v0, nz*nx*ny, Fw);
	    }
	    sf_warning("it=%d",it);
	  }
	}

	/* free memory on device */
	cudaFree(d_wlt);
	cudaFree(d_vv);
	cudaFree(d_p0);
	cudaFree(d_p1);
	cudaFree(d_szxy);
	free(v0);

    	exit (0);
}
