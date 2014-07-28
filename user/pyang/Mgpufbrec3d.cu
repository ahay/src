/* Backward reconstruction the forward modeled wavefield in 3D with GPU
NB: 2nd order FD, prepared for 3D GPU-based RTM
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
	1) Micikevicius, Paulius. "3D finite difference computation on GPUs
	using CUDA." Proceedings of 2nd Workshop on General Purpose 
	Processing on Graphics Processing Units. ACM, 2009.
	2) Symes, William W., et al. "Computational Strategies For Reverse
	time Migration." 2008 SEG Annual Meeting. Society of Exploration 
	Geophysicists, 2008.
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
#define BlockSize1 16// tile size in 1st-axis
#define BlockSize2 16// tile size in 2nd-axis

static void sf_check_gpu_error (const char *msg) {
    cudaError_t err = cudaGetLastError ();
    if (cudaSuccess != err) { sf_error ("Cuda error: %s: %s", msg, cudaGetErrorString (err)); exit(0);   }
}

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
	int nn1=n1+2;
	int nn2=n2+2;
    	if (id<ns) szxy[id]=(szbeg+id*jsz+1)+nn1*(sxbeg+id*jsx+1)+nn1*nn2*(sybeg+id*jsy+1);
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

__global__ void cuda_step_fd3d(float *p0, float *p1, float *vv, float dtz, float dtx, float dty, int n1, int n2, int n3)
/*< step forward: 3-D FD, order=2 >*/
{
    bool validr = true;
    bool validw = true;
    const int gtid1 = blockIdx.x * blockDim.x + threadIdx.x;
    const int gtid2 = blockIdx.y * blockDim.y + threadIdx.y;
    const int ltid1 = threadIdx.x;
    const int ltid2 = threadIdx.y;
    const int work1 = blockDim.x;
    const int work2 = blockDim.y;

    __shared__ float tile[BlockSize2 + 2][BlockSize1 + 2];
    const int t1 = ltid1 +1;//local thread id in the tile
    const int t2 = ltid2 +1;//local thread id in the tile
    int inIndex = 0;
    int outIndex = 0;

    const int stride2 = n1+2;
    const int stride3 = (n1+2)*(n2+2);

    // Advance inputIndex to target element in inner volume
    inIndex += (gtid2+1) * stride2 + gtid1+1;
    // Check in bounds
    if ((gtid1 >= n1 +1) ||(gtid2 >= n2 +1)) validr = false;
    if ((gtid1 >= n1) || (gtid2 >= n2)) validw = false;

    float infront, behind, current, c1, c2, c3;

    if (validr)	current = p1[inIndex];
    outIndex = inIndex;
    inIndex += stride3;
    if (validr) infront = p1[inIndex];
    inIndex += stride3;

    // Step through the zx-planes
    for (int i3 = 0 ; i3 < n3 ; i3++)
    {
        behind= current;
        current = infront;
        if (validr) infront = p1[inIndex];
        inIndex += stride3;
        outIndex+= stride3;
        __syncthreads();

        if (ltid2 < 1) // Halo above & below
        {
            tile[ltid2][t1]             = p1[outIndex -  stride2];
            tile[ltid2 + work2 + 1][t1] = p1[outIndex + work2 * stride2];
        }
        if (ltid1 < 1) // Halo left & right
        {
            tile[t2][ltid1]            = p1[outIndex -1];
            tile[t2][ltid1 + work1 +1] = p1[outIndex + work1];
        }
        tile[t2][t1] = current;
        __syncthreads();

        // Compute the output value
        c1=c2=c3=-2.0*current; 
	c1 +=tile[t2][t1-1]+ tile[t2][t1+1];
	c2 +=tile[t2-1][t1]+ tile[t2+1][t1];
	c3 +=infront + behind; 
	c1*=dtz*dtz;	
	c2*=dtx*dtx;
	c3*=dty*dty;
	
        if (validw){
	  float vtmp=vv[outIndex]; vtmp=vtmp*vtmp;
	  p0[outIndex]=2.0*p1[outIndex]-p0[outIndex]+vtmp*(c1+c2+c3);
	}
    }
}


void extend3d(float *v0, float*vv, int n1, int n2, int n3)
/*< extend 3d velocity model >*/
{
  int i1, i2, i3, nn1, nn2, nn3;

  nn1=n1+2;
  nn2=n2+2;
  nn3=n3+2;

  for(i3=0; i3<n3; i3++){
    for(i2=0; i2<n2; i2++){
      for(i1=0; i1<n1; i1++){
	vv[(i1+1)+nn1*(i2+1)+nn1*nn2*(i3+1)]=v0[i1+n1*i2+n1*n2*i3];
      }
    }
  }  

    for         (i3=0; i3<nn3; 	i3++) {
	for     (i2=0; i2<nn2; 	i2++) {
		vv[ nn1*i2+nn1*nn2*i3]=vv[1+nn1*i2+nn1*nn2*i3];
		vv[(nn1 -1)+nn1*i2+nn1*nn2]=vv[(nn1-2)+nn1*i2+nn1*nn2];
	}
    }


    for         (i3=0; i3<nn3; 	i3++) {
	    for (i1=0; i1<nn1; 	i1++) {
		vv[i1 +nn1*nn2*i3]=vv[i1+nn1 +nn1*nn2*i3];
		vv[i1+nn1*(nn2 -1)+nn1*nn2*i3]=vv[i1+nn1*(nn2-2)+nn1*nn2*i3];
	    }
    }

	for     (i2=0; i2<nn2; 	i2++) {
	    for (i1=0; i1<nn1; 	i1++) {
		vv[i1+nn1*i2 ]=vv[i1+nn1*i2 ];
		vv[i1+nn1*i2+nn1*nn2*(nn3-1 )]=vv[i1+nn1*i2+nn1*nn2*(nn3-2)];
	    }
	}
}


void window3d(float *a, float *b, int n1, int n2, int n3)
/*< window a 3d subvolume >*/
{
	int i1, i2, i3, nn1, nn2;
	nn1=n1+2;
	nn2=n2+2;
	
	for(i3=0; i3<n3; i3++)
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		a[i1+n1*i2+n1*n2*i3]=b[(i1+1)+nn1*(i2+1)+nn1*nn2*(i3+1)];
	}
}

__global__ void cuda_rw_bndr(float *bndr, float *p1, int n1, int n2, int n3, bool write)
/*< write boundaries out or read them into wavefield variables p>*/
{
	int id=threadIdx.x+blockIdx.x*blockDim.x;
	int nn1=n1+2;
	int nn2=n2+2;
	int nn3=n3+2;
	int iz, ix, iy;

	if	(id<n1)	  		// z-1:	(id+1, 		0, 	0)
	{iz=id+1; ix=0; iy=0;}
	else if (id<2*n1) 		// z-2:	(id-n1+1, 	nn2-1, 	0)
	{iz=id-n1+1; ix=nn2-1; iy=0;}
	else if (id<3*n1) 		// z-3:	(id-2*n1+1,	nn2-1,	nn3-1)
	{iz=id-2*n1+1; ix=nn2-1; iy=nn3-1;}		
	else if (id<4*n1) 		// z-4:	(id-3*n1+1,	0,	nn3-1)
	{iz=id-3*n1+1; ix=0; iy=nn3-1;}
	else if (id<4*n1+n2) 		// x-1: (0,	id-4*n1+1,	0)
	{iz=0; ix=id-4*n1+1; iy=0;}
	else if (id<4*n1+2*n2)		// x-2:	(nn1-1,	id-4*n1-n2+1,	0)
	{iz=nn1-1; ix=id-4*n1-n2+1; iy=0; }
	else if (id<4*n1+3*n2)		// x-3:	(nn1-1,	id-4*n1-2*n2+1,	nn3-1)
	{iz=nn1-1; ix=id-4*n1-2*n2+1; iy=nn3-1;}
	else if (id<4*n1+4*n2)		// x-4:	(0,	id-4*n1-3*n2+1,	nn3-1)
	{iz=0; ix=id-4*n1-3*n2+1; iy=nn3-1;}
	else if (id<4*n1+4*n2+n3)	// y-1:	(0,	0,	id-4*n1-4*n2+1)
	{iz=0; ix=0; iy=id-4*n1-4*n2+1;}
	else if (id<4*n1+4*n2+2*n3)	// y-2:	(nn1-1,	0, 	id-4*n1-4*n2-n3+1)
	{iz=nn1-1; ix=0; iy=id-4*n1-4*n2-n3+1;}
	else if (id<4*n1+4*n2+3*n3) 	// y-3:	(nn1-1,	nn2-1,	id-4*n1-4*n2-2*n3+1)
	{iz=nn1-1; ix=nn2-1; iy=id-4*n1-4*n2-2*n3+1;}
	else if (id<4*n1+4*n2+4*n3) 	// y-4:	(0,	nn2-1,	id-4*n1-4*n2-3*n3+1)
	{iz=0; ix=nn2-1; iy=id-4*n1-4*n2-3*n3+1; }

	if(write){
		bndr[id]=p1[iz+nn1*ix+nn1*nn2*iy];
	}else{
		p1[iz+nn1*ix+nn1*nn2*iy]=bndr[id];
	}
}


int main(int argc, char* argv[])
{
	bool verb;
	int nz, nx, ny, nnz, nnx, nny, ns, nt, kt, it, is, nt_h;
	int szbeg, sxbeg, sybeg, jsz, jsx, jsy;
	int *d_szxy;
	float dz, dx, dy, fm, dt, dtz, dtx, dty, phost;
	float *v0, *vv, *d_wlt, *d_vv, *d_p0, *d_p1, *h_bndr, *d_bndr, *ptr;
	sf_file Fv, Fw;

    	sf_init(argc,argv);
	Fv=sf_input("in");
	Fw=sf_output("out");

    	if (!sf_getbool("verb",&verb)) verb=false; /* verbosit2 */
    	if (!sf_histint(Fv,"n1",&nz)) sf_error("No n1= in input");
    	if (!sf_histint(Fv,"n2",&nx)) sf_error("No n2= in input");
    	if (!sf_histint(Fv,"n3",&ny)) sf_error("No n3= in input");
    	if (!sf_histfloat(Fv,"d1",&dz)) sf_error("No d1= in input");
    	if (!sf_histfloat(Fv,"d2",&dx)) sf_error("No d2= in input");
    	if (!sf_histfloat(Fv,"d3",&dy)) sf_error("No d3= in input");
   	if (!sf_getint("nt",&nt))  sf_error("nt required");
	/* total number of time steps */
    	if (!sf_getfloat("phost",&phost)) phost=0.0;
	/* phost% points on host with zero-copy pinned memory, the rest on device */
    	if (!sf_getint("kt",&kt)) sf_error("kt required");
	/* record wavefield at time kt */
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

	dtz=dt/dz;
	dtx=dt/dx;
	dty=dt/dy;
	nnz=nz+2;
	nnx=nx+2;
	nny=ny+2;
    	v0=(float*)malloc(nz*nx*ny*sizeof(float));
    	vv=(float*)malloc(nnz*nnx*nny*sizeof(float));
	sf_floatread(v0, nz*nx*ny, Fv);// read velocit2 model v0
	extend3d(v0, vv, nz, nx, ny);// init

    	cudaSetDevice(0);// initialize device, default device=0;
	sf_check_gpu_error("Failed to initialize device!");

	dim3 dimg, dimb;// grid and block size
	dimg.x=(nz+BlockSize1-1)/BlockSize1;
	dimg.y=(nx+BlockSize2-1)/BlockSize2;
	dimb.x=BlockSize1;
	dimb.y=BlockSize2;
	nt_h=0.01*phost*nt; 

	/* allocate memory on device */
	cudaMalloc(&d_wlt, nt*sizeof(float));
	cudaMalloc(&d_vv, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&d_p0, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&d_p1, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&d_szxy, ns*sizeof(int));
	cudaHostAlloc(&h_bndr, nt_h*4*(nnz+nnx+nny)*sizeof(float), cudaHostAllocMapped);	
	cudaMalloc(&d_bndr, (nt-nt_h)*4*(nnz+nnx+nny)*sizeof(float));
	sf_check_gpu_error("Failed to allocate required memory!");

	cuda_ricker_wavelet<<<(nt+511)/512, 512>>>(d_wlt, fm, dt, nt);
	cudaMemcpy(d_vv, vv, nnz*nnx*nny*sizeof(float), cudaMemcpyHostToDevice);
	cuda_set_sg<<<1, ns>>>(d_szxy, szbeg, sxbeg, sybeg, jsz, jsx, jsy, ns, nz, nx, ny);

	for(is=0; is<ns; is++){
	  cudaMemset(d_p0, 0, nnz*nnx*nny*sizeof(float));
	  cudaMemset(d_p1, 0, nnz*nnx*nny*sizeof(float));
	  for(it=0; it<nt; it++){
	    cuda_add_source<<<1,1>>>(true, d_p1, &d_wlt[it], &d_szxy[is], 1);
	    cuda_step_fd3d<<<dimg,dimb>>>(d_p0, d_p1, d_vv, dtz, dtx, dty, nz, nx, ny);
	    ptr=d_p0; d_p0=d_p1; d_p1=ptr;

	    if(it<nt_h) cudaHostGetDevicePointer(&ptr, &h_bndr[it*4*(nz+nx+ny)], 0);
	    else  ptr=&d_bndr[(it-nt_h)*4*(nz+nx+ny)];
	    cuda_rw_bndr<<<(4*(nz+nx+ny)+511)/512,512>>>(ptr, d_p0, nz, nx, ny, true);

	    sf_warning("it=%d",it);
	  }

	  ptr=d_p0; d_p0=d_p1; d_p1=ptr;
	  for(it=nt-1; it>-1; it--)
	  {
	    if(it==kt){
	      cudaMemcpy(vv, d_p1, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost);
	      window3d(v0, vv, nz, nx, ny);
	      sf_floatwrite(v0, nz*nx*ny, Fw);
	    }

	    if(it<nt_h) cudaHostGetDevicePointer(&ptr, &h_bndr[it*4*(nz+nx+ny)], 0);
	    else  ptr=&d_bndr[(it-nt_h)*4*(nz+nx+ny)];
	    cuda_rw_bndr<<<(4*(nz+nx+ny)+511)/512,512>>>(ptr, d_p0, nz, nx, ny, false);

	    cuda_step_fd3d<<<dimg,dimb>>>(d_p0, d_p1, d_vv, dtz, dtx, dty, nz, nx, ny);
	    cuda_add_source<<<1,1>>>(false, d_p1, &d_wlt[it], &d_szxy[is], 1);
	    ptr=d_p0; d_p0=d_p1; d_p1=ptr;
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
	free(vv);
	cudaFreeHost(h_bndr);
	cudaFree(d_bndr);

    	exit (0);
}
