/* 3D coherence computation using GPU	
*/

/*
  Copyright (C) 2014 Xi'an Jiaotong University, UT Austin (Pengliang Yang)
   
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

  Reference: Marfurt, Kurt J., et al. "Coherency calculations in the presence 
	of structural dip." Geophysics 64.1 (1999): 104-111.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>

extern "C" {
#include <rsf.h>
}

#define BlockSizeX 16
#define BlockSizeY 16

static const int ntw=10;// window radius in t
static const int nxw=1;// window radius in x
static const int nyw=1;// window radius in y


static float *d_u1, *d_u2;

__global__ void coh3(float *u1, float *u2, int dimx, int dimy, int dimz)
/*< C3 calculation
NB: kernel configuration <<<gridDim, blockDim, sizeofsharedmembite>>>  >*/
{
	const int ix=threadIdx.x+blockIdx.x*blockDim.x;
	const int iy=threadIdx.y+blockIdx.y*blockDim.y;
	int id=ix+iy*dimx;

	__shared__ float s_u[BlockSizeY+2*nyw][BlockSizeX+2*nxw];
	const int tx=threadIdx.x+nxw;
	const int ty=threadIdx.y+nyw;

	const int J=(2*nxw+1)*(2*nyw+1);
	const int stride=dimx*dimy;

	float cxy[J][J];
	float u[J];
	float v[J];
	for(int iz=0; iz<dimz; iz++)
	{

		for(int i2=0; i2<J; i2++)
		for(int i1=0; i1<J; i1++)
			cxy[i2][i1]=0.0;
		for(int izw=-ntw; izw<=ntw; izw++)
		if(iz+izw>=0 && iz+izw<dimz)
		{
			int idtmp=id+izw*stride;
			if(threadIdx.y<nyw)// halo above/below
			{
				s_u[threadIdx.y][tx]=(blockIdx.y)?u1[idtmp-nyw*dimx]:0.0;
				s_u[threadIdx.y+BlockSizeY+nyw][tx]=(blockIdx.y<blockDim.y-1)?u1[idtmp+BlockSizeY*dimx]:0.0;
			}
			if(threadIdx.x<nxw)// halo left/right
			{
				s_u[tx][threadIdx.x]=(blockIdx.x)?u1[idtmp-nxw]:0.0;
				s_u[tx][threadIdx.x+BlockSizeX+nxw]=(blockIdx.x<blockDim.x-1)?u1[idtmp+BlockSizeX]:0.0;
			}
		
			for(int iy1=-nyw; iy1<=nyw; iy1++)
			for(int ix1=-nxw; ix1<=nxw; ix1++)
			for(int iy2=-nyw; iy2<=nyw; iy2++)
			for(int ix2=-nxw; ix2<=nxw; ix2++)
			{
				int ixc1=ix+ix1;
				int iyc1=iy+iy1;
				int ixc2=ix+ix2;
				int iyc2=iy+iy2;
				// check in bounds
				if( (ixc1>=0)&&(ixc1<dimx)&&
				    (iyc1>=0)&&(iyc2<dimy)&&
				    (ixc2>=0)&&(ixc2<dimx)&&
				    (iyc2>=0)&&(iyc2<dimy))
				cxy[ix2+nxw+(2*nxw+1)*(iy2+nyw)][ix1+nxw+(2*nxw+1)*(iy1+nyw)]+=s_u[iy1+ty][ix1+tx]*s_u[iy2+ty][ix2+tx];
			}	
		}
		/************************ C3 calculation  ********************/
		float s, m1, m;
		int i,j,k, maxidx;

		s=m1=0;
		for(i=0; i<J; i++) 
		{
			s+=cxy[i][i];// trace{cxy}
			u[i]=1.0;//initialize u
		}

		for(k=0; k<30; k++){//iterations
			for(i=0; i<J; i++){
				v[i]=0.0;
				for(j=0; j<J; j++)
				v[i]+=cxy[i][j]*u[j];
			}

			m=fabsf(v[0]);	maxidx=0;
			for(i=0; i<J; i++){
				 maxidx=(m>fabsf(v[i]))?maxidx:i;
				 m=(m>fabsf(v[i]))?m:fabsf(v[i]);
			}
			m=v[maxidx];
			for(i=0; i<J; i++) u[i]=v[i]/m;

			if(fabsf(m-m1)<1.e-6) break;
			m1=m;
		}

		/************************ End C3 calculation *****************/
		id+=stride;
	}
	

}

int main(int argc, char *argv[])
{
    	sf_file in, out;
    	int n1, n2, n3;
	float ***u1, ***u2;

    	sf_init(argc, argv);
    	in=sf_input("in");	/* 3D seismic data volume */
   	out=sf_output("out");	/* 3D coherence volume */

    	if (!sf_histint(in,"n1",&n1)) 	sf_error("No n1= in input");
    	if (!sf_histint(in,"n2",&n2)) 	sf_error("No n2= in input");
    	if (!sf_histint(in,"n3",&n3)) 	n3=1;	/* default: n3=1 if 2D */

	u1 = sf_floatalloc3(n1, n2, n3);
	u2 = sf_floatalloc3(n1, n2, n3);
	sf_floatread(u1[0][0], n1*n2*n3, in);
	memset(u2[0][0],0, n1*n2*n3*sizeof(float));

    	cudaSetDevice(0);
    	cudaError_t err = cudaGetLastError ();
    	if (cudaSuccess != err) 
	sf_warning("Cuda error: Failed to initialize device: %s", cudaGetErrorString(err));

	cudaMalloc(&d_u1, n1*n2*n3*sizeof(float));
	cudaMalloc(&d_u2, n1*n2*n3*sizeof(float));



    	cudaMemcpy(d_u1, u1[0][0], n1*n2*n3*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemset(d_u2, 0, n1*n2*n3*sizeof(float));




    	cudaMemcpy(u2, d_u2, n1*n2*n3*sizeof(float), cudaMemcpyDeviceToHost);
	sf_floatwrite(u2[0][0], n1*n2*n3, out);

	cudaFree(d_u1);
	cudaFree(d_u2);

	free(**u1); free(*u1); free(u1);
	free(**u2); free(*u2); free(u2);

    	exit(0);
}

