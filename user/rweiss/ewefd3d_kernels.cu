/* GPU Kernel Functions used in sfewefd3d_gpu_p2p and sfewefd3d_gpu_mpi */

/*
  Authors: Robin M. Weiss and Jeffrey Shragge

  This file contains the GPU kernel functions called in the ewefd3d_gpu_mpi and
  ewefd3d_gpu_p2p modules from the Madagascar software package (http://www.reproducilitibly.org).
  These kernel functions are called from both Mewefd3d_gpu_mpi.cu and Mewefd3d_gpu_p2p.cu.
  For more information, see (Weiss and Shragge, "Solving 3D Anisotropic Elastic Wave 
  Equations on Parallel GPU Devices", GEOPHYSICS. http://software.seg.org/2012/0063)
*/

/*
  Copyright (C) 2012 University of Western Australia
  
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


// finite difference stencil coefficients are stored in constant device memory
__device__ __constant__ float C[4] = {4.0f/5.0f, -1.0f/5.0f, 4.0f/105.0f, -1.0f/280.0f};

/*
 * function expand_cpu:
 * Copy array a with dimensions _a into array b and expand data to dimensions _b and
 *         also change the order of the axes in the array
 * NOTE: axes from slowest to fastest in array a are y,x,z.
 * NOTE: axes from slowest to fastest in array b are y,z,x.
 *
 * Input:
 *      a is an input matrix with dimensions x_a y_a z_a
 *
 * Output:
 *      b is the output matrix with dimensions x_b y_b z_b
 *
 * Return:
 * 		void
 */
void expand_cpu(float *a, float *b, int nb, int x_a, int x_b, int y_a, int y_b, int z_a, int z_b){
	
	
	// copy a into center of b, swaping the x and z axes
	for (int iy = 0; iy < y_a; iy++){
		for (int iz = 0; iz < z_a; iz++){
			for (int ix = 0; ix < x_a; ix++){
				b[(iy+nb) * x_b * z_b + (iz+nb) * x_b + (ix+nb)] = a[iy * x_a * z_a + ix * z_a + iz];
			}
		}
	}
	
	// expand z direction
	for (int iy = 0; iy < y_b; iy++){
		for (int iz = 0; iz < nb; iz++){
			for (int ix = 0; ix < x_b; ix++){
				b[iy * x_b * z_b + iz * x_b + ix] = b[iy * x_b * z_b + nb * x_b + ix];
				b[iy * x_b * z_b + (z_b-iz-1) * x_b + ix] = b[iy * x_b * z_b + (z_b-nb-1) * x_b + ix];
			}
		}
	}
	
	// expand x direction
	for (int iy = 0; iy < y_b; iy++){
		for (int iz = 0; iz < z_b; iz++){
			for (int ix = 0; ix < nb; ix++){
				b[iy * x_b * z_b + iz * x_b + ix] = b[iy * x_b * z_b + iz * x_b + nb];
				b[iy * x_b * z_b + iz * x_b + (x_b-ix-1)] = b[iy * x_b * z_b + iz * x_b + (x_b-nb-1)];
			}
		}
	}
	
	// expand y direction
	for (int iy = 0; iy < nb; iy++){
		for (int iz = 0; iz < z_b; iz++){
			for (int ix = 0; ix < x_b; ix++){
				b[(y_b-iy-1) * x_b * z_b + iz * x_b + ix] = b[(y_b-nb-1) * x_b * z_b + iz * x_b + ix];
				b[iy * x_b * z_b + iz * x_b + ix] = b[nb * x_b * z_b + iz * x_b + ix];
			}
		}
	}

}


/*
 * kernel computeRo: 
 *		compute 1/ro * dt^2 for all points in this GPU's portion of the model
 *
 * launch configuration:
 *		one thread per grid point
 *		Grid 	- (ceil(fdm->nxpad/8.0f), ceil(fdm->nzpad/8.0f), ceil(nyinterior/8.0f))
 *		Block 	- (8,8,8)
 *
 * Input:
 *      d_ro		- the density for all points in this GPU's portion of the model
 *		dt			- the time step increment
 * 		nxpad		- the number of points in the x direction of the d_ro array
 *		nzpad		- the number of points in the z direction of the d_ro array
 *		nyinterior	- the number of points in the y direction of the d_ro array
 *
 * Output:
 *      d_ro		- contains 1/ro * dt^2 for all points in this GPU's portion of the model
 *
 * Return:
 * 		void
 */
__global__ void computeRo(float *d_ro, float dt, int nxpad, int nzpad, int nyinterior){
	
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int z = threadIdx.y + blockIdx.y * blockDim.y;
	int y = threadIdx.z + blockIdx.z * blockDim.z;
	
	if (x < nxpad && z < nzpad && y < nyinterior){	// ensure this thread falls inside the model
		
		int addr = (y * nxpad * nzpad) + (z * nxpad) + x;
		
		d_ro[addr] = dt * dt / d_ro[addr];
		
	}
}


/*
 * kernel dispToStrain: 
 *		apply the FD stencil to the displacement arrays to calculate strains in the model
 *
 * launch configuration:
 *		2d plane of threads
 *		Grid 	- (fdm->nxpad-2*NOP)/24.0f, (fdm->nzpad-2*NOP)/24.0f)
 *		Block 	- (24,24,1)
 *
 * Input:
 * 		nxpad			- the number of points in the x direction of the d_t* and d_uo* arrays
 *		nzpad			- the number of points in the z direction of the d_t* and d_uo* arrays
 *		nylocal			- the number of points in the y direction of the d_t* and d_uo* arrays
 *      d_uo*			- displacements for x, y, and z componenets for all points in this GPU's portion of the model
 *		idx, idy, idz	- reciprocal of the spatial sampling in the x, y, and z directions respectively
 *
 * Output:
 *      d_t*			- the strain values in resptective components for all points in this GPU's portion of the model
 *
 * Return:
 * 		void
 */
__global__ void dispToStrain(int nxpad, int nylocal, int nzpad, float *d_uox, float *d_uoy, float *d_uoz, float *d_txx, float *d_tyy, float *d_tzz, float *d_txy, float *d_tyz, float *d_tzx, float idx, float idy, float idz){
	
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int z = threadIdx.y + blockIdx.y * blockDim.y;
	
	int y = 4;	// y = 4 = NOP
	
	__shared__ float sData[32*32*3];
	float *s_uox = &sData[0];
	float *s_uoy = &sData[32*32];
	float *s_uoz = &sData[32*32*2];
	
	float uox_curr;
	float uoy_curr;	// current point in the thread front
	float uoz_curr;
	
	float4 uox_m;
	float4 uoy_m;	// 4 points in the negative y-direction of _curr
	float4 uoz_m;
	
	float4 uox_p;
	float4 uoy_p;	// 4 points in the positive y-direction of _curr
	float4 uoz_p;
	
	int nxNz = nxpad * nzpad; 						// number of elements in one y-slice in the d_ arrays
	int offset = (nxpad * 4) + 4 + (z * nxpad) + x; // location in x-z plane
	
	// load +/- 4 points in the y-direction from global memory into registers
	uox_m.w = d_uox[(nxNz * 0) + offset];
	uox_m.z = d_uox[(nxNz * 1) + offset];
	uox_m.y = d_uox[(nxNz * 2) + offset];
	uox_m.x = d_uox[(nxNz * 3) + offset];
	uox_curr = d_uox[(nxNz * 4) + offset];
	uox_p.x = d_uox[(nxNz * 5) + offset];
	uox_p.y = d_uox[(nxNz * 6) + offset];
	uox_p.z = d_uox[(nxNz * 7) + offset];
	uox_p.w = d_uox[(nxNz * 8) + offset];
	
	uoy_m.w = d_uoy[(nxNz * 0) + offset];
	uoy_m.z = d_uoy[(nxNz * 1) + offset];
	uoy_m.y = d_uoy[(nxNz * 2) + offset];
	uoy_m.x = d_uoy[(nxNz * 3) + offset];
	uoy_curr = d_uoy[(nxNz * 4) + offset];
	uoy_p.x = d_uoy[(nxNz * 5) + offset];
	uoy_p.y = d_uoy[(nxNz * 6) + offset];
	uoy_p.z = d_uoy[(nxNz * 7) + offset];
	uoy_p.w = d_uoy[(nxNz * 8) + offset];
	
	uoz_m.w = d_uoz[(nxNz * 0) + offset];
	uoz_m.z = d_uoz[(nxNz * 1) + offset];
	uoz_m.y = d_uoz[(nxNz * 2) + offset];
	uoz_m.x = d_uoz[(nxNz * 3) + offset];
	uoz_curr = d_uoz[(nxNz * 4) + offset];
	uoz_p.x = d_uoz[(nxNz * 5) + offset];
	uoz_p.y = d_uoz[(nxNz * 6) + offset];
	uoz_p.z = d_uoz[(nxNz * 7) + offset];
	uoz_p.w = d_uoz[(nxNz * 8) + offset];
	
	int localAddr = (32 * 4) + 4 + threadIdx.x + (threadIdx.y * 32);  // address in the s_ arrays for this thread
	
	float txx, tyy, tzz, txy1, txy2, tyz1, tyz2, tzx1, tzx2;
	
	// progress through the 3d volume
	while (y < nylocal - 4){
		
		int globalAddr = (nxNz * y) + offset; // global address into d_ arrays
	
		if (threadIdx.x < 4){ // load left halo into shared memory
			s_uox[localAddr - 4] = d_uox[globalAddr - 4];
			s_uoy[localAddr - 4] = d_uoy[globalAddr - 4];
			s_uoz[localAddr - 4] = d_uoz[globalAddr - 4];
		}
		if(threadIdx.y < 4){ // load top halo into shared memory
			s_uox[localAddr - 128] = d_uox[globalAddr - (4 * nxpad)];
			s_uoy[localAddr - 128] = d_uoy[globalAddr - (4 * nxpad)];
			s_uoz[localAddr - 128] = d_uoz[globalAddr - (4 * nxpad)];
		}
		if (threadIdx.x >= 20){ // load right halo into shared memory
			s_uox[localAddr + 4] = d_uox[globalAddr + 4];
			s_uoy[localAddr + 4] = d_uoy[globalAddr + 4];
			s_uoz[localAddr + 4] = d_uoz[globalAddr + 4];
		}
		if (threadIdx.y >= 20){ // load bottom halo into shared memory
			s_uox[localAddr + 128] = d_uox[globalAddr + (4 * nxpad)];
			s_uoy[localAddr + 128] = d_uoy[globalAddr + (4 * nxpad)];
			s_uoz[localAddr + 128] = d_uoz[globalAddr + (4 * nxpad)];
		}
		
		// load interior points into shared memory
		s_uox[localAddr] = uox_curr;
		s_uoy[localAddr] = uoy_curr;
		s_uoz[localAddr] = uoz_curr;
		
		__syncthreads();

		// apply FD stencils to points from shared memory and registers and store results in output arrays
		txx  = 	C[0]*(s_uox[localAddr + 1] - s_uox[localAddr - 1]) + 
				C[1]*(s_uox[localAddr + 2] - s_uox[localAddr - 2]) +
				C[2]*(s_uox[localAddr + 3] - s_uox[localAddr - 3]) +
				C[3]*(s_uox[localAddr + 4] - s_uox[localAddr - 4]);
		d_txx[globalAddr] = txx * idx;
		
		
		tyy = 	C[0]*(uoy_p.x - uoy_m.x) +
				C[1]*(uoy_p.y - uoy_m.y) +
				C[2]*(uoy_p.z - uoy_m.z) + 
				C[3]*(uoy_p.w - uoy_m.w);
		d_tyy[globalAddr] = tyy * idy;
		
		
		tzz  = 	C[0]*(s_uoz[(localAddr + 32)] - s_uoz[(localAddr - 32)]) +
				C[1]*(s_uoz[(localAddr + 64)] - s_uoz[(localAddr - 64)]) +
				C[2]*(s_uoz[(localAddr + 96)] - s_uoz[(localAddr - 96)]) +
				C[3]*(s_uoz[(localAddr + 128)] - s_uoz[(localAddr - 128)]);
		d_tzz[globalAddr] = tzz * idz;
			

		tzx1 = 	C[0]*(s_uox[(localAddr + 32)] - s_uox[(localAddr - 32)]) +
				C[1]*(s_uox[(localAddr + 64)] - s_uox[(localAddr - 64)]) +
				C[2]*(s_uox[(localAddr + 96)] - s_uox[(localAddr - 96)]) +
				C[3]*(s_uox[(localAddr + 128)] - s_uox[(localAddr - 128)]);
		tzx2 = 	C[0]*(s_uoz[localAddr + 1] - s_uoz[localAddr - 1]) +
				C[1]*(s_uoz[localAddr + 2] - s_uoz[localAddr - 2]) +
				C[2]*(s_uoz[localAddr + 3] - s_uoz[localAddr - 3]) +
				C[3]*(s_uoz[localAddr + 4] - s_uoz[localAddr - 4]);
		d_tzx[globalAddr] = tzx1 * idz + tzx2 * idx;
		

		txy1 =  C[0]*(uox_p.x - uox_m.x) +
				C[1]*(uox_p.y - uox_m.y) +
				C[2]*(uox_p.z - uox_m.z) + 
				C[3]*(uox_p.w - uox_m.w);
		txy2 = 	C[0]*(s_uoy[localAddr + 1] - s_uoy[localAddr - 1]) +
				C[1]*(s_uoy[localAddr + 2] - s_uoy[localAddr - 2]) +
				C[2]*(s_uoy[localAddr + 3] - s_uoy[localAddr - 3]) +
				C[3]*(s_uoy[localAddr + 4] - s_uoy[localAddr - 4]);
		d_txy[globalAddr] = txy1 * idy + txy2 * idx;
		
		
		tyz1 =  C[0]*(uoz_p.x - uoz_m.x) +
				C[1]*(uoz_p.y - uoz_m.y) +
				C[2]*(uoz_p.z - uoz_m.z) + 
				C[3]*(uoz_p.w - uoz_m.w);
		tyz2 = 	C[0]*(s_uoy[(localAddr + 32)] - s_uoy[(localAddr - 32)]) +
				C[1]*(s_uoy[(localAddr + 64)] - s_uoy[(localAddr - 64)]) +
				C[2]*(s_uoy[(localAddr + 96)] - s_uoy[(localAddr - 96)]) +
				C[3]*(s_uoy[(localAddr + 128)] - s_uoy[(localAddr - 128)]);
		d_tyz[globalAddr] = tyz1 * idy + tyz2 * idz;
		
		// advance the thread front by 1 in the y-direction
		y++;
		
		//rotate the data in the +/- data structures to advance the front
		uox_m.w = uox_m.z;  	uoy_m.w = uoy_m.z;		uoz_m.w = uoz_m.z;
		uox_m.z = uox_m.y;  	uoy_m.z = uoy_m.y;		uoz_m.z = uoz_m.y;
		uox_m.y = uox_m.x;  	uoy_m.y = uoy_m.x;		uoz_m.y = uoz_m.x;
		uox_m.x = uox_curr;  	uoy_m.x = uoy_curr;		uoz_m.x = uoz_curr;
		uox_curr = uox_p.x;  	uoy_curr = uoy_p.x;		uoz_curr = uoz_p.x;
		uox_p.x = uox_p.y;  	uoy_p.x = uoy_p.y;		uoz_p.x = uoz_p.y;
		uox_p.y = uox_p.z;  	uoy_p.y = uoy_p.z;		uoz_p.y = uoz_p.z;
		uox_p.z = uox_p.w;  	uoy_p.z = uoy_p.w;		uoz_p.z = uoz_p.w;
		
		// get the value y+4 away from the current position of the thread front
		uox_p.w = d_uox[(nxNz * (y+4)) + offset];
		uoy_p.w = d_uoy[(nxNz * (y+4)) + offset];
		uoz_p.w = d_uoz[(nxNz * (y+4)) + offset];
		
		__syncthreads();
	
	}
	
}


/*
 * kernel strainToStress: 
 *		multiply strain values from d_t* arrays by coefficients from d_c* arrays to calculate stress values
 *
 * launch configuration:
 *		one thread per grid point
 *		Grid 	- (ceil(fdm->nxpad/192.0f), fdm->nzpad, nyinterior)
 *		Block 	- (192,1,1)
 *
 * Input:
 *		gpuID			- the GPU this kernel is executing on
 * 		nxpad			- number of points in the x direction of the d_t* and d_c* arrays
 *		nzpad			- number of points in the z direction of the d_t* and d_c* arrays
 *		nyinterior		- number of points in the y direction of the d_c* arrays
 *      d_c*			- stiffness coefficients for all points in this GPU's portion of the model
 *		d_t*			- strain values for all points in this GPU's portion of the model
 *
 * Output:
 *      d_t*			- the stress values in resptective components for all points in this GPU's portion of the model
 *
 * Return:
 * 		void
 */
__global__ void strainToStress(int gpuID, int nxpad, int nzpad, int nyinterior, float *d_c11, float *d_c12, float *d_c13, float *d_c22, float *d_c23, float *d_c33, float *d_c44, float *d_c55, float *d_c66, float *d_txx, float *d_tyy, float *d_tzz, float *d_txy, float *d_tyz, float *d_tzx){
	
	
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int z = blockIdx.y;
	int y = blockIdx.z;
	
	if (x < nxpad){
		int cAddr = y * nzpad * nxpad + z * nxpad + x;			// stiffness arrays do not have additional ghost cells
		int tAddr;
		if (gpuID == 0){
			tAddr = y * nzpad * nxpad + z * nxpad + x;			// GPU 0 does not have any ghost cells in strain/stress arrays
		}
		else {
			tAddr = (y+4) * nzpad * nxpad + z * nxpad + x;		// all other GPUs have ghost cells in strain/stress arrays
		}
		
		// these values are all re-used, store them in registers
		float c12 = d_c12[cAddr];
		float c13 = d_c13[cAddr];
		float c22 = d_c22[cAddr];
		float c23 = d_c23[cAddr];
		
		float txx = d_txx[tAddr];
		float tyy = d_tyy[tAddr];
		float tzz = d_tzz[tAddr];
		
		d_txx[tAddr] = d_c11[cAddr] * txx
					   + c12 * tyy
					   + c13 * tzz;
	
	    d_tyy[tAddr] = c12 * txx
					 + c22 * tyy
				  	 + c23 * tzz;
	
	    d_tzz[tAddr] = c13 * txx
				     + c23 * tyy
				   + d_c33[cAddr] * tzz;
    
		// store stresses in output arrays
	    d_tyz[tAddr] = d_c44[cAddr] * d_tyz[tAddr];
	    d_tzx[tAddr] = d_c55[cAddr] * d_tzx[tAddr];
		d_txy[tAddr] = d_c66[cAddr] * d_txy[tAddr];
    
	}
	
}


/*
 * kernel freeSurf: 
 *		set all stress tensors containing a z-component to 0 in the z=0 to z=nb region of the domain
 *
 * launch configuration:
 *		one thread per grid point
 *		Grid 	- 	(ceil(fdm->nxpad/8.0f), ceil(fdm->nb/8.0f), ceil(nyinterior/8.0f))
 *		Block 	-	(8,8,8)
 *
 * Input:
 *		gpuID			- the GPU this kernel is executing on
 * 		nxpad			- the number of points in the x direction of the d_t* arrays
 *		nzpad			- the number of points in the z direction of the d_t* arrays
 *		nyinterior		- the number of points in the y direction of the d_t* arrays
 *		nb				- amount of padding added to the model in the low side of the XY-plane
 *		d_t*			- stress values involving the z-component for all points in this GPU's portion of the model
 *
 * Output:
 *      d_t*			- stress values involving the z-component for all points in this GPU's portion of the model
 *
 * Return:
 * 		void
 */
__global__ void freeSurf(int gpuID, int nxpad, int nyinterior, int nzpad, int nb, float *d_tzz, float *d_tyz, float *d_tzx){
	
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int z = threadIdx.y + blockIdx.y * blockDim.y;
	int y = threadIdx.z + blockIdx.z * blockDim.z;
	
	
	if (x < nxpad && y < nyinterior && z < nb){		// ensure this thread falls in the free region of hte model
		int tAddr;	
		if (gpuID == 0){
			tAddr = y * nzpad * nxpad + z * nxpad + x;			// GPU 0 does not have any additional halo cells in strain/stress arrays
		}
		else {
			tAddr = (y+4) * nzpad * nxpad + z * nxpad + x;		// all other GPUs have additional halo cells in strain/stress arrays
		}
		
		d_tzz[tAddr] = 0.0f;
		d_tyz[tAddr] = 0.0f;
		d_tzx[tAddr] = 0.0f;
		
	}
	
}


/*
 * kernel lint3d_bell_gpu: 
 *		inject stress or accel source from the d_ww array into the array d_uu at a x,z,y coordinate given by d_jx[], d_jz[], and d_jy[]
 *			using interpolation coefficients in the d_Sw arrays
 *
 * launch configuration:
 *		2d thread blocks with side length 2*nbell+1 and 1 block per source
 *		Grid 	- 	(ns, 1, 1);
 *		Block 	-	(2 * nbell + 1, 2 * nbell + 1, 1);
 *
 * Input:
 *		gpuID			- the GPU this kernel is executing on
 * 		nxpad			- number of points in the x direction of the d_t* and d_c* arrays
 *		nzpad			- number of points in the z direction of the d_t* and d_c* arrays
 *		nyinterior		- number of points in the y direction of the d_c* arrays
 *		it				- current iteration number
 *		nc				- number of vector components
 *		ns				- number of sources
 *		c				- component of source being injected
 *		d_ww			- values to inject for each source
 *		d_j*			- x-, y-, and z-coordinates of each source location
 *		d_Sw*			- linear interpolation coefficients
 *
 * Output:
 *      d_uu			- the input d_uu array with the source injected
 *
 * Return:
 * 		void
 */
__global__ void lint3d_bell_gpu(int gpuID, int it, int nc, int ns, int c, int nbell, int nxpad, int nyinterior, int nzpad, float *d_uu, float *d_bell, int *d_jx, int *d_jz, int *d_jy, float *d_ww, float *d_Sw000, float *d_Sw001, float *d_Sw010, float *d_Sw011, float *d_Sw100, float *d_Sw101, float *d_Sw110, float *d_Sw111){
	
	int ix = threadIdx.x; // 0 to 2*nbell + 1
	int iy = threadIdx.y; // 0 to 2*nbell + 1	
	int ia = blockIdx.x;
	
	int haloCorrection = 0;	// GPU 0 does not have any additional halo cells in stress/acceleration arrays
	
	if (gpuID != 0){
		haloCorrection = 4;	// all other GPUs have additional halo cells in stress/acceleration arrays
	}
	
	// if this thread will be injecting within its portion of the domain...
	if ((d_jy[ia] - nbell + iy) >= (gpuID * nyinterior) && ((d_jy[ia] - nbell) + iy) < (gpuID * nyinterior + nyinterior)){
		for (int iz = 0; iz < 2*nbell + 1; iz++){
			float wa = d_ww[it * nc * ns + c * ns + ia] * d_bell[(iy * (2*nbell+1) * (2*nbell+1)) + (iz * (2*nbell+1)) + ix];
			
			float wa000 = -(wa * d_Sw000[ia]);
			float wa001 = -(wa * d_Sw001[ia]);
			float wa010 = -(wa * d_Sw010[ia]);
			float wa011 = -(wa * d_Sw011[ia]);
	                                                                                                                                                                                               
			atomicAdd(&d_uu[((d_jy[ia] - (gpuID * nyinterior) - nbell) + iy + haloCorrection) * nxpad * nzpad + ((d_jx[ia] - nbell) + ix    ) + ((d_jz[ia] - nbell) + iz    ) * nxpad], wa000);
			atomicAdd(&d_uu[((d_jy[ia] - (gpuID * nyinterior) - nbell) + iy + haloCorrection) * nxpad * nzpad + ((d_jx[ia] - nbell) + ix    ) + ((d_jz[ia] - nbell) + iz + 1) * nxpad], wa001);
			atomicAdd(&d_uu[((d_jy[ia] - (gpuID * nyinterior) - nbell) + iy + haloCorrection) * nxpad * nzpad + ((d_jx[ia] - nbell) + ix + 1) + ((d_jz[ia] - nbell) + iz    ) * nxpad], wa010);
			atomicAdd(&d_uu[((d_jy[ia] - (gpuID * nyinterior) - nbell) + iy + haloCorrection) * nxpad * nzpad + ((d_jx[ia] - nbell) + ix + 1) + ((d_jz[ia] - nbell) + iz + 1) * nxpad], wa011);
	
		}
	
	}
		
	// if this thread will be injecting within its portion of the domain...
	if ((d_jy[ia] - nbell + iy + 1) >= (gpuID * nyinterior + 1) && (d_jy[ia] - nbell + iy + 1) < (gpuID * nyinterior + nyinterior)){
		for (int iz = 0; iz < 2*nbell + 1; iz++){
			float wa = d_ww[it * nc * ns + c * ns + ia] * d_bell[(iy * (2*nbell+1) * (2*nbell+1)) + (iz * (2*nbell+1)) + ix];
			
			float wa100 = -(wa * d_Sw100[ia]);
			float wa101 = -(wa * d_Sw101[ia]);
			float wa110 = -(wa * d_Sw110[ia]);
			float wa111 = -(wa * d_Sw111[ia]);
			
			atomicAdd(&d_uu[((d_jy[ia] - (gpuID * nyinterior) - nbell) + iy + 1 + haloCorrection) * nxpad * nzpad + ((d_jx[ia] - nbell) + ix    ) + ((d_jz[ia] - nbell) + iz    ) * nxpad], wa100);
			atomicAdd(&d_uu[((d_jy[ia] - (gpuID * nyinterior) - nbell) + iy + 1 + haloCorrection) * nxpad * nzpad + ((d_jx[ia] - nbell) + ix    ) + ((d_jz[ia] - nbell) + iz + 1) * nxpad], wa101);
			atomicAdd(&d_uu[((d_jy[ia] - (gpuID * nyinterior) - nbell) + iy + 1 + haloCorrection) * nxpad * nzpad + ((d_jx[ia] - nbell) + ix + 1) + ((d_jz[ia] - nbell) + iz    ) * nxpad], wa110);
			atomicAdd(&d_uu[((d_jy[ia] - (gpuID * nyinterior) - nbell) + iy + 1 + haloCorrection) * nxpad * nzpad + ((d_jx[ia] - nbell) + ix + 1) + ((d_jz[ia] - nbell) + iz + 1) * nxpad], wa111);
			
		}
	}
				
	// if the source bell crosses over a sub-domain boundary in the y direction, 
	// the threads on the high-side of the boundary have to do the work of the threads on the low-side
	if (gpuID != 0 && (d_jy[ia] - nbell < (gpuID * nyinterior)) && (d_jy[ia] + nbell >= (gpuID * nyinterior))){		
		if ((d_jy[ia] - nbell + iy) == (gpuID * nyinterior)){
			iy--;
			for (int iz = 0; iz < 2*nbell + 1; iz++){
				float wa = d_ww[it * nc * ns + c * ns + ia] * d_bell[(iy * (2*nbell+1) * (2*nbell+1)) + (iz * (2*nbell+1)) + ix];
				
				float wa100 = -1.0f * wa * d_Sw100[ia];
				float wa110 = -1.0f * wa * d_Sw110[ia];
				float wa101 = -1.0f * wa * d_Sw101[ia];
				float wa111 = -1.0f * wa * d_Sw111[ia];
	                                                                                                                        
				atomicAdd(&d_uu[4 * nxpad * nzpad + ((d_jx[ia] - nbell) + ix    ) + ((d_jz[ia] - nbell) + iz    ) * nxpad], wa100);
				atomicAdd(&d_uu[4 * nxpad * nzpad + ((d_jx[ia] - nbell) + ix    ) + ((d_jz[ia] - nbell) + iz + 1) * nxpad], wa101);
				atomicAdd(&d_uu[4 * nxpad * nzpad + ((d_jx[ia] - nbell) + ix + 1) + ((d_jz[ia] - nbell) + iz    ) * nxpad], wa110);
				atomicAdd(&d_uu[4 * nxpad * nzpad + ((d_jx[ia] - nbell) + ix + 1) + ((d_jz[ia] - nbell) + iz + 1) * nxpad], wa111);
	
			}
		}
	}

}



/*
 * kernel stressToAccel: 
 *		apply FD stencil to stress arrays (d_t*) to obtain acceleration values for grid points
 *
 * launch configuration:
 *		2d plane of threads
 *		Grid 	- (fdm->nxpad-2*NOP)/24.0f, (fdm->nzpad-2*NOP)/24.0f)
 *		Block 	- (24,24,1)
 *
 * Input:
 * 		nxpad			- number of points in the x direction of the d_t* and d_ua* arrays
 *		nzpad			- number of points in the z direction of the d_t* and d_ua* arrays
 *		nylocal			- number of points in the y direction of the d_t* and d_ua* arrays
 *		idx, idy, idz	- reciprocal of the spatial sampling in the x, y, and z directions respectively
 *		d_t*			- stresses for all grid points in this sub-domain for all vector components
 *
 * Output:
 *      d_ua*			- accelerations for all grid points in this sub-domain
 *
 * Return:
 * 		void
 */
__global__ void stressToAccel(int nxpad, int nzpad, int nylocal, float idx, float idy, float idz, float *d_txx, float *d_tyy, float *d_tzz, float *d_txy, float *d_tzx, float *d_tyz, float *d_uax, float *d_uay, float *d_uaz){
	
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int z = threadIdx.y + blockIdx.y * blockDim.y;
	
	int y = 4;	// y = NOP = 4
	
	__shared__ float sData[32*32*5];
	float *s_txx = &sData[0];
	float *s_txy = &sData[32*32];
	float *s_tyz = &sData[32*32*2];
	float *s_tzx = &sData[32*32*3];
	float *s_tzz = &sData[32*32*4];
	
	float txy_curr;
	float tyy_curr;		// current point in the thread front
	float tyz_curr;
	
	float4 txy_m;
	float4 tyy_m;	// 4 points in the negative y-direction of _curr
	float4 tyz_m;
	
	float4 txy_p;
	float4 tyy_p;	// 4 points in the positive y-direction of _curr
	float4 tyz_p;
	
	int nxNz = nxpad * nzpad; 							// number of elements in one y-slice in the d_ arrays
	int offset = (nxpad * 4) + 4 + (z * nxpad) + x; // location in x-z plane
	
	txy_m.w = d_txy[(nxNz * 0) + offset];
	txy_m.z = d_txy[(nxNz * 1) + offset];
	txy_m.y = d_txy[(nxNz * 2) + offset];
	txy_m.x = d_txy[(nxNz * 3) + offset];
	txy_curr = d_txy[(nxNz * 4) + offset];
	txy_p.x = d_txy[(nxNz * 5) + offset];
	txy_p.y = d_txy[(nxNz * 6) + offset];
	txy_p.z = d_txy[(nxNz * 7) + offset];
	txy_p.w = d_txy[(nxNz * 8) + offset];
	
	tyy_m.w = d_tyy[(nxNz * 0) + offset];
	tyy_m.z = d_tyy[(nxNz * 1) + offset];
	tyy_m.y = d_tyy[(nxNz * 2) + offset];
	tyy_m.x = d_tyy[(nxNz * 3) + offset];
	tyy_curr = d_tyy[(nxNz * 4) + offset];
	tyy_p.x = d_tyy[(nxNz * 5) + offset];
	tyy_p.y = d_tyy[(nxNz * 6) + offset];
	tyy_p.z = d_tyy[(nxNz * 7) + offset];
	tyy_p.w = d_tyy[(nxNz * 8) + offset];
	
	tyz_m.w = d_tyz[(nxNz * 0) + offset];
	tyz_m.z = d_tyz[(nxNz * 1) + offset];
	tyz_m.y = d_tyz[(nxNz * 2) + offset];
	tyz_m.x = d_tyz[(nxNz * 3) + offset];
	tyz_curr = d_tyz[(nxNz * 4) + offset];
	tyz_p.x = d_tyz[(nxNz * 5) + offset];
	tyz_p.y = d_tyz[(nxNz * 6) + offset];
	tyz_p.z = d_tyz[(nxNz * 7) + offset];
	tyz_p.w = d_tyz[(nxNz * 8) + offset];
	
	int localAddr = (32 * 4) + 4 + threadIdx.x + (threadIdx.y * 32);  // address in the s_ arrays for this thread
	
	float uax1;	// Dx(txx)
	float uax2; // Dy(txy)
	float uax3; // Dz(tzx)
	float uay1; // Dx(txy)
	float uay2; // Dy(tyy)
	float uay3; // Dz(tyz)
	float uaz1; // Dx(tzx)
	float uaz2; // Dy(tyz)
	float uaz3; // Dz(tzz)
	
	// progress through the 3d volume
	while (y < nylocal - 4){
		
		int globalAddr = (nxNz * y) + offset; // global address in d_ arrays
	
		if (threadIdx.x < 4){ 	// load left halo into shared memory
			s_txx[localAddr - 4] = d_txx[globalAddr - 4];
			s_txy[localAddr - 4] = d_txy[globalAddr - 4];
			s_tzx[localAddr - 4] = d_tzx[globalAddr - 4];
		}
		if (threadIdx.y < 4){ 	// load top halo into shared memory
			s_tzx[localAddr - 128] = d_tzx[globalAddr - (4 * nxpad)];
			s_tyz[localAddr - 128] = d_tyz[globalAddr - (4 * nxpad)];
			s_tzz[localAddr - 128] = d_tzz[globalAddr - (4 * nxpad)];
		}
		if (threadIdx.x >= 20){ // load right halo into shared memory
			s_txx[localAddr + 4] = d_txx[globalAddr + 4];
			s_txy[localAddr + 4] = d_txy[globalAddr + 4];
			s_tzx[localAddr + 4] = d_tzx[globalAddr + 4];
		}
		if (threadIdx.y >= 20){ // load bottom halo into shared memory
			s_tzx[localAddr + 128] = d_tzx[globalAddr + (4 * nxpad)];
			s_tyz[localAddr + 128] = d_tyz[globalAddr + (4 * nxpad)];
			s_tzz[localAddr + 128] = d_tzz[globalAddr + (4 * nxpad)];
		}
		
		// load interior data into shared memory
		s_txx[localAddr] = d_txx[globalAddr];
		s_tzx[localAddr] = d_tzx[globalAddr];
		s_txy[localAddr] = txy_curr;
		s_tyz[localAddr] = tyz_curr;
		s_tzz[localAddr] = d_tzz[globalAddr];
		
		__syncthreads();

		// apply FD stencils to points from shared memory
		uax1 = 	C[0]*(s_txx[localAddr + 1] - s_txx[localAddr - 1]) + 
				C[1]*(s_txx[localAddr + 2] - s_txx[localAddr - 2]) +
				C[2]*(s_txx[localAddr + 3] - s_txx[localAddr - 3]) +
				C[3]*(s_txx[localAddr + 4] - s_txx[localAddr - 4]);

		uax2 = 	C[0]*(txy_p.x - txy_m.x) +
				C[1]*(txy_p.y - txy_m.y) +
				C[2]*(txy_p.z - txy_m.z) + 
				C[3]*(txy_p.w - txy_m.w);

		uax3 = 	C[0]*(s_tzx[(localAddr + 32)] - s_tzx[(localAddr - 32)]) +
				C[1]*(s_tzx[(localAddr + 64)] - s_tzx[(localAddr - 64)]) +
				C[2]*(s_tzx[(localAddr + 96)] - s_tzx[(localAddr - 96)]) +
				C[3]*(s_tzx[(localAddr + 128)] - s_tzx[(localAddr - 128)]);

		d_uax[globalAddr] = uax1 * idx + uax2 * idy + uax3 * idz;
		

		uay1 = 	C[0]*(s_txy[localAddr + 1] - s_txy[localAddr - 1]) + 
				C[1]*(s_txy[localAddr + 2] - s_txy[localAddr - 2]) +
				C[2]*(s_txy[localAddr + 3] - s_txy[localAddr - 3]) +
				C[3]*(s_txy[localAddr + 4] - s_txy[localAddr - 4]);

		uay2 =  C[0]*(tyy_p.x - tyy_m.x) +
				C[1]*(tyy_p.y - tyy_m.y) +
				C[2]*(tyy_p.z - tyy_m.z) + 
				C[3]*(tyy_p.w - tyy_m.w);

		uay3 =  C[0]*(s_tyz[(localAddr + 32)] - s_tyz[(localAddr - 32)]) +
				C[1]*(s_tyz[(localAddr + 64)] - s_tyz[(localAddr - 64)]) +
				C[2]*(s_tyz[(localAddr + 96)] - s_tyz[(localAddr - 96)]) +
				C[3]*(s_tyz[(localAddr + 128)] - s_tyz[(localAddr - 128)]);
		
		d_uay[globalAddr] = uay1 * idx + uay2 * idy + uay3 * idz;
		
		
		uaz1 =	C[0]*(s_tzx[localAddr + 1] - s_tzx[localAddr - 1]) + 
				C[1]*(s_tzx[localAddr + 2] - s_tzx[localAddr - 2]) +
				C[2]*(s_tzx[localAddr + 3] - s_tzx[localAddr - 3]) +
				C[3]*(s_tzx[localAddr + 4] - s_tzx[localAddr - 4]);
		
		uaz2 =  C[0]*(tyz_p.x - tyz_m.x) +
				C[1]*(tyz_p.y - tyz_m.y) +
				C[2]*(tyz_p.z - tyz_m.z) + 
				C[3]*(tyz_p.w - tyz_m.w);
		
		uaz3 = 	C[0]*(s_tzz[(localAddr + 32)] - s_tzz[(localAddr - 32)]) +
				C[1]*(s_tzz[(localAddr + 64)] - s_tzz[(localAddr - 64)]) +
				C[2]*(s_tzz[(localAddr + 96)] - s_tzz[(localAddr - 96)]) +
				C[3]*(s_tzz[(localAddr + 128)] - s_tzz[(localAddr - 128)]);
		
		d_uaz[globalAddr] = uaz1 * idx + uaz2 * idy + uaz3 * idz;

		// advance the thread front by 1 in the y-direction
		y++;
		
		//rotate the data in the +/- data structures to advance the front
		txy_m.w = txy_m.z;		tyy_m.w = tyy_m.z;		tyz_m.w = tyz_m.z;
		txy_m.z = txy_m.y;		tyy_m.z = tyy_m.y;		tyz_m.z = tyz_m.y;
		txy_m.y = txy_m.x;		tyy_m.y = tyy_m.x;		tyz_m.y = tyz_m.x;
		txy_m.x = txy_curr;		tyy_m.x = tyy_curr;		tyz_m.x = tyz_curr;
		txy_curr = txy_p.x;		tyy_curr = tyy_p.x;		tyz_curr = tyz_p.x;
		txy_p.x = txy_p.y;		tyy_p.x = tyy_p.y;		tyz_p.x = tyz_p.y;
		txy_p.y = txy_p.z;		tyy_p.y = tyy_p.z;		tyz_p.y = tyz_p.z;
		txy_p.z = txy_p.w;		tyy_p.z = tyy_p.w;		tyz_p.z = tyz_p.w;
		
		// get the value y+4 away from the current position of the thread front
		txy_p.w = d_txy[(nxNz * (y+4)) + offset];
		tyy_p.w = d_tyy[(nxNz * (y+4)) + offset];
		tyz_p.w = d_tyz[(nxNz * (y+4)) + offset];
		
		__syncthreads();
	
	}
	
}



/*
 * kernel stepTime: 
 *		calculate displacements of grid points in the next time step
 *
 * launch configuration:
 *		one thread per grid point
 *		Grid 	- 	(ceil(fdm->nxpad/192.0f), fdm->nzpad, nyinterior);
 *		Block	-	(192,1,1);
 *
 * Input:
 *		gpuID			- the gpu this kernel is executing on
 * 		nxpad			- number of points in the x direction of the d_u* arrays
 *		nzpad			- number of points in the z direction of the d_u* arrays
 *		nyinterior		- number of points in the y direction of the d_u* arrays
 *		d_uo*			- displacements in the current time step for all grid points
 *		d_um*			- displacements in the previous stime step for all grid points
 *		d_ua*			- accelerations for all grid points
 *		d_ro			- density field
 *
 * Output:
 *      d_up*			- displacements in the next time step for all grid points
 *
 * Return:
 * 		void
 */
__global__ void stepTime(int gpuID, int nxpad, int nyinterior, int nzpad, float *d_ro, float *d_uox, float *d_umx, float *d_uax, float *d_upx, float *d_uoy, float *d_umy, float *d_uay, float *d_upy, float *d_uoz, float *d_umz, float *d_uaz, float *d_upz){
	
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int z = blockIdx.y;
	int y = blockIdx.z;
	
	if (x < nxpad){
		
		int addr;
		
		if (gpuID == 0){		// GPU 0 does not have any additional halo cells in arrays
			addr = y * nxpad * nzpad + z * nxpad + x;	
		}
		else {					// all other GPUs have additional halo cells in arrays
			addr = (y+4) * nxpad * nzpad + z * nxpad + x;
		}
		
		float ro = d_ro[y * nxpad * nzpad + z * nxpad + x];	// d_ro does not have any ghost cells
		
		// calculate forward time step for this grid point
		d_upx[addr] = 2.0f * d_uox[addr] - d_umx[addr] + d_uax[addr] * ro;
		d_upz[addr] = 2.0f * d_uoz[addr] - d_umz[addr] + d_uaz[addr] * ro;
		d_upy[addr] = 2.0f * d_uoy[addr] - d_umy[addr] + d_uay[addr] * ro;
		
	}
	
}



/*
 * kernel extract_gpu: 
 *		extract receiver data from the GPU without performing linear interpolation
 *
 * launch configuration:
 *		one thread per receiver
 *		Grid 	- 	(MIN(nr,ceil(nr/1024.0f)), 1, 1);
 *		Block 	-	(MIN(nr, 1024), 1, 1);
 *
 * Input:
 *		gpuID			- the gpu this kernel is executing on
 * 		nxpad			- number of points in the x direction of the d_u* arrays
 *		nzpad			- number of points in the z direction of the d_u* arrays
 *		nyinterior		- number of points in the y direction of the d_u* arrays
 *		nr				- number of receivers
 *		d_uo*			- displacements in the current time step for all grid points
 *		d_Rj*			- receiver coordinates
 *
 * Output:
 *      d_dd			- displacements at each receiver location (receiver data)
 * 							NOTE: madagascar expects components in data output to be 1=z, 2=x, 3=y
 *
 * Return:
 * 		void
 */
__global__ void extract_gpu(int gpuID, float *d_dd, int nr, int nxpad, int nyinterior, int nzpad, float *d_uoz, float *d_uox, float *d_uoy, int *d_Rjz, int *d_Rjx, int *d_Rjy){

	int rr = threadIdx.x + blockIdx.x * blockDim.x;		// rr = the receiver this thread is extracting data for
	
	int haloCorrection = 0;	// GPU 0 does not have any additional halo cells in arrays
	
	if (gpuID != 0){
		haloCorrection = 4;	// all other GPUs have additional halo cells in arrays
	}
	
	if (rr < nr && 
		d_Rjy[rr] >= (gpuID * nyinterior) && 
		d_Rjy[rr] < (gpuID * nyinterior + nyinterior)){
			
			d_dd[rr] = 			d_uoz[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection) * nxpad * nzpad + d_Rjx[rr] + d_Rjz[rr] * nxpad ];
							
			d_dd[rr+nr] = 		d_uox[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection) * nxpad * nzpad + d_Rjx[rr] + d_Rjz[rr] * nxpad ];
		                                                                                      
			d_dd[rr+nr+nr] = 	d_uoy[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection) * nxpad * nzpad + d_Rjx[rr] + d_Rjz[rr] * nxpad ];
		
	}
}



/*
 * kernel lint3d_extract_gpu: 
 *		extract receiver data from the GPU and perform linear interpolation
 *
 * launch configuration:
 *		one thread per receiver
 *		Grid 	- 	(MIN(nr,ceil(nr/1024.0f)), 1, 1);
 *		Block 	-	(MIN(nr, 1024), 1, 1);
 *
 * Input:
 *		gpuID			- the gpu this kernel is executing on
 * 		nxpad			- number of points in the x direction of the d_u* arrays
 *		nzpad			- number of points in the z direction of the d_u* arrays
 *		nyinterior		- number of points in the y direction of the d_u* arrays
 *		nr				- number of receivers
 *		d_uo*			- displacements in the current time step for all grid points
 *		d_Rj*			- receiver coordinates
 *		d_Rw*			- linear interpolation coefficients
 *
 * Output:
 *      d_dd			- displacements at each receiver location (receiver data)
 * 							NOTE: madagascar expects components in data output to be 1=z, 2=x, 3=y
 *
 * Return:
 * 		void
 */
__global__ void lint3d_extract_gpu(int gpuID, float *d_dd, int nr, int nxpad, int nyinterior, int nzpad, float *d_uoz, float *d_uox, float *d_uoy, int *d_Rjz, int *d_Rjx, int *d_Rjy, float *d_Rw000, float *d_Rw001, float *d_Rw010, float *d_Rw011, float *d_Rw100, float *d_Rw101, float *d_Rw110, float *d_Rw111){

	
	int rr = threadIdx.x + blockIdx.x * blockDim.x;		// rr = the receiver this thread is extracting data for
	
	int haloCorrection = 0;	// GPU 0 does not have any additional halo cells in stress/acceleration arrays
	
	if (gpuID != 0){
		haloCorrection = 4;	// all other GPUs have additional halo cells in stress/acceleration arrays
	}
	
	if (rr < nr && 
		d_Rjy[rr] >= (gpuID * nyinterior) && 
		d_Rjy[rr] < (gpuID * nyinterior + nyinterior)){
			
			d_dd[rr] = 			d_uoz[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 0) * nxpad * nzpad + (d_Rjx[rr] + 0) + (d_Rjz[rr] + 0) * nxpad ]  * d_Rw000[rr] +
								d_uoz[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 0) * nxpad * nzpad + (d_Rjx[rr] + 0) + (d_Rjz[rr] + 1) * nxpad ]  * d_Rw001[rr] +
								d_uoz[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 0) * nxpad * nzpad + (d_Rjx[rr] + 1) + (d_Rjz[rr] + 0) * nxpad ]  * d_Rw010[rr] +
								d_uoz[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 0) * nxpad * nzpad + (d_Rjx[rr] + 1) + (d_Rjz[rr] + 1) * nxpad ]  * d_Rw011[rr] +
								d_uoz[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 1) * nxpad * nzpad + (d_Rjx[rr] + 0) + (d_Rjz[rr] + 0) * nxpad ]  * d_Rw100[rr] +
								d_uoz[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 1) * nxpad * nzpad + (d_Rjx[rr] + 0) + (d_Rjz[rr] + 1) * nxpad ]  * d_Rw101[rr] +
								d_uoz[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 1) * nxpad * nzpad + (d_Rjx[rr] + 1) + (d_Rjz[rr] + 0) * nxpad ]  * d_Rw110[rr] +
								d_uoz[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 1) * nxpad * nzpad + (d_Rjx[rr] + 1) + (d_Rjz[rr] + 1) * nxpad ]  * d_Rw111[rr];
				
			d_dd[rr+nr] = 		d_uox[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 0) * nxpad * nzpad + (d_Rjx[rr] + 0) + (d_Rjz[rr] + 0) * nxpad ]  * d_Rw000[rr] +
								d_uox[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 0) * nxpad * nzpad + (d_Rjx[rr] + 0) + (d_Rjz[rr] + 1) * nxpad ]  * d_Rw001[rr] +
								d_uox[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 0) * nxpad * nzpad + (d_Rjx[rr] + 1) + (d_Rjz[rr] + 0) * nxpad ]  * d_Rw010[rr] +
								d_uox[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 0) * nxpad * nzpad + (d_Rjx[rr] + 1) + (d_Rjz[rr] + 1) * nxpad ]  * d_Rw011[rr] +
								d_uox[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 1) * nxpad * nzpad + (d_Rjx[rr] + 0) + (d_Rjz[rr] + 0) * nxpad ]  * d_Rw100[rr] +
								d_uox[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 1) * nxpad * nzpad + (d_Rjx[rr] + 0) + (d_Rjz[rr] + 1) * nxpad ]  * d_Rw101[rr] +
								d_uox[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 1) * nxpad * nzpad + (d_Rjx[rr] + 1) + (d_Rjz[rr] + 0) * nxpad ]  * d_Rw110[rr] +
								d_uox[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 1) * nxpad * nzpad + (d_Rjx[rr] + 1) + (d_Rjz[rr] + 1) * nxpad ]  * d_Rw111[rr];
		
			d_dd[rr+nr+nr] = 	d_uoy[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 0) * nxpad * nzpad + (d_Rjx[rr] + 0) + (d_Rjz[rr] + 0) * nxpad ]  * d_Rw000[rr] +
								d_uoy[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 0) * nxpad * nzpad + (d_Rjx[rr] + 0) + (d_Rjz[rr] + 1) * nxpad ]  * d_Rw001[rr] +
								d_uoy[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 0) * nxpad * nzpad + (d_Rjx[rr] + 1) + (d_Rjz[rr] + 0) * nxpad ]  * d_Rw010[rr] +
								d_uoy[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 0) * nxpad * nzpad + (d_Rjx[rr] + 1) + (d_Rjz[rr] + 1) * nxpad ]  * d_Rw011[rr] +
								d_uoy[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 1) * nxpad * nzpad + (d_Rjx[rr] + 0) + (d_Rjz[rr] + 0) * nxpad ]  * d_Rw100[rr] +
								d_uoy[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 1) * nxpad * nzpad + (d_Rjx[rr] + 0) + (d_Rjz[rr] + 1) * nxpad ]  * d_Rw101[rr] +
								d_uoy[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 1) * nxpad * nzpad + (d_Rjx[rr] + 1) + (d_Rjz[rr] + 0) * nxpad ]  * d_Rw110[rr] +
								d_uoy[(d_Rjy[rr] - (gpuID * nyinterior) + haloCorrection + 1) * nxpad * nzpad + (d_Rjx[rr] + 1) + (d_Rjz[rr] + 1) * nxpad ]  * d_Rw111[rr];
	}
}


/* ------------------------------------ */
/* ABSORBING BOUNDARY CONDITION KERNELS */
/* ------------------------------------ */
/*
 * kernels abcone3d_apply: 
 *		Apply one-way wave equation BC
 *
 * Input:
 *		gpuID			- the gpu this kernel is executing on
 * 		nxpad			- number of points in the x direction of the d_u* arrays
 *		nzpad			- number of points in the z direction of the d_u* arrays
 *		nyinterior		- number of points in the y direction of the d_u* arrays
 *		d_uo			- displacements in the current time step for all grid points
 *		d_um			- displacements in the previous time step for all grid points
 * 		d_b*			- boundary condition coefficients as in (Clayton and Enquist, 1977)
 *
 * Output:
 *      d_uo			- displacements at each receiver location (receiver data)
 *
 * Return:
 * 		void
 */

// Apply one-way wave equation BC in the boundary region parallel to the XY-plane
__global__ void abcone3d_apply_XY(int gpuID, int nxpad, int nyinterior, int nzpad, float *d_uo, float *d_um, float *d_bzl, float *d_bzh){
	
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int y = threadIdx.y + blockIdx.y * blockDim.y;
	
	if (x < nxpad && y < nyinterior){
		
		if (gpuID != 0){
			y += 4;		// move y up by 4 to get out of the ghost cells in the d_u arrays
		}
		
		// low inxex
		if (blockIdx.z == 0){
				
				float bzl = d_bzl[(threadIdx.y + blockIdx.y * blockDim.y) * nxpad + x];	// this array is not padded with ghost cells
				
				float um5 = d_um[y * nxpad * nzpad + 5 * nxpad + x];
				float um4 = d_um[y * nxpad * nzpad + 4 * nxpad + x];
				float um3 = d_um[y * nxpad * nzpad + 3 * nxpad + x];
				float um2 = d_um[y * nxpad * nzpad + 2 * nxpad + x];
				float um1 = d_um[y * nxpad * nzpad + 1 * nxpad + x];
				
				float uo5 = d_uo[y * nxpad * nzpad + 5 * nxpad + x];
				
				float uo4 = um5 + (um4 - uo5) * bzl;

				float uo3 = um4 + (um3 - uo4) * bzl;

				float uo2 = um3 + (um2 - uo3) * bzl;

				float uo1 = um2 + (um1 - uo2) * bzl;


				d_uo[y * nxpad * nzpad + 4 * nxpad + x] = uo4;
                d_uo[y * nxpad * nzpad + 3 * nxpad + x] = uo3;
                d_uo[y * nxpad * nzpad + 2 * nxpad + x] = uo2;
                d_uo[y * nxpad * nzpad + 1 * nxpad + x] = uo1;

		}
		else {
				
				float bzh = d_bzh[(threadIdx.y + blockIdx.y * blockDim.y) * nxpad + x];
				
				float um6 = d_um[y * nxpad * nzpad + (nzpad-6) * nxpad + x];
				float um5 = d_um[y * nxpad * nzpad + (nzpad-5) * nxpad + x];
				float um4 = d_um[y * nxpad * nzpad + (nzpad-4) * nxpad + x];
				float um3 = d_um[y * nxpad * nzpad + (nzpad-3) * nxpad + x];
				float um2 = d_um[y * nxpad * nzpad + (nzpad-2) * nxpad + x];
				
				float uo6 = d_uo[y * nxpad * nzpad + (nzpad-6) * nxpad + x];
				
				float uo5 = um6 + (um5 - uo6) * bzh;
				
				float uo4 = um5 + (um4 - uo5) * bzh;
				
				float uo3 = um4 + (um3 - uo4) * bzh;
				
				float uo2 = um3 + (um2 - uo3) * bzh;
			
				d_uo[y * nxpad * nzpad + (nzpad-5) * nxpad + x] = uo5;
				d_uo[y * nxpad * nzpad + (nzpad-4) * nxpad + x] = uo4;
				d_uo[y * nxpad * nzpad + (nzpad-3) * nxpad + x] = uo3;
				d_uo[y * nxpad * nzpad + (nzpad-2) * nxpad + x] = uo2;
			
		}
		
	}
}

// Apply one-way wave equation BC in the boundary region parallel to the ZY-plane
__global__ void abcone3d_apply_ZY(int gpuID, int nxpad, int nyinterior, int nzpad, float *d_uo, float *d_um, float *d_bxl, float *d_bxh){
	
	int z = threadIdx.z + blockIdx.z * blockDim.z;
	int y = threadIdx.y + blockIdx.y * blockDim.y;
	
	if (z < nzpad && y < nyinterior){
		
		if (gpuID != 0){
			y += 4;		// move y up by 4 to get out of the ghost cells in the d_u arrays
		}
		
		if (blockIdx.x == 0){
			float um5 = d_um[y * nxpad * nzpad + z * nxpad + 5];
			float um4 = d_um[y * nxpad * nzpad + z * nxpad + 4];
			float um3 = d_um[y * nxpad * nzpad + z * nxpad + 3];
			float um2 = d_um[y * nxpad * nzpad + z * nxpad + 2];
			float um1 = d_um[y * nxpad * nzpad + z * nxpad + 1];
			
			float bxl = d_bxl[(threadIdx.y + blockIdx.y * blockDim.y) * nzpad + z];
			
			float uo5 = d_uo[y * nxpad * nzpad + z * nxpad + 5];
			
			float uo4 = um5 + (um4 - uo5) * bxl;
				
			float uo3 = um4 + (um3 - uo4) * bxl;
				
			float uo2 = um3 + (um2 - uo3) * bxl;
				
			float uo1 = um2 + (um1 - uo2) * bxl;
			
			d_uo[y * nxpad * nzpad + z * nxpad + 4] = uo4;
			d_uo[y * nxpad * nzpad + z * nxpad + 3] = uo3;
			d_uo[y * nxpad * nzpad + z * nxpad + 2] = uo2;
			d_uo[y * nxpad * nzpad + z * nxpad + 1] = uo1;	
		}
		else {
			float um6 = d_um[y * nxpad * nzpad + z * nxpad + (nxpad-6)];
			float um5 = d_um[y * nxpad * nzpad + z * nxpad + (nxpad-5)];
			float um4 = d_um[y * nxpad * nzpad + z * nxpad + (nxpad-4)];
			float um3 = d_um[y * nxpad * nzpad + z * nxpad + (nxpad-3)];
			float um2 = d_um[y * nxpad * nzpad + z * nxpad + (nxpad-2)];
			
			float bxh = d_bxh[(threadIdx.y + blockIdx.y * blockDim.y) * nzpad + z];
			
			float uo6 = d_uo[y * nxpad * nzpad + z * nxpad + (nxpad-6)];
			
			float uo5 = um6 + (um5 - uo6) * bxh;
				
			float uo4 = um5 + (um4 - uo5) * bxh;
				
			float uo3 = um4 + (um3 - uo4) * bxh;
				
			float uo2 = um3 + (um2 - uo3) * bxh;
			
			d_uo[y * nxpad * nzpad + z * nxpad + (nxpad-5)] = uo5;
			d_uo[y * nxpad * nzpad + z * nxpad + (nxpad-4)] = uo4;
			d_uo[y * nxpad * nzpad + z * nxpad + (nxpad-3)] = uo3;
			d_uo[y * nxpad * nzpad + z * nxpad + (nxpad-2)] = uo2;
		}
		
	}
}

// Apply one-way wave equation BC in the boundary region in the low Y-indicies parallel to the XZ-plane
__global__ void abcone3d_apply_XZ_low(int nxpad, int nzpad, float *d_uo, float *d_um, float *d_byl){
	
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int z = threadIdx.z + blockIdx.z * blockDim.z;
	
	if (x < nxpad && z < nzpad){
			
			float um5 = d_um[5 * nxpad * nzpad + z * nxpad + x];
			float um4 = d_um[4 * nxpad * nzpad + z * nxpad + x];
			float um3 = d_um[3 * nxpad * nzpad + z * nxpad + x];
			float um2 = d_um[2 * nxpad * nzpad + z * nxpad + x];
			float um1 = d_um[1 * nxpad * nzpad + z * nxpad + x];
			
			float byl = d_byl[x * nzpad + z];
			
			float uo5 = d_uo[5 * nxpad * nzpad + z * nxpad + x];
			
			float uo4 = um5 + (um4 - uo5) * byl;
			
			float uo3 = um4 + (um3 - uo4) * byl;
			
			float uo2 = um3 + (um2 - uo3) * byl;
			
			float uo1 = um2 + (um1 - uo2) * byl;
			
			d_uo[4 * nxpad * nzpad + z * nxpad + x] = uo4;
			d_uo[3 * nxpad * nzpad + z * nxpad + x] = uo3;
			d_uo[2 * nxpad * nzpad + z * nxpad + x] = uo2;
			d_uo[1 * nxpad * nzpad + z * nxpad + x] = uo1;
	}
	
}

// Apply one-way wave equation BC in the boundary region in the high Y-indicies parallel to the XZ-plane
__global__ void abcone3d_apply_XZ_high(int nxpad, int nylocal, int nzpad, float *d_uo, float *d_um, float *d_byh){
	
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int z = threadIdx.z + blockIdx.z * blockDim.z;
	
	if (x < nxpad && z < nzpad){

			float um6 = d_um[(nylocal-6) * nxpad * nzpad + z * nxpad + x];
			float um5 = d_um[(nylocal-5) * nxpad * nzpad + z * nxpad + x];
			float um4 = d_um[(nylocal-4) * nxpad * nzpad + z * nxpad + x];
			float um3 = d_um[(nylocal-3) * nxpad * nzpad + z * nxpad + x];
			float um2 = d_um[(nylocal-2) * nxpad * nzpad + z * nxpad + x];

			float byh = d_byh[x * nzpad + z];
			
			float uo6 = d_uo[(nylocal-6) * nxpad * nzpad + z * nxpad + x];
			
			float uo5 = um6 + (um5 - uo6) * byh;
			
			float uo4 = um5 + (um4 - uo5) * byh;
			
			float uo3 = um4 + (um3 - uo4) * byh;
			
			float uo2 = um3 + (um2 - uo3) * byh;
	
			d_uo[(nylocal-5) * nxpad * nzpad + z * nxpad + x] = uo5;
			d_uo[(nylocal-4) * nxpad * nzpad + z * nxpad + x] = uo4;
			d_uo[(nylocal-3) * nxpad * nzpad + z * nxpad + x] = uo3;
			d_uo[(nylocal-2) * nxpad * nzpad + z * nxpad + x] = uo2;
			
	}
}


/* ----------------------------------------------------- */
/* EXPONENTIAL-DAMPING SPONGE BOUNDARY CONDITION KERNELS */
/* ----------------------------------------------------- */
/*
 * kernels sponge3d_apply: 
 *		Apply the exponential-damping sponge BC
 *
 * Input:
 *		gpuID			- the gpu this kernel is executing on
 * 		nxpad			- number of points in the x direction of the d_u* arrays
 *		nzpad			- number of points in the z direction of the d_u* arrays
 *		nyinterior		- number of points in the y direction of the d_u* arrays
 *		d_uu			- displacements in the current time step for all grid points
 * 		spo				- damping coefficient
 *
 * Output:
 *      d_uu			- displacements at each receiver location (receiver data)
 *
 * Return:
 * 		void
 */
// Apply the exponential-damping sponge BC in the boundary region parallel to the XY-plane
__global__ void sponge3d_apply_XY(int gpuID, float *d_uu, int nxpad, int nyinterior, int nzpad, int nb){
	
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int y = threadIdx.y + blockIdx.y * blockDim.y;
	
	if (x < nxpad && y < nyinterior){
		
		if (gpuID != 0){
			y += 4;		// ghost cells in d_uu arrays need to be skipped over
		}
		
		for (int z = 0; z < nb; z++){
			// calculate the sponge coeff
			float spo = (nb - z - 1) / (sqrt(2.0) * 4.0f * nb);
			spo = exp(-spo * spo);

			d_uu[y * nxpad * nzpad + z * nxpad + x] *= spo;
			d_uu[y * nxpad * nzpad + (nzpad - z - 1) * nxpad + x] *= spo;
		}
		
	}
	
}

// Apply the exponential-damping sponge BC in the boundary region parallel to the ZY-plane
__global__ void sponge3d_apply_ZY(int gpuID, float *d_uu, int nxpad, int nyinterior, int nzpad, int nb, int nx, float spo){
	
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int z = threadIdx.y + blockIdx.y * blockDim.y;
	int y = threadIdx.z + blockIdx.z * blockDim.z;
	
	if (x < nb && y < nyinterior && z < nzpad){
		if (gpuID != 0){
			y += 4;
		}
		
		// calculate the sponge coeff
		spo = x / spo;
		spo = exp(-spo * spo);
		
		d_uu[y * nxpad * nzpad + z * nxpad + (nb - x - 1)] *= spo;
		d_uu[y * nxpad * nzpad + z * nxpad + (x + nb + nx)] *= spo;
		
	}
	
}

// Apply the exponential-damping sponge BC in the boundary region in the low Y-indicies parallel to the XZ-plane
__global__ void sponge3d_apply_XZ_low(float *d_uu, int nxpad, int nzpad, int nb){
	
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int z = threadIdx.z + blockIdx.z * blockDim.z;
	
	if (x < nxpad && z < nzpad){
		
		for (int y = 0; y < nb; y++){
			// calculate the sponge coeff
			float spo = (nb - y - 1) / (sqrt(2.0) * 4.0f * nb);
			spo = exp(-spo * spo);

			d_uu[y * nxpad * nzpad + z * nxpad + x] *= spo;
		}
		
	}
	
}

// Apply the exponential-damping sponge BC in the boundary region in the high Y-indicies parallel to the XZ-plane
__global__ void sponge3d_apply_XZ_high(float *d_uu, int nxpad, int nylocal, int nzpad, int nb){
	
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int z = threadIdx.z + blockIdx.z * blockDim.z;
	
	if (x < nxpad && z < nzpad){
		
		for (int y = 0; y < nb; y++){
			// calculate the sponge coeff
			float spo = (nb - y - 1) / (sqrt(2.0) * 4.0f * nb);
			spo = exp(-spo * spo);

			d_uu[(nylocal - y - 1) * nxpad * nzpad + z * nxpad + x] *= spo;
		}
		
	}
	
}

