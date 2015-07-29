/* GPU Kernel Functions used in sfewefd2d_gpu */

/*
  Authors: Robin M. Weiss and Jeffrey Shragge

  This file contains the GPU kernel functions called in the ewefd2d_gpu module from
  the Madagascar software package (http://www.reproducibility.org).  The calling
  functions for these kernels can be found in the file Mewefd2d_gpu.cu.  For more 
  information, see (Weiss and Shragge, "Solving 3D Anisotropic Elastic Wave 
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
 * Copy array a with dimensions *_a into array b and expand data to dimensions *_b and
 *         also change the order of the axes in the array
 * NOTE: axes from slowest to fastest in array a are x,z.
 * NOTE: axes from slowest to fastest in array b are z,x.
 *
 * Input:
 *      a is an input matrix with dimensions x_a z_a
 *
 * Output:
 *      b is the output matrix with dimensions x_b z_b
 *
 * Return:
 * 		void
 */
void expand_cpu(float *a, float *b, int nb, int x_a, int x_b, int z_a, int z_b){
	
	
	// copy a into center of b, flipping the x and z axes
	for (int ix = 0; ix < x_a; ix++){
		for (int iz = 0; iz < z_a; iz++){
			b[(iz+nb) * x_b + (ix+nb)] = a[ix * z_a + iz];
		}
	}
	
	// expand z direction
	for (int ix = 0; ix < x_b; ix++){
		for (int iz = 0; iz < nb; iz++){
			b[iz * x_b + ix] = b[nb * x_b + ix];
			b[(z_b-iz-1) * x_b + ix] = b[(z_b-nb-1) * x_b + ix];
		}
	}
	
	// expand x direction
	for (int ix = 0; ix < nb; ix++){
		for (int iz = 0; iz < z_b; iz++){
			b[iz * x_b + ix] = b[iz * x_b + nb];
			b[iz * x_b + (x_b-ix-1)] = b[iz * x_b + (x_b-nb-1)];
		}
	}

}



/*
 * kernel computeRo: 
 *		compute 1/ro * dt^2 for all points in the model
 *
 * launch configuration:
 *		one thread per grid point
 *		Grid 	- (ceil(fdm->nxpad/16.0f), ceil(fdm->nzpad/16.0f))
 *		Block 	- (16,16)
 *
 * Input:
 *      d_ro		- density for all points in the model
 *		dt			- time step increment
 * 		nxpad		- number of points in the x direction of the d_ro array
 *		nzpad		- number of points in the z direction of the d_ro array
 *
 * Output:
 *      d_ro		- contains 1/ro * dt^2 for all points in the model
 *
 * Return:
 * 		void
 */
__global__ void computeRo(float *d_ro, float dt, int nxpad, int nzpad, int nop){
	
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int z = threadIdx.y + blockIdx.y * blockDim.y;
	
	if (x < nxpad && z < nzpad){	// ensure this thread falls inside the model dimensions
		
		int addr = z * nxpad + x;
	
		d_ro[addr] = dt * dt / d_ro[addr];
		
	}
}


/*
 * kernel dispToStrain_strainToStress: 
 *		apply the FD stencil to the displacement arrays to calculate strain in the model
 *      then multiply by stiffness coefficients to obtain stresses in the model
 *
 * launch configuration:
 *		one thread per grid point
 *		Grid 	- (ceil(fdm->nxpad/32.0f), ceil(fdm->nzpad/32.0f))
 *		Block 	- (32,32)
 *
 * Input:
 * 		nxpad			- number of points in the x direction of the d_t*, d_uo*, and d_c* arrays
 *		nzpad			- number of points in the z direction of the d_t*, d_uo*, and d_c* arrays
 *      d_uo*			- displacements for x and z componenets for all points in the model
 *      d_c*			- stiffness coefficients for all points in the model
 *		idx, idz		- reciprocal of the spatial sampling in the x and z directions
 *
 * Output:
 *      d_t*			- the stress values in resptective components for all points in the model
 *
 * Return:
 * 		void
 */
__global__ void dispToStrain_strainToStress(float *d_txx, float *d_tzz, float *d_tzx, float *d_uox, float *d_uoz, float *d_c11, float *d_c33, float *d_c55, float *d_c13, float idx, float idz, int nxpad, int nzpad, int nop){
	
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int z = threadIdx.y + blockIdx.y * blockDim.y;
	
	
	if (x < nxpad && z < nzpad){	// ensure this thread falls inside the model dimensions
		int i;
		int addr = z * nxpad + x;	// the address of the point this thread is responsible for

		float txx, tzz, tzx;
		
		// store stiffness coeffients in registers
		float c11 = d_c11[addr];
		float c13 = d_c13[addr];
		float c33 = d_c33[addr];
		float c55 = d_c55[addr];
		
		// if this thread is atleast nop points away from the edge of the model, apply the FD stencil
		if (x >= nop && x < nxpad-nop && z >= nop && z < nzpad-nop){
			
			txx = 0;
			tzz = 0;
			float tzx1 = 0;
			float tzx2 = 0;
			for (i = 4; i > 0; i--){
				txx  += C[i-1]*(d_uox[z * nxpad + (x+i)] - d_uox[z * nxpad + (x-i)]);
				tzz  += C[i-1]*(d_uoz[(z+i) * nxpad + x] - d_uoz[(z-i) * nxpad + x]);
				tzx1 += C[i-1]*(d_uoz[z * nxpad + (x+i)] - d_uoz[z * nxpad + (x-i)]);
				tzx2 += C[i-1]*(d_uox[(z+i) * nxpad + x] - d_uox[(z-i) * nxpad + x]);
			}
			
			txx *= idx;
			tzz *= idz;
			tzx = tzx1 * idx + tzx2 * idz;
		}
		else {
			txx = d_txx[addr];
			tzz = d_tzz[addr];
			tzx = d_tzx[addr];
		}
		
		// multiple the resulting strain value by stiffness coefficient to obtain stresses for this point
		d_txx[addr] = c11 * txx + 
					  c13 * tzz;
			
		d_tzz[addr] = c13 * txx + 
					  c33 * tzz; 
	
		d_tzx[addr] = c55 * tzx;
		
	}
	
}


/*
 * kernel freeSurf: 
 *		set all stress tensors containing a z-component to 0 in the z=0 to z=nb region of the domain
 *
 * launch configuration:
 *		one thread per grid point
 *		Grid 	- 	(ceil(fdm->nxpad/16.0f), ceil(fdm->nb/16.0f))
 *		Block 	-	(16,16)
 *
 * Input:
 * 		nxpad			- the number of points in the x direction of the d_t* arrays
 *		nzpad			- the number of points in the z direction of the d_t* arrays
 *		nb				- amount of padding added to the model in the low side of the XY-plane
 *		d_t*			- stress values involving the z-component for all points in the model
 *
 * Output:
 *      d_t*			- stress values involving the z-component for all points in the model
 *
 * Return:
 * 		void
 */
__global__ void freeSurf(float *d_tzz, float *d_tzx, int nxpad, int nb){
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int z = threadIdx.y + blockIdx.y * blockDim.y;
	
	if (x < nxpad && z < nb){	// ensure this thread falls within the free area
		int addr = z * nxpad + x;
		
		d_tzz[addr] = 0;
		d_tzx[addr] = 0;
	
	}
}



/*
 * kernel lint2d_bell_gpu: 
 *		inject stress or accel source from the d_ww array into the array d_uu at a x,z coordinate given by d_jx[] and d_jz[]
 *			using interpolation coefficients in the d_Sw arrays
 *
 * launch configuration:
 *		2d square of threads with side length 2*nbell+1 and 1 block per source
 *		Grid 	- 	(ns, 1, 1);
 *		Block 	-	(2 * nbell + 1, 2 * nbell + 1, 1);
 *
 * Input:
 *		d_uu			- acceleration/stress array to inject source wavelet into
 * 		nxpad			- number of points in the x direction of the d_uu array
 *		it				- current iteration number
 *		nc				- number of vector components
 *		ns				- number of sources
 *		c				- component of source wavelet being injected
 *		d_ww			- the source wavelet(s)
 *		d_j*			- x- and z-coordinates of each source ia
 *		d_Sw*			- linear interpolation coefficients for each source ia
 *
 * Output:
 *      d_uu			- accleration/stress array with source wavelet injected
 *
 * Return:
 * 		void
 */
__global__ void lint2d_bell_gpu(float *d_uu, float *d_ww, float *d_Sw00, float *d_Sw01, float *d_Sw10, float *d_Sw11, float *d_bell, int *d_jz, int *d_jx, int it, int nc, int ns, int c, int nbell, int nxpad){

	int ix = threadIdx.x; // 0 to 2*nbell + 1
	int iz = threadIdx.y; // 0 to 2*nbell + 1
	int ia = blockIdx.x;
	
	float wa = d_ww[it * nc * ns + c * ns + ia] * d_bell[(iz * (2*nbell+1)) + ix];
	
	atomicAdd(&d_uu[((d_jz[ia] - nbell) + iz    ) * nxpad + ((d_jx[ia] - nbell) + ix    )], (-(wa * d_Sw00[ia])));
	atomicAdd(&d_uu[((d_jz[ia] - nbell) + iz + 1) * nxpad + ((d_jx[ia] - nbell) + ix    )], (-(wa * d_Sw01[ia])));
	atomicAdd(&d_uu[((d_jz[ia] - nbell) + iz    ) * nxpad + ((d_jx[ia] - nbell) + ix + 1)], (-(wa * d_Sw10[ia])));
	atomicAdd(&d_uu[((d_jz[ia] - nbell) + iz + 1) * nxpad + ((d_jx[ia] - nbell) + ix + 1)], (-(wa * d_Sw11[ia])));
	
}



/*
 * kernel stressToAccel: 
 *		apply FD stencil to stress arrays (d_t*) to obtain acceleration values for grid points
 *
 * launch configuration:
 *		one thread per grid point
 *		Grid 	- 	(ceil((fdm->nxpad-(2*NOP))/32.0f), ceil((fdm->nzpad-(2*NOP))/32.0f))
 *		Block 	-	(32,32)
 *
 * Input:
 * 		nxpad			- number of points in the x direction of the d_t* and d_ua* arrays
 *		nzpad			- number of points in the z direction of the d_t* and d_ua* arrays
 *		idx, idz		- reciprocal of the spatial sampling in the x and z directions respectively
 *		d_t*			- stresses for all grid points in the model for all vector components
 *
 * Output:
 *      d_ua*			- accelerations at all grid points in the model
 *
 * Return:
 * 		void
 */
__global__ void stressToAcceleration(float *d_uax, float *d_uaz, float* d_txx, float *d_tzz, float *d_tzx, float idx, float idz, int nxpad, int nzpad){
	
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int z = threadIdx.y + blockIdx.y * blockDim.y;
	
	if (x < nxpad - 8 && z < nzpad - 8){ 	// ensure this thread falls inside the domain
		int offset = (nxpad * 4) + 4 + x;
		
		// apply FD stencil to d_txx and d_tzx
		float uax = (C[3]*(d_txx[offset + 4 + (z * nxpad)] - d_txx[offset - 4 + (z * nxpad)]) +
					 C[2]*(d_txx[offset + 3 + (z * nxpad)] - d_txx[offset - 3 + (z * nxpad)]) +
					 C[1]*(d_txx[offset + 2 + (z * nxpad)] - d_txx[offset - 2 + (z * nxpad)]) +
					 C[0]*(d_txx[offset + 1 + (z * nxpad)] - d_txx[offset - 1 + (z * nxpad)]))*idx;
		uax += 		(C[3]*(d_tzx[offset + ((z + 4) * nxpad)] - d_tzx[offset + ((z - 4) * nxpad)]) +
					 C[2]*(d_tzx[offset + ((z + 3) * nxpad)] - d_tzx[offset + ((z - 3) * nxpad)]) +
					 C[1]*(d_tzx[offset + ((z + 2) * nxpad)] - d_tzx[offset + ((z - 2) * nxpad)]) +
					 C[0]*(d_tzx[offset + ((z + 1) * nxpad)] - d_tzx[offset + ((z - 1) * nxpad)]))*idz;
		
		// apply FD stencil to d_tzx and d_tzz
		float uaz = (C[3]*(d_tzx[offset + 4 + (z * nxpad)] - d_tzx[offset - 4 + (z * nxpad)]) +
					 C[2]*(d_tzx[offset + 3 + (z * nxpad)] - d_tzx[offset - 3 + (z * nxpad)]) +
					 C[1]*(d_tzx[offset + 2 + (z * nxpad)] - d_tzx[offset - 2 + (z * nxpad)]) +
				  	 C[0]*(d_tzx[offset + 1 + (z * nxpad)] - d_tzx[offset - 1 + (z * nxpad)]))*idx;
		uaz += 		(C[3]*(d_tzz[offset + ((z + 4) * nxpad)] - d_tzz[offset + ((z - 4) * nxpad)]) +
					 C[2]*(d_tzz[offset + ((z + 3) * nxpad)] - d_tzz[offset + ((z - 3) * nxpad)]) +
					 C[1]*(d_tzz[offset + ((z + 2) * nxpad)] - d_tzz[offset + ((z - 2) * nxpad)]) +
					 C[0]*(d_tzz[offset + ((z + 1) * nxpad)] - d_tzz[offset + ((z - 1) * nxpad)]))*idz;
		
		// store calculated acceleration in d_ua arrays
		d_uaz[offset + (z * nxpad)] = uaz;
		d_uax[offset + (z * nxpad)] = uax;
		
	}
}


/*
 * kernel stepTime: 
 *		calculate displacements for all grid points in the next time step
 *
 * launch configuration:
 *		one thread per grid point
 *		Grid 	- 		(ceil(fdm->nxpad/16.0f), ceil(fdm->nzpad/12.0f))
 *		Block	-		(16,12)
 *
 * Input:
 * 		nxpad			- number of points in the x direction of the d_u* arrays
 *		nzpad			- number of points in the z direction of the d_u* arrays
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
__global__ void stepTime(float *d_upz, float *d_uoz, float *d_umz, float *d_uaz, float *d_upx, float *d_uox, float *d_umx, float *d_uax, float *d_ro, int nxpad, int nzpad){
	
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int z = threadIdx.y + blockIdx.y * blockDim.y;
	
	if (x < nxpad && z < nzpad){	// ensure thread falls inside the domain
		
		int addr = z * nxpad + x;
		
		float ro = d_ro[addr];
	
		float uoz = d_uoz[addr];
		float umz = d_umz[addr];
		float uaz = d_uaz[addr];
		
		float uox = d_uox[addr];
		float umx = d_umx[addr];
		float uax = d_uax[addr];
		
		d_upz[addr] = 2.0f * uoz - umz + uaz * ro;
		d_upx[addr] = 2.0f * uox - umx + uax * ro;

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
 * 		nxpad			- number of points in the x direction of the d_u* arrays
 *		nzpad			- number of points in the z direction of the d_u* arrays
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
__global__ void abcone2d_apply_LR_gpu(float *d_uo, float *d_um, float *d_bxl, float *d_bxh, int nxpad, int nzpad){
	
	int iz = threadIdx.y + blockIdx.y * blockDim.y;
	
	if (iz < nzpad){
		if (blockIdx.x == 0){	// DO LEFT BOUNDARY
				
			float bxl = d_bxl[iz];
			
			float um5 = d_um[iz * nxpad + 5];
			float um4 = d_um[iz * nxpad + 4];
			float um3 = d_um[iz * nxpad + 3];
			float um2 = d_um[iz * nxpad + 2];
			float um1 = d_um[iz * nxpad + 1];
			
			float uo5 = d_uo[iz * nxpad + 5];
			
			float uo4 = um5 + (um4 - uo5) * bxl;
			
			float uo3 = um4 + (um3 - uo4) * bxl;
			
			float uo2 = um3 + (um2 - uo3) * bxl;
			
			float uo1 = um2 + (um1 - uo2) * bxl;
			
			d_uo[iz * nxpad + 4] = uo4;
			d_uo[iz * nxpad + 3] = uo3;
			d_uo[iz * nxpad + 2] = uo2;
			d_uo[iz * nxpad + 1] = uo1;
				
		}
		else {	// DO RIGHT BOUNDARY
				
			float bxh = d_bxh[iz];
			
			float um6 = d_um[iz * nxpad + nxpad-6];
			float um5 = d_um[iz * nxpad + nxpad-5];
			float um4 = d_um[iz * nxpad + nxpad-4];
			float um3 = d_um[iz * nxpad + nxpad-3];
			float um2 = d_um[iz * nxpad + nxpad-2];
			
			float uo6 = d_uo[iz * nxpad + nxpad-6];
			
			float uo5 = um6 + (um5 - uo6) * bxh;
			                                                                                                                          
			float uo4 = um5 + (um4 - uo5) * bxh;
			                                                                                                                         
			float uo3 = um4 + (um3 - uo4) * bxh;
			                                                                                                                           
			float uo2 = um3 + (um2 - uo3) * bxh;
			
			
			d_uo[iz * nxpad + nxpad-5] = uo5;
			d_uo[iz * nxpad + nxpad-4] = uo4;
			d_uo[iz * nxpad + nxpad-3] = uo3;
			d_uo[iz * nxpad + nxpad-2] = uo2;
				
		}
	}
}


__global__ void abcone2d_apply_TB_gpu(float *d_uo, float *d_um, float *d_bzl, float *d_bzh, int nxpad, int nzpad, int fsrf){
	
	int ix = threadIdx.x + blockIdx.x * blockDim.x;
	
	if (ix < nxpad){
		if (blockIdx.y == 0 && fsrf == 0){	// DO TOP BOUNDARY
				
			float bzl = d_bzl[ix];
			
			float um5 = d_um[5 * nxpad + ix];
			float um4 = d_um[4 * nxpad + ix];
			float um3 = d_um[3 * nxpad + ix];
			float um2 = d_um[2 * nxpad + ix];
			float um1 = d_um[1 * nxpad + ix];
			
			float uo5 = d_uo[5 * nxpad + ix];
			
			float uo4 = um5 + (um4 - uo5) * bzl;
			                                                                                                          
			float uo3 = um4 + (um3 - uo4) * bzl;
			                                                                                                         
			float uo2 = um3 + (um2 - uo3) * bzl;
			                                                                                                         
			float uo1 = um2 + (um1 - uo2) * bzl;
			
			
			d_uo[4 * nxpad + ix] = uo4;
			d_uo[3 * nxpad + ix] = uo3;
			d_uo[2 * nxpad + ix] = uo2;
			d_uo[1 * nxpad + ix] = uo1;
				
		}
		else if (blockIdx.y == 1){	// DO BOTTOM BOUNDARY
			
			float bzh = d_bzh[ix];
			
			float um6 = d_um[(nzpad-6) * nxpad + ix];
			float um5 = d_um[(nzpad-5) * nxpad + ix];
			float um4 = d_um[(nzpad-4) * nxpad + ix];
			float um3 = d_um[(nzpad-3) * nxpad + ix];
			float um2 = d_um[(nzpad-2) * nxpad + ix];
			
			float uo6 = d_uo[(nzpad-6) * nxpad + ix];
			
			float uo5 = um6 + (um5 - uo6) * bzh;
			                                                                                                                                 
			float uo4 = um5 + (um4 - uo5) * bzh;
			                                                                                                                                
			float uo3 = um4 + (um3 - uo4) * bzh;
			                                                                                                                                 
			float uo2 = um3 + (um2 - uo3) * bzh;
			
			
			d_uo[(nzpad-5) * nxpad + ix] = uo5;
			d_uo[(nzpad-4) * nxpad + ix] = uo4;
			d_uo[(nzpad-3) * nxpad + ix] = uo3;
			d_uo[(nzpad-2) * nxpad + ix] = uo2;
		}
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
 * 		nx			- number of points in the x direction of the model domain
 *		nz			- number of points in the z direction of the model domain
 *		nb			- amount of padding added to each dimensions of the model
 *		d_spo		- damping sponge coefficients
 *		d_uu		- displacements in the current time step for all grid points
 *
 * Output:
 *      d_uu			- displacements at each receiver location (receiver data)
 *
 * Return:
 * 		void
 */
__global__ void sponge2d_apply_TB_gpu(float *d_uu, float *d_spo, int nxpad, int nb, int nz){
	
	int ix = threadIdx.x + blockIdx.x * blockDim.x;
	int iz = blockIdx.y;
	
	// do top and bottom
	if (ix < nxpad){
		float w = d_spo[nb - iz - 1];
		d_uu[iz * nxpad + ix] *= w;
		d_uu[((nb - iz - 1) + nb + nz) * nxpad + ix] *= w;
	}
	
}


__global__ void sponge2d_apply_LR_gpu(float *d_uu, float *d_spo, int nxpad, int nb, int nx){
	
	int ix = threadIdx.x + blockIdx.x * blockDim.x;
	int iz = blockIdx.y;
	
	// do left and right
	if (ix < nb){
		float w = d_spo[nb - ix - 1];
		d_uu[iz * nxpad + ix] *= w;
		d_uu[iz * nxpad + (((nb - 1) - ix) + nb + nx)] *= w;
	}
	
}



/*
 * kernel lint2d_extract_gpu: 
 *		extract receiver data from the GPU and perform linear interpolation
 *
 * launch configuration:
 *		one thread per receiver
 *		Grid 	- 	(MIN(nr,ceil(nr/1024.0f)), 1, 1);
 *		Block 	-	(MIN(nr, 1024), 1, 1);
 *
 * Input:
 * 		nxpad			- number of points in the x direction of the d_u* arrays
 *		nzpad			- number of points in the z direction of the d_u* arrays
 *		nr				- number of receivers
 *		d_uo*			- displacements in the current time step for all grid points
 *		d_Rj*			- receiver coordinates
 *		d_Rw*			- linear interpolation coefficients
 *
 * Output:
 *      d_dd			- displacements at each receiver location (receiver data)
 * 							NOTE: madagascar expects components in data output to be 1=z, 2=x
 *
 * Return:
 * 		void
 */
__global__ void lint2d_extract_gpu(float *d_dd, int nr, int nxpad, float *d_uoz, float* d_uox, int *d_Rjz, int *d_Rjx, float *d_Rw00, float *d_Rw01, float *d_Rw10, float *d_Rw11){
	
	
	int rr = threadIdx.x + blockIdx.x * blockDim.x;
	
	if (rr < nr){
		d_dd[rr] = d_uoz[d_Rjz[rr] * nxpad + (d_Rjx[rr])] * d_Rw00[rr] +
		    	 d_uoz[(d_Rjz[rr]+1) * nxpad + d_Rjx[rr]] * d_Rw01[rr] +
		    	 d_uoz[d_Rjz[rr] * nxpad + (d_Rjx[rr]+1)] * d_Rw10[rr] +
		    	 d_uoz[(d_Rjz[rr]+1) * nxpad + (d_Rjx[rr]+1)] * d_Rw11[rr];
		
		d_dd[rr + nr] = d_uox[d_Rjz[rr] * nxpad + (d_Rjx[rr])] * d_Rw00[rr] +
		    	 	  d_uox[(d_Rjz[rr]+1) * nxpad + d_Rjx[rr]] * d_Rw01[rr] +
		    	 	  d_uox[ d_Rjz[rr] * nxpad + (d_Rjx[rr]+1)] * d_Rw10[rr] +
		    	 	  d_uox[(d_Rjz[rr]+1) * nxpad + (d_Rjx[rr]+1)] * d_Rw11[rr];
	}
	
	
}

