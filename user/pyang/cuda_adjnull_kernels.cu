/* Claerbout-style adjoint zeroing using CUDA
Note: #include "cuda_adjnull_kernels.cu" in main function if needed */
/*
  Copyright (C) 2014 Xi'an Jiaotong University, Pengliang Yang
  
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
#include "cuda_def.cu"

void cuda_adjnull (bool adj /* adjoint flag */, 
		 bool add /* addition flag */, 
		 int nx   /* size of x */, 
		 int ny   /* size of y */, 
		 float* x, 
		 float* y) 
/*< Zeros out the output (unless add is true). 
  Useful first step for any linear operator. >*/
{
    if(add) return;
    
    if(adj) 	cudaMemset(x, 0, nx*sizeof(float));
    else	cudaMemset(y, 0, ny*sizeof(float));
}
