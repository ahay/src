/* Simple weight operator using CUDA */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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
#include "cuda_adjnull_kernels.cu"

static float* w;

void cuda_weight_init(float *w1)
/*< initialize >*/
{
    w = w1;
}

void cuda_weight_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy)
/*< linear operator >*/
{
    if (ny!=nx) sf_error("%s: size mismatch: %d != %d",__FILE__,ny,nx);

    cuda_adjnull (adj, add, nx, ny, xx, yy);
    if(adj) cuda_mul(true, xx, yy, w, nx);
    else cuda_mul(true, yy, xx, w, nx);
}


