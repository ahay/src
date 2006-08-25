/* Patch operation (in core) */
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

#include <rsf.h>
/*^*/

#include "patch.h"

static int dim, ipatch;
static int *npatch, *nwind, *nwall, ii[SF_MAX_DIM], jj[SF_MAX_DIM]; 

void patch_init(int dim_in     /* number of dimensions */, 
		int* npatch_in /* number of patches [dim] */, 
		int* nwall_in  /* data size [dim] */, 
		int* nwind_in  /* patch size [dim] */)
/*< initialize >*/
{
    dim = dim_in;
    npatch = npatch_in; 
    nwall = nwall_in; 
    nwind = nwind_in; 
    ipatch = 0;
}

void patch_lop (bool adj, bool add, int nx, int ny, 
		float* wall, float* wind)
/*< patch operator >*/
{
    int i, j, shift;
 
    sf_adjnull (adj, add, nx, ny, wall, wind);
    sf_line2cart(dim, npatch, ipatch, jj);  
    for(i = 0; i < dim; i++) {
	if(npatch[i] == 1) {
	    jj[i] = 0;
	} else if (jj[i] == npatch[i]-1) {
	    jj[i] = nwall[i] - nwind[i];
	} else {	    
	    jj[i] = jj[i]*(nwall[i] - nwind[i])/(npatch[i] - 1.0);
	}
    }

    /* shift till the patch start */
    shift = sf_cart2line(dim, nwall, jj); 
    for(i = 0; i < ny; i++) {
	sf_line2cart(dim, nwind, i, ii); 
	j = sf_cart2line(dim, nwall, ii) + shift;   
	if (adj) wall[j] += wind[i];
	else     wind[i] += wall[j];
    }
}

void patch_close(void)
/*< Move to next patch. >*/
{
    ipatch++; /* increase patch counter */
}
