/* Make wall weighting from window weighting. */
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

#include "mkwallwt.h"
#include "patch.h"

void mkwallwt(int dim     /* number of dimensions */, 
	      int* npatch /* number of patches [dim] */, 
	      int* nwall  /* data size [dim] */, 
	      int* nwind  /* patch size [dim] */, 
	      float* windwt /* window weighting (input) */, 
	      float* wallwt /* wall weighting (output) */)
/*< make wall weight >*/
{
    int i, j, ip, np, n, nw;

    np = 1; 
    n = 1;
    nw = 1;

    for (j=0; j < dim; j++) {
	np *= npatch[j];
	n *= nwall[j];
	nw *= nwind[j];
    }

    for (i = 0; i < n; i++) {
	wallwt[i] = 0.;
    }

    patch_init(dim, npatch, nwall, nwind);
  
    for (ip=0; ip < np; ip++) {
	patch_lop(true, true, n, nw, wallwt, windwt);
	patch_close ();
    }

    for (i = 0; i < n; i++) {
	if ( wallwt[i] != 0.) wallwt[i] = 1. / wallwt[i];
    }
}
