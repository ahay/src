/* Parcel operator for testing patching */
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

#include "parcel.h"
#include "patch.h"

static int np, nw;

void parcel_init(int dim /* number of dimensions */, 
		 int *npatch /* number of patches [dim] */, 
		 int *nwall  /* data size [dim] */, 
		 int *nwind  /* patch size [dim] */)
/*< initialize >*/
{
    int i;

    patch_init (dim, npatch, nwall, nwind);
    np = 1;
    nw = 1;  
    for (i=0; i < dim; i++) {
	np *= npatch[i]; /* compute number of patches */
	nw *= nwind[i];  /* compute window size */
    }
}

void parcel_lop(bool adj, bool add, int n, int mw, float* wall, float* wind)
/*< parcel operator >*/
{
    int ip;

    if (mw != np*nw) sf_error("%s: wrong dimensions",__FILE__);
  
    sf_adjnull (adj, add, n, mw, wall, wind);
    for (ip=0; ip < np; ip++) {
	patch_lop (adj, true, n, nw, wall, wind+ip*nw);
	patch_close ();
    }
}
