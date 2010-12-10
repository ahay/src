/* Convolution in patches with a local filter */
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
#include "helix_icai.h"

#include "loconvol_internal.h"

static sf_filter aa;

void loconvol_internal_init(sf_filter aa_in)
/*< initialize with the first filter >*/
{
    aa = aa_in;
}

void loconvol_internal_lop(bool adj, bool add, int nx, int ny,
		  float *xx, float *yy)
/*< convolve >*/
{
	helix_icai_init(aa,1);
    aa++;

    helix_icai_lop(adj, add, nx, ny, xx, yy);
}
