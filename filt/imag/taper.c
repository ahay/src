/* Tapering */
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
#include <math.h>
#include <rsf.h>

#include "taper.h"

void taper2(int ntx, int nty      /* taper lengths */, 
	    bool begx, bool begy  /* taper in the beginning  */, 
	    int nx, int ny        /* dimensions */, 
	    float** tt            /* [nx][ny] tapered array (in and out) */)
/*< 2-D taper >*/
{
    int it,ix,iy;
    float gain;

    for (it=0; it < ntx; it++) {
	gain = sinf(0.5*SF_PI*it/ntx);
	gain *= gain;
	for (iy=0; iy < ny; iy++) {
	    if (begx) tt[it][iy] *= gain;
	    tt[nx-it-1][iy] *= gain;
	}
    }

    for (it=0; it < nty; it++) {
	gain = sinf(0.5*SF_PI*it/nty);
	gain *= gain;
	for (ix=0; ix < nx; ix++) {
	    if (begy) tt[ix][it] *= gain;
	    tt[ix][ny-it-1] *= gain;
	}
    }
}
