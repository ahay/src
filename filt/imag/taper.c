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

void taper3(int nt3, int nt2, int nt1      /* taper lengths */, 
	    bool b3, bool b2, bool b1  /* taper in the beginning  */, 
	    int n3, int n2, int n1        /* dimensions */, 
	    float*** tt            /* [n3][n2][n1] tapered array (in and out) */)
/*< 3-D taper >*/
{
    int it,i1,i2,i3;
    float gain;

    for (it=0; it < nt3; it++) {
	gain = sinf(0.5*SF_PI*it/nt3);
	gain *= gain;
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		if (b3) tt[it][i2][i1] *= gain;
		tt[n3-it-1][i2][i1] *= gain;
	    }
	}
    }

    for (it=0; it < nt2; it++) {
	gain = sinf(0.5*SF_PI*it/nt2);
	gain *= gain;
	for (i3=0; i3 < n3; i3++) {
	    for (i1=0; i1 < n1; i1++) {
		if (b2) tt[i3][it][i1] *= gain;
		tt[i3][n2-it-1][i1] *= gain;
	    }
	}
    }

    for (it=0; it < nt1; it++) {
	gain = sinf(0.5*SF_PI*it/nt1);
	gain *= gain;
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		if (b1) tt[i3][i2][it] *= gain;
		tt[i3][i2][n1-it-1] *= gain;
	    }
	}
    }
}
